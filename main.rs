// SPDX-License-Identifier: GPL-3.0-or-later
// Rust implementation of bwakit/bwa-postalt.js (r985).
// Original algorithm: Heng Li / BWA. See COPYING and README.md.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread;

const BIN_SHIFT: u32 = 13;

type Res<T> = Result<T, Box<dyn std::error::Error>>;

fn parse_int(s: &str) -> i64 { s.parse::<i64>().unwrap_or(0) }



fn walk_cigar(cigar: &str, mut action: impl FnMut(u8, i64)) {
    let b = cigar.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if !b[i].is_ascii_digit() { i += 1; continue; }
        let start = i;
        while i < b.len() && b[i].is_ascii_digit() { i += 1; }
        if i == b.len() { break; }
        let op = b[i];
        if matches!(op, b'M' | b'I' | b'D' | b'S' | b'H' | b'N') {
            action(op, parse_int(&cigar[start..i]));
        }
        i += 1;
    }
}

#[derive(Clone, Debug)]
struct CigarOp { op: u8, len: i64 }

#[derive(Clone)]
struct AltMap {
    // half-open interval on ALT query
    start: i64,
    end: i64,
    query_len: i64,
    primary: String,
    rev: bool,
    // In original JS: pos = SAM POS - 1; stored offset = pos - 1.
    offset: i64,
    cigar: Vec<CigarOp>,
}

#[derive(Clone)]
struct PriMap { start: i64, end: i64, alt: String }

trait Bounds { fn bounds(&self) -> (i64, i64); }
impl Bounds for AltMap { fn bounds(&self) -> (i64, i64) { (self.start, self.end) } }
impl Bounds for PriMap { fn bounds(&self) -> (i64, i64) { (self.start, self.end) } }



struct IntervalIndex<T> {
    entries: Vec<T>,
    first: Vec<Option<usize>>,
    preceding: Vec<Option<usize>>,
    max_bin: i64,
}
impl<T: Bounds> IntervalIndex<T> {
    fn new(mut entries: Vec<T>) -> Self {
        entries.sort_by_key(|m| m.bounds().0);
        let max_bin = entries.iter().map(|m| (m.bounds().1 - 1) >> BIN_SHIFT).max().unwrap_or(-1);
        let mut first = vec![None; (max_bin + 1).max(0) as usize];
        for (i, m) in entries.iter().enumerate() {
            let (b, e) = m.bounds();
            if e <= 0 { continue; }
            let a = (b >> BIN_SHIFT).max(0) as usize;
            let z = ((e - 1) >> BIN_SHIFT) as usize;
            for slot in first.iter_mut().take(z.saturating_add(1)).skip(a) {
                if slot.is_none() { *slot = Some(i); }
            }
        }
        let mut preceding = Vec::with_capacity(first.len());
        let mut previous = None;
        for &v in &first {
            if v.is_some() { previous = v; }
            preceding.push(previous);
        }
        Self { entries, first, preceding, max_bin }
    }
    fn overlaps(&self, b: i64, e: i64, mut fun: impl FnMut(&T)) {
        if self.entries.is_empty() || b > e || (b >> BIN_SHIFT) > self.max_bin { return; }
        let x = (b >> BIN_SHIFT).max(0) as usize;
        let from = self.first.get(x).copied().flatten().or_else(|| {
            // Match the upstream fallback (search backward from the end bin - 1).
            let z = ((e - 1) >> BIN_SHIFT) - 1;
            if z < 0 { None } else {
                self.preceding.get((z as usize).min(self.preceding.len().saturating_sub(1)))
                    .copied().flatten()
            }
        }).unwrap_or(0);
        for m in self.entries.iter().skip(from) {
            let (s, t) = m.bounds();
            if s >= e { break; }
            if t > b { fun(m); }
        }
    }
}

struct Reference {
    alt: HashMap<String, IntervalIndex<AltMap>>,
    pri: HashMap<String, IntervalIndex<PriMap>>,
    alt_names: HashSet<String>,
    hla_genes: HashSet<String>,
    hla_chr: Option<String>,
}

fn hla_gene(s: &str) -> Option<&str> {
    if !s.starts_with("HLA-") { return None; }
    let star = s.find('*')?;
    if star < 5 || !s[star+1..].bytes().next()?.is_ascii_digit() { return None; }
    if s[..star].bytes().any(|c| c.is_ascii_whitespace()) { return None; }
    Some(&s[..star])
}

fn load_reference(path: &Path) -> Res<Reference> {
    let f = BufReader::new(File::open(path)?);
    let mut alts: HashMap<String, Vec<AltMap>> = HashMap::new();
    let mut pris: HashMap<String, Vec<PriMap>> = HashMap::new();
    let mut alt_names = HashSet::new();
    let mut hla_genes = HashSet::new();
    let mut hla_chr = None;
    for line in f.lines() {
        let line = line?;
        if line.starts_with('@') { continue; }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 11 { continue; }
        alt_names.insert(t[0].to_string());
        let pos = parse_int(t[3]) - 1;
        let flag = parse_int(t[1]);
        if flag & 4 != 0 || t[2] == "*" { continue; }
        if let Some(g) = hla_gene(t[0]) {
            hla_genes.insert(g.to_owned());
            hla_chr = Some(t[2].to_owned());
        }
        let mut cigar: Vec<CigarOp> = Vec::new();
        let mut l_qaln = 0;
        let mut l_tlen = 0;
        let mut l_qclip = 0;
        walk_cigar(t[5], |op, len| {
            cigar.push(CigarOp { op: if op == b'H' { b'S' } else { op }, len });
            match op {
                b'M' => { l_qaln += len; l_tlen += len; }
                b'I' => l_qaln += len,
                b'S' | b'H' => l_qclip += len,
                b'D' | b'N' => l_tlen += len,
                _ => {}
            }
        });
        if cigar.is_empty() { continue; }
        let rev = flag & 16 != 0;
        let j = if rev { cigar.len() - 1 } else { 0 };
        let start = if cigar[j].op == b'S' { cigar[j].len } else { 0 };
        alts.entry(t[0].to_string()).or_default().push(AltMap {
            start, end: start + l_qaln, query_len: l_qaln + l_qclip,
            primary: t[2].to_owned(), rev, offset: pos - 1, cigar,
        });
        pris.entry(t[2].to_owned()).or_default().push(PriMap {
            start: pos, end: pos + l_tlen, alt: t[0].to_owned(),
        });
    }
    Ok(Reference {
        alt: alts.into_iter().map(|(k,v)| (k, IntervalIndex::new(v))).collect(),
        pri: pris.into_iter().map(|(k,v)| (k, IntervalIndex::new(v))).collect(),
        alt_names, hla_genes, hla_chr,
    })
}

fn cigar_to_pos(cigar: &[CigarOp], pos: i64) -> i64 {
    let (mut x, mut y) = (0, 0);
    for c in cigar {
        match c.op {
            b'M' => {
                if y <= pos && pos < y + c.len { return x + (pos - y); }
                x += c.len;
                y += c.len;
            }
            b'D' => x += c.len,
            b'I' => {
                if y <= pos && pos < y + c.len { return x; }
                y += c.len;
            }
            b'S' | b'H' => {
                if y <= pos && pos < y + c.len { return -1; }
                y += c.len;
            }
            _ => {} // upstream cigar2pos ignores N
        }
    }
    -1
}

#[derive(Clone)]
struct Lift { ctg: String, rev: bool, start: i64, end: i64 }

struct Hit {
    ctg: String,
    start: i64,
    end: i64,
    rev: bool,
    cigar: String,
    nm: i64,
    hard: bool,
    score: i64,
    pctg: String,
    pstart: i64,
    pend: i64,
    original: usize,
    group: usize,
    lifted: Vec<Lift>,
}

fn parse_hit(ctg: &str, signed_pos: &str, cigar: &str, nm: &str) -> Hit {
    let start = parse_int(signed_pos.get(1..).unwrap_or("0")) - 1;
    let rev = signed_pos.starts_with('-');
    let mut matches = 0;
    let mut ins = 0;
    let mut dels = 0;
    let mut skipped = 0;
    let mut n_ins = 0;
    let mut n_del = 0;
    let mut hard = false;
    walk_cigar(cigar, |op,len| match op {
        b'M' => matches += len,
        b'D' => { n_del += 1; dels += len; }
        b'I' => { n_ins += 1; ins += len; }
        b'N' => skipped += len,
        b'H' => hard = true,
        _ => {}
    });
    let nm = parse_int(nm).max(ins + dels);
    let score = matches - 5 * (nm - dels - ins) - 6 * (n_del + n_ins) - (dels + ins);
    let end = start + matches + dels + skipped;
    Hit {
        ctg: ctg.to_owned(), start, end, rev, cigar: cigar.to_owned(), nm, hard, score,
        pctg: ctg.to_owned(), pstart: start, pend: end, original: 0, group: 0, lifted: Vec::new(),
    }
}



struct Record {
    raw: String,
    fields: Option<Vec<String>>,
    flag: i64,
}
impl Record {
    fn new(raw: String) -> Self {
        let flag = parse_int(raw.split('\t').nth(1).unwrap_or("0"));
        Self { raw, fields: None, flag }
    }
    fn from_fields(fields: Vec<String>) -> Self {
        let flag = parse_int(fields.get(1).map(String::as_str).unwrap_or("0"));
        Self { raw: String::new(), fields: Some(fields), flag }
    }
    fn field(&self, index: usize) -> &str {
        if let Some(ref fields) = self.fields {
            fields.get(index).map(String::as_str).unwrap_or("")
        } else {
            self.raw.split('\t').nth(index).unwrap_or("")
        }
    }
    fn name(&self) -> &str { self.field(0) }
    fn tag(&self, prefix: &str) -> Option<&str> {
        if let Some(ref fields) = self.fields {
            fields.iter().skip(11).find_map(|f| f.strip_prefix(prefix))
        } else {
            self.raw.split('\t').skip(11).find_map(|f| f.strip_prefix(prefix))
        }
    }
    fn ensure_fields(&mut self) -> &mut Vec<String> {
        if self.fields.is_none() {
            self.fields = Some(self.raw.split('\t').map(str::to_owned).collect());
            self.raw.clear();
        }
        self.fields.as_mut().unwrap()
    }
    fn set_mapq(&mut self, mapq: i64) {
        if self.field(4) != mapq.to_string() {
            let v = self.ensure_fields();
            if v.len() > 4 { v[4] = mapq.to_string(); }
        }
    }
    fn add_tag(&mut self, tag: String) { self.ensure_fields().push(tag); }
    fn print(&self, out: &mut impl Write) -> io::Result<()> {
        if let Some(ref fields) = self.fields {
            if let Some(first) = fields.first() {
                out.write_all(first.as_bytes())?;
                for f in fields.iter().skip(1) {
                    out.write_all(b"\t")?;
                    out.write_all(f.as_bytes())?;
                }
            }
        } else { out.write_all(self.raw.as_bytes())?; }
        out.write_all(b"\n")
    }
}

fn collect_hla(idx: &HashMap<String, IntervalIndex<PriMap>>,
               ctg: &str, start: i64, end: i64, hla: &mut HashSet<String>) {
    if let Some(index) = idx.get(ctg) {
        index.overlaps(start, end, |m| {
            if let Some(g) = hla_gene(&m.alt) { hla.insert(g.to_owned()); }
        });
    }
}

struct PostAlt<'a> { reference: &'a Reference, min_pa_ratio: f64 }
impl PostAlt<'_> {
    fn process(&self, mut rec: Record, buffered: &mut Vec<Record>, hla: &mut HashSet<String>) {
        let flag = rec.flag;
        if flag & 4 != 0 { buffered.push(rec); return; }
        // Avoid field allocation and CIGAR scoring for the overwhelmingly common
        // read that cannot participate in ALT processing or HLA extraction.
        // HLA-chromosome reads still take the general path for HLA hit tracking.
        if !self.reference.alt_names.contains(rec.field(2))
            && rec.tag("XA:Z:").is_none()
            && self.reference.hla_chr.as_deref() != Some(rec.field(2)) {
            buffered.push(rec);
            return;
        }
        let rname = rec.field(2).to_owned();
        let orig_mapq = parse_int(rec.field(4));
        let pos = rec.field(3).to_owned();
        let cigar = rec.field(5).to_owned();
        let nm = rec.tag("NM:i:").unwrap_or("0").to_owned();
        let signed = format!("{}{}", if flag & 16 != 0 {'-'} else {'+'}, pos);
        let initial = parse_hit(&rname, &signed, &cigar, &nm);
        if self.reference.hla_chr.as_deref() == Some(rname.as_str()) {
            collect_hla(&self.reference.pri, &rname, initial.start, initial.end, hla);
        }
        if initial.hard { buffered.push(rec); return; }
        let mut hits = vec![initial];
        if let Some(xa) = rec.tag("XA:Z:") {
            let xa = xa.split_whitespace().next().unwrap_or("");
            for entry in xa.split(';').filter(|s| !s.is_empty()) {
                let mut a = entry.split(',');
                if let (Some(ctg), Some(spos), Some(cigar), Some(nm)) =
                    (a.next(), a.next(), a.next(), a.next()) {
                    hits.push(parse_hit(ctg, spos, cigar, nm));
                }
            }
        }
        if !hits.iter().any(|h| self.reference.alt_names.contains(&h.ctg)) {
            buffered.push(rec); return;
        }

        // Lift all ALT hits. Upstream grouping deliberately uses only the FIRST
        // lift when an ALT has several possible ALT-to-primary mappings.
        let mut n_rpt_lifted = 0;
        let mut rpt_lifted: Option<Lift> = None;
        for (i, h) in hits.iter_mut().enumerate() {
            if let Some(idx) = self.reference.alt.get(&h.ctg) {
                idx.overlaps(h.start, h.end, |a| {
                    let (s,e) = if !a.rev {
                        (cigar_to_pos(&a.cigar, h.start), cigar_to_pos(&a.cigar, h.end - 1) + 1)
                    } else {
                        (cigar_to_pos(&a.cigar, a.query_len - h.end),
                         cigar_to_pos(&a.cigar, a.query_len - h.start - 1) + 1)
                    };
                    if s < 0 || e < 0 { return; }
                    h.lifted.push(Lift { ctg: a.primary.clone(), rev: h.rev != a.rev,
                                         start: s + a.offset, end: e + a.offset });
                    if i == 0 { n_rpt_lifted += 1; }
                });
                if i == 0 && n_rpt_lifted == 1 { rpt_lifted = h.lifted.first().cloned(); }
            }
        }
        for (i,h) in hits.iter_mut().enumerate() {
            h.original = i;
            if let Some(l) = h.lifted.first() {
                h.pctg = l.ctg.clone(); h.pstart = l.start; h.pend = l.end;
            }
        }
        hits.sort_by(|a,b| a.pctg.cmp(&b.pctg).then(a.pstart.cmp(&b.pstart)));
        let mut groups = 0usize;
        let mut last_chr: Option<String> = None;
        let mut end = 0;
        for h in &mut hits {
            if last_chr.as_deref() != Some(h.pctg.as_str()) { groups += 1; end = 0; }
            else if h.pstart >= end { groups += 1; }
            h.group = groups - 1;
            last_chr = Some(h.pctg.clone());
            end = end.max(h.pend);
        }
        let reported_i = hits.iter().position(|h| h.original == 0).unwrap_or(0);
        let reported_g = hits[reported_i].group;
        let n_group0 = hits.iter().filter(|h| h.group == reported_g).count();

        let mut mapq = orig_mapq;
        if n_group0 > 1 {
            let mut group_max = vec![(i64::MIN, 0usize); groups];
            for h in &hits {
                if group_max[h.group].0 < h.score {
                    group_max[h.group] = (h.score, h.group);
                }
            }
            // Stable sort, exactly like JS Array.sort with comparator y[0]-x[0].
            group_max.sort_by(|a,b| b.0.cmp(&a.0));
            mapq = if group_max[0].1 == reported_g {
                if group_max.len() == 1 { 60 } else { 6 * (group_max[0].0 - group_max[1].0) }
            } else { 0 };
            mapq = mapq.min(60);
            mapq = if !self.reference.alt.contains_key(&rname) {
                mapq.min(orig_mapq)
            } else {
                mapq.max(orig_mapq)
            };
        }

        if self.reference.hla_chr.as_deref() == Some(hits[reported_i].pctg.as_str()) {
            let mut start = 1i64 << 30;
            let mut stop = 0;
            for h in hits.iter().filter(|h| h.group == reported_g) {
                start = start.min(h.pstart);
                stop = stop.max(h.pend);
            }
            if let Some(ref chr) = self.reference.hla_chr {
                collect_hla(&self.reference.pri, chr, start, stop, hla);
            }
        }

        if n_rpt_lifted <= 1 {
            for s in buffered.iter_mut() {
                let is_overlap = if let Some(ref l) = rpt_lifted {
                    if s.field(2) != l.ctg || ((s.flag & 16 != 0) != l.rev) { false }
                    else {
                        let b = parse_int(s.field(3)) - 1;
                        let mut e = b;
                        // Deliberately use the CURRENT record's CIGAR (t[5]),
                        // as the upstream JS does, rather than the buffered s[5].
                        walk_cigar(&cigar, |op,len| {
                            if matches!(op, b'M'|b'D'|b'N') { e += len; }
                        });
                        b < l.end && l.start < e
                    }
                } else { false };
                let om = s.tag("om:i:").map(parse_int).unwrap_or(-1);
                let pa = s.tag("pa:f:").and_then(|f| f.parse::<f64>().ok()).unwrap_or(10.0);
                if is_overlap {
                    if om > 0 { s.set_mapq(om); }
                    let adjusted_mapq = parse_int(s.field(4)).min(mapq);
                    s.set_mapq(adjusted_mapq);
                } else if pa < self.min_pa_ratio {
                    if om < 0 {
                        let old_mapq = s.field(4).to_string();
                        s.add_tag(format!("om:i:{}", old_mapq));
                    }
                    s.set_mapq(0);
                }
            }
        }

        let mut lifted_strs = Vec::with_capacity(hits.len());
        for h in &hits {
            let mut lifted = String::new();
            for l in &h.lifted {
                use std::fmt::Write as _;
                let _ = write!(lifted, "{},{},{},{};", l.ctg, l.start, l.end,
                               if l.rev { '-' } else { '+' });
            }
            lifted_strs.push(lifted);
        }
        let rg = rec.tag("RG:Z:").map(|s| format!("RG:Z:{}", s));
        rec.set_mapq(mapq);
        if n_group0 > 1 { rec.add_tag(format!("om:i:{}", orig_mapq)); }
        if !lifted_strs[reported_i].is_empty() {
            rec.add_tag(format!("lt:Z:{}", lifted_strs[reported_i]));
        }

        let qname = rec.field(0).to_owned();
        let rnext = rec.field(6).to_owned();
        let pnext = rec.field(7).to_owned();
        let tlen = rec.field(8).to_owned();
        let seq = rec.field(9).to_owned();
        let qual = rec.field(10).to_owned();
        let reported_rev = hits[reported_i].rev;
        buffered.push(rec);
        let mut revseq: Option<String> = None;
        let mut revqual: Option<String> = None;
        for (i,h) in hits.iter().enumerate() {
            if i == reported_i || h.group != reported_g || !self.reference.alt.contains_key(&h.ctg) {
                continue;
            }
            let same_strand = h.rev == reported_rev;
            if !same_strand && revseq.is_none() {
                revseq = Some(reverse_complement(&seq));
                revqual = Some(qual.chars().rev().collect());
            }
            let mut fields = vec![
                qname.clone(),
                (if same_strand { flag | 0x800 } else { (flag ^ 0x10) | 0x800 }).to_string(),
                h.ctg.clone(), (h.start+1).to_string(), mapq.to_string(), h.cigar.clone(),
                if rnext == "=" && h.ctg != rname { rname.clone() } else { rnext.clone() },
                pnext.clone(), tlen.clone(),
                if same_strand { seq.clone() } else { revseq.as_ref().unwrap().clone() },
                if same_strand { qual.clone() } else { revqual.as_ref().unwrap().clone() },
            ];
            fields.push(format!("NM:i:{}", h.nm));
            if !lifted_strs[i].is_empty() { fields.push(format!("lt:Z:{}", lifted_strs[i])); }
            if let Some(ref rg) = rg { fields.push(rg.clone()); }
            buffered.push(Record::from_fields(fields));
        }
    }
}

fn reverse_complement(seq: &str) -> String {
    fn rc(c: u8) -> u8 {
        match c {
            b'W'|b'w' => c, b'S'|b's' => c,
            b'A'=>b'T', b'a'=>b't', b'T'=>b'A', b't'=>b'a', b'U'=>b'A', b'u'=>b'a',
            b'G'=>b'C', b'g'=>b'c', b'C'=>b'G', b'c'=>b'g',
            b'Y'=>b'R', b'y'=>b'r', b'R'=>b'Y', b'r'=>b'y',
            b'K'=>b'M', b'k'=>b'm', b'M'=>b'K', b'm'=>b'k',
            b'B'=>b'V', b'b'=>b'v', b'V'=>b'B', b'v'=>b'b',
            b'D'=>b'H', b'd'=>b'h', b'H'=>b'D', b'h'=>b'd',
            b'N'=>b'N', b'n'=>b'n', _=>c,
        }
    }
    seq.as_bytes().iter().rev().map(|&b| rc(b) as char).collect()
}

struct Group {
    recs: Vec<Record>,
    hla: HashSet<String>,
}
impl Group {
    fn new() -> Self { Self { recs: Vec::with_capacity(8), hla: HashSet::new() } }
    fn flush(&mut self, output: &mut impl Write,
             files: &mut Option<HashMap<String, BufWriter<File>>>) -> io::Result<()> {
        for rec in &self.recs { rec.print(output)?; }
        if !self.recs.is_empty() && !self.hla.is_empty() {
            if let Some(files) = files {
                let first = &self.recs[0];
                let name = format!("{}/{}{}", first.name(), (first.flag >> 6) & 3,
                                   if first.flag & 16 != 0 { '-' } else { '+' });
                let seq = first.field(9);
                let qual = first.field(10);
                for gene in &self.hla {
                    if let Some(out) = files.get_mut(gene) {
                        writeln!(out, "@{}\n{}\n+\n{}", name, seq, qual)?;
                    }
                }
            }
        }
        self.recs.clear();
        self.hla.clear();
        Ok(())
    }
}

struct Args { alt: String, sam: Option<String>, hla_prefix: Option<String>, ratio: f64, threads: usize }
fn usage() {
    eprintln!("Usage: postalt-rs [-t THREADS] [-p HLA_FASTQ_PREFIX] [-r MIN_PA_RATIO] <alt.sam> [aln.sam]\n\
               Or:    bwa-mem3 mem ... | postalt-rs ref.fa.alt | samtools sort -n ...\n\
               -p prefix: write per-HLA FASTQs (forces serial mode for exact compatibility)\n\
               -t INT: processing workers (default 1; use 4 for parallel processing)\n\
               -r FLOAT: reduce MAPQ for non-overlapping hits when pa < ratio (default 1)\n\
               -v: print upstream-compatible version marker r985");
}
fn parse_args() -> Res<Args> {
    let mut args = env::args().skip(1).peekable();
    let mut hla_prefix = None;
    let mut ratio: f64 = 1.0;
    let mut threads: usize = 1;
    let mut pos = Vec::new();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-v" | "--version" => { println!("r985"); std::process::exit(0); }
            "-h" | "--help" => { usage(); std::process::exit(0); }
            "-p" => hla_prefix = Some(args.next().ok_or("Missing argument to -p")?),
            "-r" => ratio = args.next().ok_or("Missing argument to -r")?.parse()?,
            "-t" | "--threads" => threads = args.next().ok_or("Missing argument to -t")?.parse()?,
            "--" => { pos.extend(args); break; }
            _ if arg.starts_with("-p") && arg.len()>2 => hla_prefix=Some(arg[2..].to_owned()),
            _ if arg.starts_with("-r") && arg.len()>2 => ratio=arg[2..].parse()?,
            _ if arg.starts_with("-t") && arg.len()>2 => threads=arg[2..].parse()?,
            _ if arg.starts_with("--threads=") => threads=arg[10..].parse()?,
            _ => pos.push(arg),
        }
    }
    if pos.is_empty() || pos.len() > 2 {
        usage(); return Err("Expected <alt.sam> [aln.sam]".into());
    }
    if threads == 0 || threads > 64 { return Err("-t/--threads must be between 1 and 64".into()); }
    Ok(Args { alt: pos[0].clone(), sam: pos.get(1).cloned(), hla_prefix, ratio: ratio.min(1.0), threads })
}



fn process_serial(
    processor: &PostAlt<'_>,
    input: impl BufRead,
    output: &mut impl Write,
    hla_files: &mut Option<HashMap<String, BufWriter<File>>>,
) -> io::Result<()> {
    let mut group = Group::new();
    for line in input.lines() {
        let mut line = line?;
        if line.ends_with('\r') { line.pop(); }
        if line.starts_with('@') {
            output.write_all(line.as_bytes())?;
            output.write_all(b"\n")?;
            continue;
        }
        if line.is_empty() { continue; }
        let rec = Record::new(line);
        if let Some(first) = group.recs.first() {
            if first.name() != rec.name() || (first.flag & 0xc0) != (rec.flag & 0xc0) {
                group.flush(output, hla_files)?;
            }
        }
        processor.process(rec, &mut group.recs, &mut group.hla);
    }
    group.flush(output, hla_files)?;
    output.flush()
}



enum BatchPart {
    Header(String),
    Group(Vec<Record>),
}
struct WorkBatch { id: usize, parts: Vec<BatchPart> }

fn compute_batch(processor: &PostAlt<'_>, parts: Vec<BatchPart>) -> io::Result<Vec<u8>> {
    let mut output = Vec::with_capacity(1 << 20);
    let mut hla_files = None; // Parallel mode is disabled if -p was requested.
    for part in parts {
        match part {
            BatchPart::Header(line) => {
                output.extend_from_slice(line.as_bytes());
                output.push(b'\n');
            }
            BatchPart::Group(reads) => {
                let mut group = Group::new();
                for rec in reads {
                    processor.process(rec, &mut group.recs, &mut group.hla);
                }
                group.flush(&mut output, &mut hla_files)?;
            }
        }
    }
    Ok(output)
}

const BATCH_GROUP_LIMIT: usize = 256;
const BATCH_BYTE_LIMIT: usize = 1 << 20;

struct BatchBuilder {
    parts: Vec<BatchPart>,
    groups: usize,
    bytes: usize,
}
impl BatchBuilder {
    fn new() -> Self { Self { parts: Vec::new(), groups: 0, bytes: 0 } }
    fn push_group(&mut self, group: Vec<Record>) {
        self.bytes += group.iter().map(|r| r.raw.len()).sum::<usize>();
        self.groups += 1;
        self.parts.push(BatchPart::Group(group));
    }
    fn push_header(&mut self, line: String) {
        self.bytes += line.len();
        self.parts.push(BatchPart::Header(line));
    }
    fn ready(&self, max_groups: usize, max_bytes: usize) -> bool {
        self.groups >= max_groups || self.bytes >= max_bytes
    }
    fn take(&mut self, id: usize) -> WorkBatch {
        let parts = std::mem::take(&mut self.parts);
        self.groups = 0;
        self.bytes = 0;
        WorkBatch { id, parts }
    }
}

fn channel_closed() -> io::Error {
    io::Error::new(io::ErrorKind::BrokenPipe, "parallel postalt worker/writer disconnected")
}




fn dispatch(
    batch: WorkBatch,
    jobs: &[SyncSender<WorkBatch>],
    acks: &Receiver<()>,
    outstanding: &mut usize,
    max_in_flight: usize,
) -> io::Result<()> {
    while *outstanding >= max_in_flight {
        acks.recv().map_err(|_| channel_closed())?;
        *outstanding -= 1;
    }
    let worker_id = batch.id % jobs.len();
    jobs[worker_id].send(batch).map_err(|_| channel_closed())?;
    *outstanding += 1;
    Ok(())
}




fn process_parallel<R: BufRead, W: Write + Send>(
    processor: &PostAlt<'_>,
    input: R,
    output: W,
    workers: usize,
    max_groups: usize,
    max_bytes: usize,
) -> Res<W> {
    assert!(workers > 1 && max_groups > 0 && max_bytes > 0);
    thread::scope(|scope| -> Res<W> {
        let (result_tx, result_rx) = mpsc::sync_channel::<(usize, io::Result<Vec<u8>>)>(workers * 2);
        let (ack_tx, ack_rx) = mpsc::channel::<()>();
        let writer = scope.spawn(move || -> io::Result<(W, usize)> {
            let mut output = output;
            let mut pending = BTreeMap::<usize, Vec<u8>>::new();
            let mut next_id = 0usize;
            for (id, result) in result_rx {
                pending.insert(id, result?);
                while let Some(bytes) = pending.remove(&next_id) {
                    output.write_all(&bytes)?;
                    next_id += 1;
                    let _ = ack_tx.send(());
                }
            }
            if !pending.is_empty() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "missing parallel batch"));
            }
            output.flush()?;
            Ok((output, next_id))
        });

        // Dedicated bounded job channel for each worker avoids a mutex around
        // the receiver. Round-robin scheduling ensures finite reordering.
        let mut senders = Vec::with_capacity(workers);
        let mut worker_handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (tx, rx) = mpsc::sync_channel::<WorkBatch>(2);
            let result_tx = result_tx.clone();
            let handle = scope.spawn(move || {
                while let Ok(job) = rx.recv() {
                    let result = compute_batch(processor, job.parts);
                    if result_tx.send((job.id, result)).is_err() { break; }
                }
            });
            senders.push(tx);
            worker_handles.push(handle);
        }
        drop(result_tx);

        // Keep input grouping exact: the same QNAME and FLAG 0xc0 segment must
        // never be split across jobs. Headers are forwarded in original order.
        let read_result = (|| -> io::Result<usize> {
            let mut builder = BatchBuilder::new();
            let mut current: Vec<Record> = Vec::new();
            let mut sent = 0usize;
            let mut outstanding = 0usize;
            let max_in_flight = workers * 3;
            for line in input.lines() {
                let mut line = line?;
                if line.ends_with('\r') { line.pop(); }
                if line.starts_with('@') {
                    // Valid SAM has all headers before alignment records.
                    if !current.is_empty() {
                        builder.push_group(std::mem::take(&mut current));
                    }
                    builder.push_header(line);
                } else if !line.is_empty() {
                    let rec = Record::new(line);
                    if let Some(first) = current.first() {
                        if first.name() != rec.name() || (first.flag & 0xc0) != (rec.flag & 0xc0) {
                            builder.push_group(std::mem::take(&mut current));
                        }
                    }
                    current.push(rec);
                }
                if builder.ready(max_groups, max_bytes) {
                    dispatch(builder.take(sent), &senders, &ack_rx,
                             &mut outstanding, max_in_flight)?;
                    sent += 1;
                }
            }
            if !current.is_empty() { builder.push_group(current); }
            if !builder.parts.is_empty() {
                dispatch(builder.take(sent), &senders, &ack_rx,
                         &mut outstanding, max_in_flight)?;
                sent += 1;
            }
            Ok(sent)
        })();
        drop(senders);
        let mut worker_panicked = false;
        for handle in worker_handles {
            if handle.join().is_err() { worker_panicked = true; }
        }
        let writer_result = writer.join().map_err(|_| "parallel output writer panicked")?;
        let (output, completed) = writer_result?;
        if worker_panicked { return Err("postalt worker panicked".into()); }
        let submitted = read_result?;
        if submitted != completed { return Err("parallel output lost one or more batches".into()); }
        Ok(output)
    })
}

fn run() -> Res<()> {
    let opt = parse_args()?;
    let reference = load_reference(Path::new(&opt.alt))?;
    let processor = PostAlt { reference: &reference, min_pa_ratio: opt.ratio };
    let mut hla_files = if let Some(ref prefix) = opt.hla_prefix {
        let mut files = HashMap::new();
        for gene in &reference.hla_genes {
            let filename = format!("{}.{}.fq", prefix, gene);
            files.insert(gene.clone(), BufWriter::new(File::create(filename)?));
        }
        Some(files)
    } else { None };
    let input: Box<dyn BufRead> = if let Some(path) = opt.sam {
        if path == "-" { Box::new(BufReader::with_capacity(1 << 20, io::stdin())) }
        else { Box::new(BufReader::with_capacity(1 << 20, File::open(path)?)) }
    } else { Box::new(BufReader::with_capacity(1 << 20, io::stdin())) };
    if opt.threads > 1 && opt.hla_prefix.is_none() {
        let output = BufWriter::with_capacity(1 << 20, io::stdout());
        let _ = process_parallel(&processor, input, output, opt.threads,
                                 BATCH_GROUP_LIMIT, BATCH_BYTE_LIMIT)?;
    } else {
        if opt.threads > 1 && opt.hla_prefix.is_some() {
            eprintln!("postalt-rs: -p HLA extraction uses serial mode to preserve exact output");
        }
        let stdout = io::stdout();
        let mut output = BufWriter::with_capacity(1 << 20, stdout.lock());
        process_serial(&processor, input, &mut output, &mut hla_files)?;
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("postalt-rs: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cigar_translation() {
        let cigar = vec![CigarOp{op:b'M',len:10}, CigarOp{op:b'I',len:3},
                         CigarOp{op:b'M',len:10}, CigarOp{op:b'D',len:4}, CigarOp{op:b'M',len:5}];
        assert_eq!(cigar_to_pos(&cigar, 2), 2);
        assert_eq!(cigar_to_pos(&cigar, 11), 10);
        assert_eq!(cigar_to_pos(&cigar, 13), 10);
        assert_eq!(cigar_to_pos(&cigar, 14), 11);
        assert_eq!(cigar_to_pos(&cigar, 30), -1);
    }
    #[test]
    fn hit_scores() {
        let h = parse_hit("ALT", "+11", "10M2I20M3D5M", "7");
        assert_eq!(h.start, 10);
        assert_eq!(h.end, 48);
        assert_eq!(h.score, 35 - 5*(7-5) - 6*2 - 5);
    }
    #[test]
    fn rc_iupac() {
        assert_eq!(reverse_complement("ACGTNryk"), "mryNACGT");
    }
    #[test]
    fn lifted_alt_and_reverse_supplementary() {
        let p = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mini.alt.sam"));
        let reference = load_reference(p).unwrap();
        let handler = PostAlt { reference: &reference, min_pa_ratio: 1.0 };
        let mut buf = Vec::new();
        let mut hla = HashSet::new();
        let q = "rB\t0\tALT1\t12\t20\t20M\t*\t0\t0\tACGTTGCAAAACCCCCGGGG\tIIIIIIIIIIIIIIIIIIII\tNM:i:0\tXA:Z:ALT1,-13,20M,0;";
        handler.process(Record::new(q.to_owned()), &mut buf, &mut hla);
        assert_eq!(buf.len(), 2);
        assert_eq!(buf[0].field(4), "60");
        assert_eq!(buf[0].tag("om:i:"), Some("20"));
        assert_eq!(buf[0].tag("lt:Z:"), Some("chr1,110,130,+;"));
        assert_eq!(buf[1].field(1), "2064");
        assert_eq!(buf[1].field(3), "13");
        assert_eq!(buf[1].tag("lt:Z:"), Some("chr1,111,131,-;"));
        assert_eq!(buf[1].field(9), reverse_complement("ACGTTGCAAAACCCCCGGGG").as_str());
    }
    #[test]
    fn reverse_alt_lifts_to_primary() {
        let p = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mini.alt.sam"));
        let reference = load_reference(p).unwrap();
        let handler = PostAlt { reference: &reference, min_pa_ratio: 1.0 };
        let mut buf = Vec::new();
        let mut hla = HashSet::new();
        let q = "rF\t0\tALT2\t12\t35\t20M\t*\t0\t0\tACGTTGCAAAACCCCCGGGG\tIIIIIIIIIIIIIIIIIIII\tNM:i:0";
        handler.process(Record::new(q.to_owned()), &mut buf, &mut hla);
        assert_eq!(buf.len(), 1);
        assert_eq!(buf[0].field(4), "35");
        assert_eq!(buf[0].tag("lt:Z:"), Some("chr2,268,288,-;"));
    }
    #[test]
    fn index_overlap() {
        let idx = IntervalIndex::new(vec![
            PriMap{start:0,end:100,alt:"a".into()},
            PriMap{start:10000,end:10100,alt:"b".into()}]);
        let mut v = Vec::new();
        idx.overlaps(50,80,|m| v.push(m.alt.clone()));
        assert_eq!(v, vec!["a"]);
        v.clear(); idx.overlaps(9000,10050,|m| v.push(m.alt.clone()));
        assert_eq!(v, vec!["b"]);
    }
}

// Regression tests for reverse-complement behaviour expected from k8 v1.2.
#[cfg(test)]
mod rc_regression_tests {
    use super::reverse_complement;

    #[test]
    fn canonical_bases_and_unknown_symbols() {
        assert_eq!(reverse_complement("ACGTN"), "NACGT");
        assert_eq!(reverse_complement("A=CG.TN"), "NA.CG=T");
    }

    #[test]
    fn iupac_bases_and_case() {
        assert_eq!(reverse_complement("RYKMBDHVWS"), "SWBDHVKMRY");
        assert_eq!(reverse_complement("acgt"), "acgt");
    }
}




#[cfg(test)]
mod parallel_regression_tests {
    use super::*;
    use std::io::Cursor;

    fn run_both(sam: &[u8], workers: usize, max_groups: usize, max_bytes: usize) {
        let reference = load_reference(Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mini.alt.sam"
        ))).unwrap();
        let postalt = PostAlt { reference: &reference, min_pa_ratio: 1.0 };
        let mut serial = Vec::<u8>::new();
        process_serial(&postalt, Cursor::new(sam), &mut serial, &mut None).unwrap();
        let parallel = process_parallel(
            &postalt, Cursor::new(sam), Vec::<u8>::new(),
            workers, max_groups, max_bytes,
        ).unwrap();
        assert_eq!(serial, parallel,
                   "parallel mismatch: workers={}, groups={}, bytes={}",
                   workers, max_groups, max_bytes);
    }

    #[test]
    fn fixture_byte_identical_at_many_batch_boundaries() {
        let sam = include_bytes!("../tests/fixtures/mini.sam");
        for workers in [2, 3, 4] {
            for groups in [1, 2, 5, 256] {
                run_both(sam, workers, groups, 512);
            }
        }
    }

    #[test]
    fn large_stream_crosses_many_batches_and_retains_header_order() {
        let sam = include_str!("../tests/fixtures/mini.sam");
        let mut lines = sam.lines();
        let mut big = String::new();
        let mut records = Vec::new();
        for l in lines.by_ref() {
            if l.starts_with('@') { big.push_str(l); big.push('\n'); }
            else { records.push(l.to_owned()); break; }
        }
        records.extend(lines.map(str::to_owned));
        for n in 0..320 {
            for row in &records {
                // Ensure distinct QNAMEs across repetitions, retaining the
                // same name for corresponding primary + ALT records.
                let (qname, rest) = row.split_once('\t').unwrap();
                big.push_str(&format!("{}_{}\t{}\n", n, qname, rest));
            }
        }
        for workers in [2, 4] {
            for groups in [1, 7, 128] {
                run_both(big.as_bytes(), workers, groups, 4000);
            }
        }
    }

    #[test]
    fn empty_sam_and_headers_only_do_not_deadlock() {
        run_both(b"", 4, 1, 1);
        run_both(b"@HD\tVN:1.6\n@SQ\tSN:ALT1\tLN:100\n", 4, 1, 1);
    }
}
