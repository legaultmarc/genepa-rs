#![allow(dead_code)]

use std::fmt;
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::io::{BufReader, BufRead};
use std::fs::File;


#[derive(Debug)]
pub struct VarFieldIdx {
  // Denotes field indices in a delimited file containing variants. 
  pub delimiter: char,
  pub name: usize,
  pub chrom: usize,
  pub pos: usize,
  pub a1: usize,
  pub a2: usize,
}


pub struct DelimitedVariantsReader {
    iter: Box<std::io::Lines<std::io::BufReader<std::fs::File>>>,
    delim: char,
    idx: VarFieldIdx
}


impl DelimitedVariantsReader {
    pub fn new(filename: &str, delim: char, has_header: bool, idx: VarFieldIdx)
        -> DelimitedVariantsReader
    {
        // Read the file.
        let f = File::open(filename)
            .expect(&format!("Couldn't open file: {:?}", filename));

        let mut iter = BufReader::new(f).lines();

        // Skip header if needed.
        if has_header {
            iter.next();
        }

        DelimitedVariantsReader {
            iter: Box::new(iter),
            delim: delim,
            idx: idx
        }
    }
}


impl Iterator for DelimitedVariantsReader {
    type Item = OrderedAllelesVariant;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.iter.next() {
            // Split the string.
            let line = s.unwrap();
            let fields = Vec::from_iter(line.split(self.delim));

            // Upper alleles.
            let a1 = fields[self.idx.a1].to_string().to_uppercase();
            let a2 = fields[self.idx.a2].to_string().to_uppercase();

            // Parse the variant.
            let v = Variant::new(
                fields[self.idx.name].to_string(),
                fields[self.idx.chrom].to_string(),
                fields[self.idx.pos].parse().unwrap(),   
                (a1.clone(), a2)
            );

            let a1_idx = if &v.alleles.0 == &a1 { 0 } else { 1 };

            return Some(
                OrderedAllelesVariant { variant: v, a1_idx: a1_idx }
            );
        }
        None
    }
}


// This is simply a struct where the a1 allele is identified.
// It can be used for arbitrary cases where we need to remember the order
// of the alleles (e.g. which one is deleterious, which one is the coded
// allele in a statistical model, etc.)
//
// For instance, it is used when creating an iterator of variants from
// a tab-delimited file so that it is possible to known which field
// the alleles came from if needed.
pub struct OrderedAllelesVariant {
    pub variant: Variant,
    pub a1_idx: u8,
}


#[derive(PartialEq, Clone, Hash, Debug)]
pub struct Chromosome {
    pub name: String
}


impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}


#[derive(Clone, Debug)]
#[repr(C)]
pub struct Variant {
    pub name: String,
    pub chrom: Chromosome,
    pub position: u32,
    pub alleles: (String, String)
}


impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "<Variant chr{}:{}_{:?}>", self.chrom, self.position,
               self.alleles)
    }
}

impl Variant {

    pub fn new(
        name: String,
        chrom: String,
        pos: u32,
        alleles: (String, String)
    ) -> Variant {

        // We use ordered uppercase alleles.
        let uc_alleles = match alleles {
            (a1, a2) => order_alleles(a1.to_uppercase(), a2.to_uppercase())
        };

        Variant {
            name: name,
            chrom: Chromosome { name: chrom },
            position: pos,
            alleles: uc_alleles
        }

    }

    pub fn alleles_ambiguous(&self) -> bool {
        match &self.alleles {
            (a1, a2) => {
                (a1 == "C" && a2 == "G") ||
                (a1 == "A" && a2 == "T")
            }
        }
    }

    pub fn primitive_locus_eq(&self, chrom: &String, pos: u32) -> bool {
        self.chrom.name == *chrom && self.position == pos
    }

    pub fn locus_eq(&self, other: &Variant) -> bool {
        self.primitive_locus_eq(&other.chrom.name, other.position)
    }

    pub fn alleles_eq(&self, other: &Variant) -> bool {
        self.alleles == other.alleles
    }

    pub fn complement_alleles(&mut self) {
        match &self.alleles {
            (a1, a2) => {
                let new_a1 = complement(&a1);
                let new_a2 = complement(&a2);

                self.alleles = order_alleles(new_a1, new_a2);
            }
        }
    }

    pub fn alleles_set(&self) -> HashSet<String> {
        let mut alleles = HashSet::new();
        alleles.insert(self.alleles.0.clone());
        alleles.insert(self.alleles.1.clone());

        alleles
    }

    pub fn get_hash(&self) -> u64 {
        let mut s = DefaultHasher::new();
        self.hash(&mut s);
        s.finish()
    }

}


impl PartialEq for Variant {
    // It would have been possible to use hash instead but this is likely
    // to be faster.
    fn eq(&self, other: &Variant) -> bool {
        let locus_match = self.locus_eq(other);
        let mut alleles_match = false;

        // We need the alleles to be either equal or complemented to equal.
        let v_alleles = self.alleles_set();
        let mut o_alleles = other.alleles_set();

        if v_alleles == o_alleles {
            alleles_match = true;
        }

        // Look at complementary.
        o_alleles = HashSet::from_iter(o_alleles.iter().map(|s| {
            complement(&s)
        }));

        if v_alleles == o_alleles {
            alleles_match = true;
        }

        locus_match && alleles_match
    }
}

impl Eq for Variant {}

impl Hash for Variant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // We create a hashable version of the alleles (as a standardized
        // string) that includes alleles from the other strand as well.
        let mut alleles = vec![
            self.alleles.0.clone(),
            self.alleles.1.clone(),
            complement(&self.alleles.0),
            complement(&self.alleles.1)
        ];

        alleles.sort();
        alleles.dedup();

        let hashable_alleles = alleles.join(",");

        self.chrom.hash(state);
        self.position.hash(state);
        hashable_alleles.hash(state);
    }
}


pub fn complement(s: &String) -> String {
    String::from_iter(s.chars().map(|c| {
        match c {
            'T' => 'A',
            'A' => 'T',
            'G' => 'C',
            'C' => 'G',
            _ => c
        }
    }))
}


#[derive(Debug)]
pub struct Genotypes {
    pub variant: Variant,
    pub genotypes: Vec<Option<u8>>,
    coded_idx: u8
}


impl Genotypes {
    pub fn new(variant: Variant, genotypes: Vec<Option<u8>>, coded_allele: &str)
        -> Genotypes {
            // Find the coded allele.
            let coded = String::from(coded_allele).to_uppercase();
            let coded_idx: u8 = if variant.alleles.0 == coded {
                0
            } else if variant.alleles.1 == coded {
                1
            } else {
                panic!("Coded allele `{}` is not an allele of `{}`",
                       coded_allele, &variant);
            };

            Genotypes { variant, genotypes, coded_idx }
    }

    pub fn coded_freq(&self) -> f64 {
        let n = self.genotypes.len() as f64;
        let sum: f64 = self.genotypes.iter()
            .filter_map(|g| *g)
            .map(f64::from)
            .sum();

        sum / (2.0 * n)
    }

    pub fn maf(&self) -> f64 {
        let freq = self.coded_freq();
        freq.min(1.0 - freq)
    }
}

impl PartialEq for Genotypes {
    fn eq(&self, other: &Genotypes) -> bool {
        (self.variant == other.variant) &&
        (self.genotypes == other.genotypes) &&
        (self.coded_idx == other.coded_idx)
    }
}


fn order_alleles(a1: String, a2: String) -> (String, String) {
    if a1.len() == a2.len() {
        // Order alphabetically.
        if a1 < a2 {
            return (a1, a2);
        }

        else {
            return (a2, a1);
        }
    }

    else {
        // Order wrt length.
        if a1.len() < a2.len() {
            return (a1, a2);
        }

        else {
            return (a2, a1);
        }

    }
}