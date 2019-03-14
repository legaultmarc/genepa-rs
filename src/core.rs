#![allow(dead_code)]

use std::fmt;
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;


#[derive(PartialEq, Clone, Hash)]
pub struct Chromosome {
    name: String
}


impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}


#[derive(Clone)]
#[repr(C)]
pub struct Variant {
    name: String,
    chrom: Chromosome,
    position: u32,
    alleles: (String, String)
}


impl fmt::Display for Variant {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "<Variant chr{}:{}_{:?}>", self.chrom, self.position,
               self.alleles)
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
