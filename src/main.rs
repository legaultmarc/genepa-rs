mod core;
mod plink;

use crate::plink::PlinkReader;

fn main() {
    // For fun read a plink file and compute all MAFs.
    let reader = PlinkReader::new("/Users/legaultmarc/projects/StatGen/grs/test_data/big");

    for g in reader {
        println!("{}: {}", g.variant, g.maf());
    }
}
