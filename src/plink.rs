/**
 * Utilities to read plink files.
 */

use crate::core::{VarFieldIdx, DelimitedVariantsReader};

pub struct BimReader;
impl BimReader {
    pub fn new(filename: &str) -> DelimitedVariantsReader {
        let idx = VarFieldIdx {
            delimiter: '\t',
            name: 1,
            chrom: 0,
            pos: 3,
            a1: 4,
            a2: 5
        };

        DelimitedVariantsReader::new(filename, '\t', false, idx)
    }
}

