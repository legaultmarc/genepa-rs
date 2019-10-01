#![allow(dead_code)]

/**
 * Utilities to read plink files.
 */

use std::iter::{FromIterator};
use std::path::Path;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Write, SeekFrom, Seek};
use std::fs::{File, OpenOptions};

use crate::core::{VarFieldIdx, DelimitedVariantsReader, Variant, Genotypes,
                  Chromosome};


struct BimIndex {
    filename: String,
    n_variants: u32

}

impl BimIndex {
    pub fn get_or_create_bim_index(filename: &str) -> BimIndex {
        let f = File::open(filename)
            .expect(&format!("Could not read BIM: `{}`", filename));
        let buf_reader = BufReader::new(f);

        let output_filename = String::from(filename).replace(".bim", ".bimidx.gz");

        if Path::new(&output_filename).is_file() {
            // Unfortunately, we have to know the number of variants so it
            // is necessary to count the lines.
            let n_variants = buf_reader.lines().count() as u32;

            return BimIndex { filename: output_filename, n_variants };
        }

        let output = OpenOptions::new()
            .write(true)
            .create(true)
            .open(&output_filename)
            .expect("Can't create BIM index file.");

        // Spawn the bgzip process to directly write to it.
        let mut bgzip = Command::new("bgzip")
            .stdin(Stdio::piped())
            .stdout(Stdio::from(output))
            .spawn()
            .unwrap();

        let mut bgzip_stdin = bgzip.stdin.as_mut()
            .expect("Could not get bgzip stdin.");

        // Write the index to disk.
        let mut n_variants: usize = 0;
        for (i, line) in buf_reader.lines().enumerate() {
            // TODO we could build a name index at the same time.
            write!(&mut bgzip_stdin, "{}\t{}\n", line.unwrap().as_str(), i)
                .expect("Failed writing line to BIM index.");

            n_variants = i;
        };

        match bgzip.wait() {
            Ok(status) => {
                 if !status.success() {
                     panic!("Error creating index: bgzip process returned with \
                            an error");
                 }
            },
            Err(e) => panic!("Error executing bgzip: {:?}", e)
        };

        // Tabix
        let tabix_output = Command::new("tabix")
                .args(&["-s", "1", "-b", "4", "-e", "4", &output_filename])
                .output()
                .expect("Tabix failed.");

        if !tabix_output.status.success() {
            panic!("Tabix returned an error, could not build BIM index.");
        }

        BimIndex {
            filename: output_filename,
            n_variants: n_variants as u32
        }
    }

    // Returns a vector of index, variant, coded_allele
    fn _run_tabix(&self, region: &str) -> Vec<(u32, Variant, String)> {
        let tabix = Command::new("tabix")
            .arg(&self.filename)
            .arg(region)
            .output()
            .expect("Couldn't spawn tabix for BIM variant query.");

        if !tabix.status.success() {
            panic!("Error searching the BIM index using tabix.");
        }

        String::from_utf8(tabix.stdout)
            .unwrap()
            .lines()
            .map(|line| {
                // Parse a variant.
                let vec = Vec::from_iter(line.split('\t'));

                let chrom: String = vec[0].to_string();
                let name: String = vec[1].to_string();
                let pos: u32 = vec[3].to_string().parse().unwrap();
                let a1: String = vec[4].to_string();
                let a2: String = vec[5].to_string();

                let variant = Variant::new(name, chrom, pos, (a1.clone(), a2));

                let idx: u32 = vec[6].to_string().parse().unwrap();

                (idx, variant, a1)
            })
            .collect()
    }

    fn get_region_index_and_coded(&self, chrom: &str, start: u32, end: u32)
        -> Vec<(u32, Variant, String)> {
            let region = format!("{}:{}-{}", chrom, start, end);
            self._run_tabix(&region)
        }

    fn get_variant_index_and_coded(&self, v: &Variant) -> Option<(u32, String)> {
        let region = format!("{}:{}-{}", v.chrom.name, v.position, v.position);

        let matches: Vec<(u32, Variant, String)> = self._run_tabix(&region)
            .into_iter()
            .filter(|(_, observed, _)| {
                observed == v
            })
            .collect();

        match matches.len() {
            0 => None,
            1 => {
                let mtch = &matches[0];
                // Returns index and a1.
                Some((mtch.0, mtch.2.clone()))
            },
            _ => panic!("There are duplicate variants in the bim file.")
        }
    }
}

// Read a fam into a vector of sample IDs.
fn read_fam(filename: &str) -> Vec<String> {
    let f = File::open(filename).expect("Could not open FAM");
    let reader = BufReader::new(f);

    reader
        .lines()
        .map(|l| {
            let line = l.unwrap();
            let vec = Vec::from_iter(line.split('\t'));

            vec[0].to_string()  // Return the sample id
        })
        .collect()
}


pub struct PlinkReader {
    bim_reader: DelimitedVariantsReader,
    bim_index: BimIndex,
    samples: Vec<String>,
    bed_reader: BedReader<BufReader<File>>
}

impl PlinkReader {
    pub fn new(prefix: &str) -> PlinkReader {
        // Get or create the index for the bim.
        let bim_filename = format!("{}.bim", &prefix);
        let bim_index = BimIndex::get_or_create_bim_index(&bim_filename);
        let bim_reader = BimReader::new(&bim_filename);

        let fam_filename = format!("{}.fam", &prefix);
        let samples = read_fam(&fam_filename);

        let n_samples = samples.len() as u32;

        let bed_filename = format!("{}.bed", &prefix);
        let bed_reader = BedReader::new(
            &bed_filename, n_samples, bim_index.n_variants
        );

        PlinkReader {bim_reader, bim_index, samples, bed_reader}
    }

    fn _seek_to_idx(&mut self, idx: u32) {
        let actual_seek = 3 + self.bed_reader._chunk_size * idx as usize;
        self.bed_reader.reader.seek(SeekFrom::Start(actual_seek as u64))
            .expect("Could not seek in BED");
    }

    fn _seek_and_read_to_idx(&mut self, idx: u32) -> Vec<Option<u8>> {
        self._seek_to_idx(idx);
        self.bed_reader._read_variant_chunk()
    }

    pub fn get_variant_genotypes(&mut self, v: &Variant) -> Option<Genotypes> {
        // Find in the index.
        let res = self.bim_index.get_variant_index_and_coded(&v);
        if let Some((idx, coded)) = res {
            // Read variant
            let geno_vec = self._seek_and_read_to_idx(idx);
            return Some(Genotypes::new(v.clone(), geno_vec, &coded));
        }
        None
    }

    pub fn get_variants_in_region(&mut self, chrom: &Chromosome, start: u32,
                                  end: u32)
        -> Vec<Genotypes>
    {
        // Do a region query on the BIM index.
        self.bim_index.get_region_index_and_coded(&chrom.name, start, end)
            .into_iter()
            .map(|(idx, v, coded)| {
                // For every index, variant and coded, read the genotypes.
                let geno_vec = self._seek_and_read_to_idx(idx);

                Genotypes::new(v.clone(), geno_vec, &coded)
            })
            .collect()
    }
}


impl Iterator for PlinkReader {
    type Item = Genotypes;

    fn next(&mut self) -> Option<Self::Item> {
        // Check if we're at the end of the file.
        // let file = self.bed_reader.reader.get_mut();
        // let tell = file.seek(SeekFrom::Current(0)).unwrap();
        // let chunk_size = self.bed_reader._chunk_size as u64;
        
        match self.bim_reader.next() {
            // oav is ordered alleles variant.
            Some(ref oav) => {
                let geno_vec = self.bed_reader._read_variant_chunk();

                let coded_allele =  if oav.a1_idx == 0 {
                    &oav.variant.alleles.0
                } else if oav.a1_idx == 1 {
                    &oav.variant.alleles.1
                } else {
                    panic!("Problem with the ordered allele variant index");
                };

                Some(Genotypes::new(
                    oav.variant.to_owned(),
                    geno_vec,
                    &coded_allele)
                )
            }
            None => None
        }
    }
}


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


struct BedReader<T: BufRead> {
    reader: T,
    n_samples: u32,
    n_variants: u32,
    _chunk_size: usize
}

impl BedReader<BufReader<File>> {
    pub fn new(filename: &str, n_samples: u32, n_variants: u32)
        -> BedReader<BufReader<File>> 
    {
        let f = File::open(filename).unwrap();
        BedReader::new_from_reader(BufReader::new(f), n_samples, n_variants)
    }

    fn get_chunk_size(n_samples: u32) -> usize {
        (f64::from(n_samples) / 4.0).ceil() as usize
    }
}

impl<T: BufRead> BedReader<T> {
    pub fn new_from_reader(reader: T, n_samples: u32, n_variants: u32)
        -> BedReader<T>
    {
        let mut bed_reader = BedReader {
            reader,
            n_samples,
            n_variants,
            _chunk_size: BedReader::get_chunk_size(n_samples)
        };

        if !&bed_reader._verify_magic_number() {
            panic!("The provided file is not in the BED format (according to \
                   the magic number)");
        }

        bed_reader
    }

    fn _read_variant_chunk(&mut self) -> Vec<Option<u8>> {
        let n_samples = self.n_samples as usize;

        let mut buf_vec: Vec<u8> = vec![0; self._chunk_size];
        self.reader.read_exact(&mut buf_vec)
            .expect("Could not read bytes.");

        let mask: u8 = 0b11;
        let mut genotypes: Vec<Option<u8>> = buf_vec
          .iter()
          .map(|b| {
            // Every byte has the information on up to 4 samples.
            // DD CC BB AA
            // We use bitshifts to bring the current sample to the lowest bits
            // and the mask to extract them.
            let cur_geno: Vec<Option<u8>> = (0..=6).step_by(2).map(|shft| {
                let coded_geno = (b >> shft as u8) & mask;
                match coded_geno {
                    0 => Some(2), // Homo A1
                    1 => None,    // NA
                    2 => Some(1), // Hetero
                    3 => Some(0),  // Homo A2
                    _ => panic!("Unexpected value in bed file.")
                }
            }).collect();
            cur_geno
          })
          .flatten()
          .collect();

        // It is possible that the last data is not relevant.
        if genotypes.len() > n_samples {
            genotypes.truncate(n_samples)
        }

        genotypes
    }

    fn _verify_magic_number(&mut self) -> bool {
        // Make sure the first 3 bytes are 0x6c, 0x1b, 0x01.
        let mut first_3_bytes = [0; 3];
        self.reader.read_exact(&mut first_3_bytes).unwrap();

        (
            first_3_bytes[0] == 0x6c &&
            first_3_bytes[1] == 0x1b &&
            first_3_bytes[2] == 0x01
        )
    }
}


#[cfg(test)]
mod tests {

    use std::io::BufReader;
    use super::*;

    fn get_example_bed() -> BufReader<&'static [u8]> {
        let bed = include_bytes!(
            "../test_data/common_extracted_1kg.missing.bed"
        );
        let slice: &[u8] = bed;
        BufReader::new(slice)
    }

    fn read_genotypes(filename: &str) -> Vec<Option<u8>> {
        let expct_f = File::open(filename).expect(
            "Could not read file containing genotypes."
        );

        BufReader::new(expct_f).lines().map(|l| {
            let line = l.unwrap();
            match line.as_str() {
                "0.0" => Some(0 as u8),
                "1.0" => Some(1 as u8),
                "2.0" => Some(2 as u8),
                "nan" => None,
                x => { panic!("Unexpected symbol: `{}`", x); }
            }
        }).collect()
    }

    #[test]
    fn test_constructor_from_reader() {
        let bytes_reader = get_example_bed();
        BedReader::new_from_reader(
            bytes_reader,
            503,
            11158
        );
    }

    #[test]
    fn test_constructor_from_file() {
        BedReader::new(
            "test_data/common_extracted_1kg.missing.bed",
            503,
            11158
        );
    }

    #[test]
    fn test_create_bim_index() {
        // TODO
        BimIndex::get_or_create_bim_index(
            "test_data/common_extracted_1kg.missing.bim"
        );
    }

    #[test]
    fn test_read_variant_genotypes() {
        let mut reader = PlinkReader::new(
            "test_data/common_extracted_1kg.missing"
        );

        let v = Variant::new(
            "rs1610216".to_string(),
            "16".to_string(),
            56642284,
            ("G".to_string(),
             "A".to_string())
        );

        let expct_geno_vec = read_genotypes(
            "test_data/rs1610216_16_56642284_a_g.genotypes.txt"
        );

        let expct_geno = Genotypes::new(v, expct_geno_vec, "G");

        assert_eq!(
            expct_geno,
            reader.get_variant_genotypes(&expct_geno.variant).unwrap()
        );
    }

/*
    #[test]
    fn cur() {
        let mut bed = BedReader::new_from_reader(
            get_example_bed(),
            503,
            11158
        );
        bed._verify_magic_number();
        let genotypes = bed._read_variant_chunk();
        println!("{:?}", genotypes);
        assert!(false);
    }
*/
}