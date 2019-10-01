use ndarray::ArrayViewMut;

use crate::core::Genotypes;

pub fn compute_ld(mut g: Genotypes, mut other_genotypes: Vec<Genotypes>, r2: bool)
    -> Vec<f64>
{
    let n_samples = g.genotypes.len();
    let n_variants = other_genotypes.len();

    let mut other_geno_data: Vec<&Option<u8>> = other_genotypes.iter_mut()
        .flat_map(|g| &g.genotypes)
        .collect();

    let geno_arr = ArrayViewMut::from_shape(
        (1, n_samples),
        &mut g.genotypes
    ).unwrap();

    let geno_arr_others = ArrayViewMut::from_shape(
        (n_samples, n_variants),
        &mut other_geno_data
    ).unwrap();

    // Example computing maf.
    // acc is n_alleles, n_samples
    let (n_alleles, n_samples) = geno_arr.fold((0, 0), |acc, x| {
        match x {
            Some(geno) => {
                (
                    (acc.0 + *geno as u64),
                    acc.1 + 1
                )
            },
            _ => acc
        }
    });

    println!("{:?}", n_alleles);
    println!("{:?}", n_samples);
    println!("{:?}", n_alleles as f64 / (2.0 * n_samples as f64));

    // TODO
    vec![3.2]
}

#[cfg(test)]
mod tests {
    use crate::plink::PlinkReader;
    use crate::core::{Chromosome, Variant};
    use super::*;

    #[test]
    fn test_test() {
        let mut plink = PlinkReader::new(
            "./test_data/common_extracted_1kg.missing"
        );

        let v = Variant::new(
            "rs77883301".to_string(),
            "16".to_string(),
            56647136,
            ("G".to_string(), "A".to_string())
        );
        let g = plink.get_variant_genotypes(&v).unwrap();

        let chr = Chromosome { name: "16".to_string() };
        let other_geno = plink.get_variants_in_region(
            &chr,
            56651160,
            56899006
        );

        compute_ld(g, other_geno, true);

        assert!(false);
    }
}