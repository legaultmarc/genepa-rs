mod core;

use crate::core::Variant;
use crate::core::complement;

fn main() {

    let v1 = Variant::new(
        String::from("rs12345"),
        String::from("1"),
        12345,
        (String::from("G"), String::from("C"))
    );

    let v1_2 = Variant::new(
        String::from("rs12345"),
        String::from("1"),
        12345,
        (String::from("C"), String::from("G"))
    );

    let v2 = Variant::new(
        String::from("rs9471841"),
        String::from("15"),
        9414141,
        (String::from("T"), String::from("G"))
    );

    let mut v3 = Variant::new(
        String::from("rs841719831"),
        String::from("20"),
        519731741,
        (String::from("T"), String::from("TCC"))
    );

    println!("{}", v1);
    println!("{}", v1_2);
    println!("{}", v2);
    println!("{}", v3);
    println!("{:?}", v1.alleles_ambiguous());
    println!("{:?}", v1_2.alleles_ambiguous());
    println!("{:?}", v2.alleles_ambiguous());

    v3.complement_alleles();
    println!("{}", v3);

    assert_eq!(
        String::from("CAGCATGNNNX"),
        complement(&String::from("GTCGTACNNNX"))
    );

    let mut v3_1 = v3.clone();
    v3_1.complement_alleles();
    println!("{}, {} : equal? {:?}", v3, v3_1, v3 == v3_1);
    println!("{}, {} : equal? {:?}", v3, v2, v2 == v3);

    println!("{:?} {:?}", v3.get_hash(), v3_1.get_hash());

}
