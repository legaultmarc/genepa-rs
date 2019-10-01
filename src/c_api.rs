use std::f32::NAN;
use std::os::raw::{c_char, c_uint, c_float};
use std::ffi::CStr;

use crate::core::{Variant, Genotypes};
use crate::plink::PlinkReader;


fn c_string_to_string(ptr: *const c_char) -> String {
    unsafe {
        CStr::from_ptr(ptr).to_string_lossy().into_owned()
    }
}


#[no_mangle]
pub extern "C" fn variant_new(
        name: *const c_char,
        chrom: *const c_char,
        position: c_uint,
        a1: *const c_char,
        a2: *const c_char
    ) -> *mut Variant {

    let v = Box::new(
        Variant::new(
            c_string_to_string(name),
            c_string_to_string(chrom),
            position,
            (c_string_to_string(a1), c_string_to_string(a2))
        )
    );

    Box::into_raw(v)
}

#[no_mangle]
pub extern "C" fn variant_free(ptr: *mut Variant) {
    if ptr.is_null() { return; }

    unsafe {
        // We're letting the object go out of scope so it will get cleaned by
        // rust.
        Box::from_raw(ptr);
    }
}


#[no_mangle]
pub extern "C" fn variant_print(ptr: *mut Variant) {
    unsafe {
        println!("{}", *ptr);
    }
}


#[no_mangle]
pub extern "C" fn variant_complement_alleles(ptr: *mut Variant) {
    unsafe {
        match ptr.as_mut() {
            Some(v) => v.complement_alleles(),
            None => println!("Complement alleles failed")
        }
    }
}


/**
 * Genotypes bindings.
 **/
#[no_mangle]
pub extern "C" fn genotypes_print(ptr: *const Genotypes) {
    let geno: &Genotypes = unsafe {
        ptr.as_ref().unwrap()
    };

    println!("{}", geno.variant);
    println!("{:?}", geno);
}


#[no_mangle]
pub extern "C" fn genotypes_get_variant(ptr: *const Genotypes)
    -> *const Variant {
    let variant = unsafe {
        ptr.as_ref().unwrap().variant.clone()
    };

    Box::into_raw(Box::new(variant))
}


#[no_mangle]
pub extern "C" fn genotypes_get_genotypes(ptr: *const Genotypes)
    -> *const c_float {

    let genotypes = unsafe {
        &ptr.as_ref().unwrap().genotypes
    };

    // We convert to a vector of c float
    let v: Vec<c_float> = genotypes
        .iter()
        .map(|g| match g {
            Some(g_byte) => c_float::from(*g_byte),
            None => std::f32::NAN
        })
        .collect();

    let buf = v.into_boxed_slice();
    let data = buf.as_ptr();
    std::mem::forget(buf);

    data

}

#[no_mangle]
pub extern "C" fn genotypes_free(ptr: *mut Genotypes) {
    if ptr.is_null() { return; }

    unsafe {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub extern "C" fn genotypes_data_free(ptr: *mut c_float) {
    if ptr.is_null() { return; }

    unsafe {
        Box::from_raw(ptr);
    }
}

#[no_mangle]
pub extern "C" fn genotypes_len(ptr: *mut Genotypes) -> usize {
    unsafe {
        ptr.as_ref().unwrap().genotypes.len()
    }
}

#[no_mangle]
pub extern "C" fn genotypes_maf(ptr: *mut Genotypes) -> f64 {
    unsafe {
        ptr.as_ref().unwrap().maf()
    }
}


/**
 * Plink reader bindings.
 **/

#[no_mangle]
pub extern "C" fn plink_reader_new(
    prefix: *const c_char
) -> *mut PlinkReader {

    let prefix_str = unsafe {
        CStr::from_ptr(prefix).to_str().expect(
            "Could not parse str from pointer (plink reader new)."
        )
    };

    let reader = Box::new(PlinkReader::new(prefix_str));

    Box::into_raw(reader)

}


#[no_mangle]
pub extern "C" fn plink_reader_free(ptr: *mut PlinkReader) {
    if ptr.is_null() { return; }

    unsafe {
        Box::from_raw(ptr);
    }
}


#[no_mangle]
pub extern "C" fn plink_reader_next(ptr: *mut PlinkReader) -> *mut Genotypes {
    if ptr.is_null() {
        println!("Got a null pointer (plink_reader_next).");
        return std::ptr::null_mut();
    }

    let reader = unsafe {
        ptr.as_mut()
    }.unwrap();

    match reader.next() {
        Some(g) => Box::into_raw(Box::new(g)),
        None => std::ptr::null_mut()
    }

}


/*
#[no_mangle]
pub extern "C" fn plink_reader_print(ptr: *mut PlinkReader) {
    println!("<PlinkReader - rs instance>");
}
*/


