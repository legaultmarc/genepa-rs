use std::os::raw::{c_char, c_uint};
use std::ffi::CStr;
use std::ptr;

use crate::core::Variant;


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
