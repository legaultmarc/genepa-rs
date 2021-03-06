mod core;
mod c_api;

pub mod plink;
pub mod utils;

pub use crate::c_api::*;
pub use crate::core::{Variant, OrderedAllelesVariant};
