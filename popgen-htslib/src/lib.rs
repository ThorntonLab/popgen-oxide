//! Adapter types for [`rust_htslib`]

use popgen::AlleleID;
use std::ffi::c_void;

/// Error type
#[non_exhaustive]
#[derive(Debug)]
pub enum Error {
    // NOTE: this is a bad name...
    /// Encapsulation of errors from [`rust_htslib`]
    RecordError(rust_htslib::errors::Error),
    /// Integer error codes from the htslib C API
    ErrorCode(i32),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RecordError(e) => write!(f, "{e}"),
            Self::ErrorCode(code) => write!(f, "error code: {code}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<rust_htslib::errors::Error> for Error {
    fn from(value: rust_htslib::errors::Error) -> Self {
        Self::RecordError(value)
    }
}

struct GenotypesAdapter<T>(T, *mut c_void);

impl<T, I> Iterator for GenotypesAdapter<T>
where
    T: Iterator<Item = I>,
{
    type Item = I;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}

impl<T> DoubleEndedIterator for GenotypesAdapter<T>
where
    T: DoubleEndedIterator,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.0.next_back()
    }
}

impl<T> ExactSizeIterator for GenotypesAdapter<T>
where
    T: ExactSizeIterator,
{
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<T> Drop for GenotypesAdapter<T> {
    fn drop(&mut self) {
        // Safety: we control construction of this type, and we only construct it using a buffer allocated by htslib C.
        unsafe {
            rust_htslib::htslib::free(self.1);
        }
    }
}

/// Iterator over genotypes in a [`Record`](rust_htslib::bcf::Record).
///
/// The iterator emits iterators over [`Option`] of [`popgen::AlleleID`] in each genotype.
/// The [`Option::None`] variant implies missing data.
///
///
/// Use [`Iterator::flatten`] to convert the return value into an iterator over all the
/// alleles called in the record.
pub fn bcf_record_to_genotypes_iter_adapter(
    record: &rust_htslib::bcf::Record,
) -> Result<
    impl Iterator<Item = impl DoubleEndedIterator<Item = Option<AlleleID>> + ExactSizeIterator> + '_,
    Error,
> {
    // NOTE: this implementation does not rely on the rust_htslib safe API
    // because the iterator types defined there cannot be properly flattened/aggregated.
    // Basically, the borrow checker prevents this API from being written.
    // Therefore, we work at the level of htslib C types.
    // Here, the borrow checker correctly notes that the iterator is tied to the
    // lifetime of the borrowed (rust-side) Record.
    use rust_htslib::htslib;

    // The dst and ndst parameters to bcf_get_format_values and related functions require a buffer with provided pointer and length.
    // This is malloc'd/realloc'd within htslib, but we need to hold onto the pointer ourselves.

    let mut gt_buf = (std::ptr::null_mut(), 0);
    // SAFETY: we are not using this pointer after the header is dropped.
    // (see docs of the header.as_ptr fn in the unsafe block below)
    assert!(!unsafe { record.header().as_ptr() }.is_null());

    let n_alleles_in_record = {
        // See docs at https://github.com/samtools/htslib/blob/7c5e3e7ebcdf90c8f96afd8a06d75ffa5603e417/htslib/vcf.h#L1135
        // We ask htslib to parse the GT field.
        // Within htslib, #define bcf_get_format_int32 is provided to omit the last argument.
        let format_values = unsafe {
            htslib::bcf_get_format_values(
                // The header struct bcf_hdr_t.
                record.header().as_ptr(),
                // The line/record bcf1_t.
                record.inner,
                // We want the field tagged "GT".
                c"GT".as_ptr(),
                // We pass the buffer from earlier.
                &mut gt_buf.0,
                &mut gt_buf.1,
                // We want this parsed as a collection of integers, where an integer is one of the VCF datatypes.
                htslib::BCF_HT_INT as i32,
            )
        };

        // The return value is either a negative error code or the number of values written.
        match format_values {
            -1 => {
                return Err(rust_htslib::errors::Error::BcfUndefinedTag {
                    tag: String::from("GT"),
                }
                .into())
            }
            -2 => {
                return Err(rust_htslib::errors::Error::BcfUnexpectedType {
                    tag: String::from("GT"),
                }
                .into())
            }
            -3 => {
                return Err(rust_htslib::errors::Error::BcfMissingTag {
                    tag: String::from("GT"),
                    record: record.desc(),
                }
                .into())
            }
            -4 => return Err(rust_htslib::errors::Error::BcfAllocationError.into()),
            other if other < 0 => return Err(Error::ErrorCode(other)),
            ret => ret,
        }
    };

    // https://github.com/samtools/htslib/blob/7c5e3e7ebcdf90c8f96afd8a06d75ffa5603e417/htslib/vcf.h#L1049
    // https://github.com/samtools/htslib/blob/7c5e3e7ebcdf90c8f96afd8a06d75ffa5603e417/htslib/vcf.h#L166
    // We need field n on the struct bcf_fmt_t, describing the maximum number of alleles per sample.
    let fmt_inner =
        unsafe { htslib::bcf_get_fmt(record.header().as_ptr(), record.inner, c"GT".as_ptr()) };
    assert!(!fmt_inner.is_null());

    assert!(!gt_buf.0.is_null());

    // Safety: We called for the parsing of i32 from htslib, and if the pointer is not null, it points to malloc'd memory.
    // As long as htslib correctly returns n_alleles_in_record, this slice is valid.
    let gt_iter_iter = unsafe {
        std::slice::from_raw_parts(
            gt_buf.0.cast_const().cast::<i32>(),
            n_alleles_in_record as usize,
        )
    }
    // See https://github.com/samtools/htslib/blob/7c5e3e7ebcdf90c8f96afd8a06d75ffa5603e417/vcf.c#L6157
    // If fmt_inner is not null, then it points to a newly or previously initialized bcf_fmt_t.
    .chunks(unsafe { fmt_inner.as_ref().unwrap() }.n as usize)
    .map(|s| {
        // The magic number below is
        // rust_htslib::bcf::record::VECTOR_END_INTEGER,
        // which is private.
        // It indicates the position where data ends for this sample.
        s.split(|v| *v == i32::MIN + 1)
            .next()
            .unwrap()
            .iter()
            .map(|a| {
                // This replicates code in rust_htslib::bcf::record::GenotypeAllele.
                if a > &0 {
                    Some(AlleleID::from(((*a >> 1) - 1) as usize))
                } else {
                    None
                }
            })
    });

    Ok(GenotypesAdapter(gt_iter_iter, gt_buf.0))
}
