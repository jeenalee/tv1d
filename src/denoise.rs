//! Total Variation 1D Denoising Algorithms.
//!

#![deny(missing_docs)]

use num;
use std::cmp;
use std::fmt;
use std::iter;
use std::ops;
use utils::*;

/// Denoise an array of data points given a `lambda`, which decides
/// the degree of denoising.
///
/// A very large `lambda` will create a denoised output, in which
/// every value is equivalent to the average of the input data points.
/// A very small `lambda` (i.e. 0) will not denoise the input at all,
/// creating an output that is the same as the input.
///
/// This algorithm was first described by Laurent Condat in 2013 in
/// his paper "A Direct Algorithm for 1D Total Variation Denoising."
pub fn condat_denoise<T>(input: &Vec<T>, lambda: T) -> Vec<T>
    where T: num::Num + num::FromPrimitive + cmp::PartialOrd
    + ops::Neg + ops::Add<Output=T> + ops::Neg<Output=T>
    + ops::AddAssign<T> + Sized + fmt::Debug + Copy
{
    assert!(input.len() > 0,
            "Input list should have at least one value.");

    let width = input.len();
    let mut output = Vec::new();

    // k: current sample location
    let mut k = 0;
    // k0: beginning of current segment
    let mut k0 = 0;

    let twolambda = T::from_u8(2).expect("Unable to transform `2` to T.") * lambda;
    let minlambda = -lambda;
    // umin and umax are used for keeping track of previous data
    // points we have seen.
    let mut umin = lambda;
    let mut umax = minlambda;
    // boundaries of the segment's value
    let mut segment_lower_bound = input[0] - lambda;
    let mut segment_upper_bound = input[0] + lambda;

    // last position where umax = -lambda
    let mut kplus = 0;
    // last position where umin = lambda
    let mut kminus = 0;

    loop {
        if k == (width - 1) {
            if umin < num::zero() {
                // if segment_lower_bound is too high,
                // jump down.
                output.extend(iter::repeat(segment_lower_bound).take(kminus - k0 + 1));
                k0 = kminus + 1;
                sync_values(k0, &mut [&mut k, &mut kminus]);
                segment_lower_bound = input[kminus];
                umin = lambda;
                umax = segment_lower_bound + umin - segment_upper_bound;
            } else if umax > num::zero() {
                // if segment_upper_bound is too low,
                // jump up.
                output.extend(iter::repeat(segment_upper_bound).take(kplus - k0 + 1));
                k0 = kplus + 1;
                sync_values(k0, &mut [&mut k, &mut kplus]);
                segment_upper_bound = input[kplus];
                umax = minlambda;
                umin = segment_upper_bound + umax - segment_lower_bound;
            } else {
                // segment_lower_bound and segment_upper_bound is appropriate.
                segment_lower_bound += umin /
                                       num::FromPrimitive::from_usize(k - k0 + 1)
                    .expect("Unable to convert usize to num::FromPrimitive.");
                output.extend(iter::repeat(segment_lower_bound).take(k - k0 + 1));
                return output;
            }
        } else {
            umin += input[k + 1] - segment_lower_bound; // input 0 and 1
            umax += input[k + 1] - segment_upper_bound; // input 0 and 1
            if umin < minlambda {
                output.extend(iter::repeat(segment_lower_bound).take(kminus - k0 + 1));
                k0 = kminus + 1;
                sync_values(k0, &mut [&mut k, &mut kminus, &mut kplus]);
                segment_lower_bound = input[kplus];
                segment_upper_bound = segment_lower_bound + twolambda;
                umin = lambda;
                umax = minlambda;
            } else if umax > lambda {
                output.extend(iter::repeat(segment_upper_bound).take(kplus - k0 + 1));
                k0 = kplus + 1;
                sync_values(k0, &mut [&mut k, &mut kminus, &mut kplus]);
                segment_upper_bound = input[kplus];
                segment_lower_bound = segment_upper_bound - twolambda;
                umin = lambda;
                umax = minlambda;
            } else {
                k += 1;
                if umin >= lambda {
                    kminus = k;
                    segment_lower_bound = (segment_lower_bound + (umin - lambda)) /
                                          num::FromPrimitive::from_usize(kminus - k0 + 1)
                        .expect("Unable to convert usize to num::FromPrimitive.");
                    umin = lambda;
                }
                if umax <= minlambda {
                    kplus = k;
                    segment_upper_bound += (umax + lambda) /
                                           num::FromPrimitive::from_usize(kplus - k0 + 1)
                        .expect("Unable to convert usize to num::FromPrimitive.");
                    umax = minlambda;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tautstring_test_input_output_length() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = tautstring_denoise(&input, 0.0);
        assert_eq!(input.len(), output.len());
    }

    #[test]
    fn tautstring_test_zero_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = tautstring_denoise(&input, 0.0);
        let output_expected = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.0001);
        }
    }

    #[test]
    fn tautstring_test_large_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = tautstring_denoise(&input, 100.0);
        // The expected output is taken from the Laurent Condat's C
        // implementation.
        let output_expected = vec![4.925, 4.925, 4.925, 4.925, 4.925, 4.925, 4.925, 4.925];
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.0001);
        }
    }

    #[test]
    fn tautstring_test_moderate_lambda() {
        let input = vec![111.0, 422.1, 145.2, 248.2, 871.4, 675.2, 436.2, 310.1];
        let output = tautstring_denoise(&input, 5.0);
        // The expected output is taken from the Laurent Condat's C
        // implementation.
        let output_expected = vec![116.0, 412.100006, 155.199997, 248.199997, 861.400024,
                                   675.200012, 436.200012, 315.100006];
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.0001);
        }
    }

    #[test]
    fn condat_test_input_output_length() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat_denoise(&input, 0.0);
        assert_eq!(input.len(), output.len());
    }

    #[test]
    fn condat_test_zero_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat_denoise(&input, 0.0);
        let output_expected = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        assert_eq!(output, output_expected);
    }

    #[test]
    fn condat_test_large_lambda() {
        let input = vec![111.0, 422.1, 145.2, 248.2, 871.4, 675.2, 436.2, 310.1];
        let output = condat_denoise(&input, 10000.0);
        // The expected output is taken from the Laurent Condat's C
        // implementation.
        let output_expected = vec![402.425049, 402.425049, 402.425049, 402.425049, 402.425049,
                                   402.425049, 402.425049, 402.425049];
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.0001);
        }
    }

    #[test]
    fn condat_test_moderate_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat_denoise(&input, 3.0);
        // The expected output is taken from the Laurent Condat's C
        // implementation.
        let output_expected = vec![3.050000, 3.050000, 4.933333, 4.933333, 4.933333, 5.200000,
                                   6.200000, 7.100000];
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.000001);
        }
    }

    #[test]
    #[should_panic]
    fn condat_test_negative_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        condat_denoise(&input, -1.0);
    }

    #[test]
    #[should_panic]
    fn condat_test_empty_input() {
        let input = vec![];
        condat_denoise(&input, 1.0);
    }
}
