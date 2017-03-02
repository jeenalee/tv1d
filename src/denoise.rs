//! Total Variation 1D Denoising Algorithms.
//!

#![deny(missing_docs)]

use num;
use std::cmp;
use std::fmt;
use std::iter;
use std::ops;
use utils::*;

// TODO change tautstring to be generic
/// TODO
/// Note: numerical blow-up if floattype=float for N>=10^6
/// because the algorithm is based on the running sum of the signal
/// values.
// pub fn tautstring_denoise(input: &Vec<f64>, lambda: f64) -> Vec<f64> {
//     assert!(input.len() > 0, "Input list should have at least one value.");
//     let mut output = vec![0.0; input.len()];
//     let width = input.len() + 1;

//     let mut index_low = vec![0; width];
//     let mut slope_low = vec![0.0; width];
//     let mut index_up = vec![0; width];
//     let mut slope_up = vec![0.0; width];
//     let mut index = vec![0; width];

//     let mut z = vec![0.0; width];
//     // Lower and upper boundaries
//     let mut y_low = vec![0.0; width];
//     let mut y_up = vec![0.0; width];
//     let mut s_low = 0;
//     let mut c_low = 0;
//     let mut s_up = 0;
//     let mut c_up = 0;
//     let mut c = 0;
//     let mut i = 2;

//     // The boundaries are first defined by the input values.
//     y_low[1] = input[0] - lambda;
//     y_up[1] = input[0] + lambda;
//     // Based on the input values, get culmulative sums.
//     while i < width {
//         y_low[i] = y_low[i - 1] + input[i - 1];
//         y_up[i] = y_up[i - 1] + input[i - 1];
//         i += 1;
//     }

//     // The last value's lower boundary is elevated by lambda.
//     y_low[width - 1] += lambda;
//     // The last value's lower boundary is decreased by lambda.
//     y_up[width - 1] -= lambda;

//     slope_low[0] = f64::INFINITY;
//     slope_up[0] = f64::NEG_INFINITY;

//     // z is first set to be the first lower boundary.
//     z[0] = y_low[0];

//     for i in 1..width {
//         c_low += 1;
//         c_up += 1;

//         index_low[c_low] = i;
//         index_up[c_up] = i;
//         slope_low[c_low] = y_low[i] - y_low[i - 1];

//         // c_low is too large. Decrease it by one.
//         while (c_low > s_low + 1) && (slope_low[cmp::max(s_low, c_low - 1)]
// <= slope_low[c_low]) {
//             c_low -= 1;
//             index_low[c_low] = i;
//             if c_low > s_low + 1 {
//                 slope_low[c_low] = (y_low[i] - y_low[index_low[c_low - 1]]) /
//                     (i as f64 - index_low[c_low - 1] as f64);
//             } else {
//                 slope_low[c_low] = (y_low[i] - z[c]) /
//                     (i as f64 - index[c] as f64);
//             }
//         }

//         slope_up[c_up] = y_up[i] - y_up[i - 1];
//         while (c_up > s_up + 1) && (slope_up[cmp::max(c_up - 1, s_up)] >= slope_up[c_up]) {
//             c_up -= 1;
//             index_up[c_up] = i;
//             if c_up > s_up+1 {
//                 slope_up[c_up] = (y_up[i] - y_up[index_up[c_up - 1]]) /
//                     (i as f64 - index_up[c_up - 1] as f64);
//             } else {
//                 slope_up[c_up] = (y_up[i] - z[c]) / (i as f64 - index[c] as f64);
//             }
//         }
//         while (c_low == s_low + 1) &&
//             (c_up > s_up+1) &&
//             (slope_low[c_low] >= slope_up[s_up + 1]) {
//                 c += 1;
//                 s_up += 1;
//                 index[c] = index_up[s_up];
//                 z[c] = y_up[index[c]];
//                 index_low[s_low] = index[c];
//                 slope_low[c_low] = (y_low[i] - z[c]) / (i as f64 - index[c] as f64);
//             }
//         while (c_up == s_up + 1) &&
//             (c_low > s_low + 1) &&
//             (slope_up[c_up] <= slope_low[s_low + 1]) {
//                 c += 1;
//                 s_low += 1;
//                 index[c] = index_low[s_low];
//                 z[c] = y_low[index[c]];
//                 index_up[s_up] = index[c];
//                 slope_up[c_up] = (y_up[i] - z[c]) / (i as f64 - index[c] as f64);
//             }
//     }
//     for i in 1..(c_low - s_low + 1) {
//         index[c + i] = index_low[s_low + i];
//         z[c + i]=y_low[index[c + i]];
//     }
//     c += c_low - s_low;
//     let mut j = 0;
//     let mut a;
//     i = 1;

//     while i <= c {
//         a = (z[i] - z[i - 1]) / (index[i] as f64 - index[i - 1] as f64);
//         while j < index[i] {
//             output[j] = a;
//             j += 1;
//         }
//         i += 1;
//     }
//     output
// }
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
