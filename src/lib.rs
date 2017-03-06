//! A collection of total variation denoising algorithms for 1D data.

#![deny(missing_docs)]
#![deny(missing_debug_implementations)]

extern crate num;

mod utils;

use std::cmp;
use std::iter;
use std::ops;

/// Denoises the input values based on a tautstring algorithm by
/// Davies P. and Kovac A. in 2001 in the paper ["Local extremes, runs, strings
/// and
/// multiresolution"](https://pure.tue.nl/ws/files/2200362/Metis214115.pdf).
///
/// The algorithm was implemented by Condat L. in C, which was then
/// implemented in Rust here. The algorithm can be understood as
/// assuming the input values as a string of data points, which is
/// transformed taut.
///
/// Note that this algorithm is based on the running sum of the input
/// values. If the input values are very large or there are many of
/// them, this function will fail. The input must be a slice of
/// floats, not integers.
///
/// `lambda` closer to `0` means the denoised
/// output will resemble the input more. As `lambda` increases, the
/// denoised output values become closer to the average of the input
/// values.
///
/// # Panics
/// Panics if input vector's length is `0`.
///
/// # Examples
///
/// With `lambda` equal to `0`, the denoised output will be the same as the input:
///
/// ```
/// let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
///
/// let denoised_with_zero_lambda = tautstring(input, 0.0);
/// assert_eq!(denoised_with_zero_lambda, vec![1.0, 2.0, 3.0, 4.0, 5.0])
/// ```
///
/// With larger `lambda`, the denoised output becomes closer to the average of the input:
///
/// ```
/// let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
///
/// let denoised_with_larger_lambda = tautstring(input, 10.0);
/// assert_eq!(denoised_with_larger_lambda, vec![3.0, 3.0, 3.0, 3.0, 3.0]);
/// ```
///
pub fn tautstring<T>(input: &[T], lambda: T) -> Vec<T>
    where T: num::Num + num::FromPrimitive + cmp::PartialOrd
    + ops::AddAssign<T> + ops::SubAssign<T>  + num::Float + num::ToPrimitive
{
    assert!(input.len() > 0,
            "Input list should have at least one value.");
    let mut output = vec![num::zero(); input.len()];
    let width = input.len() + 1;

    // Vectors for keeping track of indices.
    let mut index = vec![num::zero(); width];
    let mut index_low = vec![num::zero(); width];
    let mut index_up = vec![num::zero(); width];

    // `slope_low` and `slope_up` is used to store the slope between
    // consecutive input values.
    let mut slope_low = vec![num::zero(); width];
    let mut slope_up = vec![num::zero(); width];

    // `z` stores either `lower_boundary` or `upper_boundary`
    // throughout the program, which will be used as the denoised
    // output at the end of the program.
    let mut z = vec![num::zero(); width];

    // `lower_bound` and `upper_bound` first stores the
    // cumulative sums of the input values. This will be used to find
    // slopes between each input points, and later in the denoising
    // step will be used as the denoised output.
    let mut lower_bound = vec![num::zero(); width];
    let mut upper_bound = vec![num::zero(); width];

    let mut s_low = num::zero();
    let mut c_low = 0;
    let mut s_up = 0;
    let mut c_up = 0;
    let mut c = 0;

    // First define `lower_bound` and `upper_bound` by the
    // first input value and lambda.
    lower_bound[1] = input[0] - lambda;
    upper_bound[1] = input[0] + lambda;

    // Generic tracker.
    let mut i = 2;

    // Get culmulative sum of the input values.
    while i < width {
        lower_bound[i] = lower_bound[i - 1] + input[i - 1];
        upper_bound[i] = upper_bound[i - 1] + input[i - 1];
        i += 1;
    }

    lower_bound[width - 1] += lambda;
    upper_bound[width - 1] -= lambda;

    slope_low[0] = num::Float::infinity();
    slope_up[0] = num::Float::neg_infinity();

    // `z` is first set to be the first lower bound.
    z[0] = lower_bound[0];

    for i in 1..width {
        c_low += 1;
        c_up += 1;

        index_low[c_low] = i;
        index_up[c_up] = i;
        slope_low[c_low] = lower_bound[i] - lower_bound[i - 1];

        while (c_low > s_low + 1) && (slope_low[cmp::max(s_low, c_low - 1)] <= slope_low[c_low]) {
            c_low -= 1;
            index_low[c_low] = i;
            if c_low > s_low + 1 {
                slope_low[c_low] = (lower_bound[i] - lower_bound[index_low[c_low - 1]]) /
                                   num::FromPrimitive::from_usize(i - index_low[c_low - 1])
                    .expect("Unable to convert usize to num::FromPrimitive.");
            } else {
                slope_low[c_low] = (lower_bound[i] - z[c]) /
                                   num::FromPrimitive::from_usize(i - index[c])
                    .expect("Unable to convert usize to num::FromPrimitive.");
            }
        }

        slope_up[c_up] = upper_bound[i] - upper_bound[i - 1];
        while (c_up > s_up + 1) && (slope_up[cmp::max(c_up - 1, s_up)] >= slope_up[c_up]) {
            c_up -= 1;
            index_up[c_up] = i;
            if c_up > s_up + 1 {
                slope_up[c_up] = (upper_bound[i] - upper_bound[index_up[c_up - 1]]) /
                                 num::FromPrimitive::from_usize(i - index_up[c_up - 1])
                    .expect("Unable to convert usize to num::FromPrimitive.");
            } else {
                slope_up[c_up] = (upper_bound[i] - z[c]) /
                                 num::FromPrimitive::from_usize(i - index[c])
                    .expect("Unable to convert usize to num::FromPrimitive.");
            }
        }
        while (c_low == s_low + 1) && (c_up > s_up + 1) &&
              (slope_low[c_low] >= slope_up[s_up + 1]) {
            c += 1;
            s_up += 1;
            index[c] = index_up[s_up];
            z[c] = upper_bound[index[c]];
            index_low[s_low] = index[c];
            slope_low[c_low] = (lower_bound[i] - z[c]) /
                               num::FromPrimitive::from_usize(i - index[c])
                .expect("Unable to convert usize to num::FromPrimitive.");;
        }
        while (c_up == s_up + 1) && (c_low > s_low + 1) &&
              (slope_up[c_up] <= slope_low[s_low + 1]) {
            c += 1;
            s_low += 1;
            index[c] = index_low[s_low];
            z[c] = lower_bound[index[c]];
            index_up[s_up] = index[c];
            slope_up[c_up] = (upper_bound[i] - z[c]) /
                             num::FromPrimitive::from_usize(i - index[c])
                .expect("Unable to convert usize to num::FromPrimitive.");;
        }
    }

    for i in 1..(c_low - s_low + 1) {
        index[c + i] = index_low[s_low + i];
        z[c + i] = lower_bound[index[c + i]];
    }
    c += c_low - s_low;

    // Finally, write the denoised output.
    let mut output_index = 0;
    let mut denoised_output;
    i = 1;
    while i <= c {
        denoised_output = (z[i] - z[i - 1]) /
                          num::FromPrimitive::from_usize(index[i] - index[i - 1])
            .expect("Unable to convert usize to num::FromPrimitive.");;
        while output_index < index[i] {
            output[output_index] = denoised_output;
            output_index += 1;
        }
        i += 1;
    }
    output
}


/// Denoises the input values based on a non-iterative algorithm
/// described by Condat L. in 2013 in the paper ["A Direct
/// Algorithm for 1D Total Variation
/// Denoising"](https://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/publis/Condat-fast_TV-SPL-2013.pdf).
///
/// `lambda` closer to `0` means the denoised output will resemble the
/// input more. As `lambda` increases, the denoised output values
/// become closer to the average of the input values.
///
/// # Panics
/// Panics if input vector's length is `0`.
///
/// # Examples
///
/// With `lambda` equal to `0`, the denoised output will be the same as the input:
///
/// ```
/// let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
///
/// let denoised_with_zero_lambda = condat(input, 0.0);
/// assert_eq!(denoised_with_zero_lambda, vec![1.0, 2.0, 3.0, 4.0, 5.0])
/// ```
///
/// With larger `lambda`, the denoised output becomes closer to the average of the input:
///
/// ```
/// let input = vec![1.0, 2.0, 3.0, 4.0, 5.0];
///
/// let denoised_with_larger_lambda = condat(input, 10.0);
/// assert_eq!(denoised_with_larger_lambda, vec![3.0, 3.0, 3.0, 3.0, 3.0]);
/// ```
///
pub fn condat<T>(input: &[T], lambda: T) -> Vec<T>
    where T: num::Num + num::FromPrimitive
    + cmp::PartialOrd + ops::Neg<Output=T> + ops::AddAssign<T> + Copy
{
    assert!(input.len() > 0,
            "Input list should have at least one value.");

    let width = input.len();
    let mut output = Vec::with_capacity(width);

    // `current_input_index` is the location of the element the
    // program is currently inspecting.
    let mut current_input_index = 0;

    // `segment_start` is the index for the beginning of current
    // segment.
    let mut segment_start = 0;

    let twolambda = T::from_u8(2).expect("Unable to transform `2` to T.") * lambda;
    let minlambda = -lambda;

    // `umin` and `umax` are used for keeping track of previous data
    // points we have seen, and will be used to decide whether
    // `segment_lower_bound` and `segment_upper_bound` should be
    // adjusted.
    let mut umin = lambda;
    let mut umax = minlambda;

    // `segment_lower_bound` and `segment_upper_bound` are the
    // Boundaries of the segment's value.
    let mut segment_lower_bound = input[0] - lambda;
    let mut segment_upper_bound = input[0] + lambda;

    // The last position where `umax = -lambda`.
    let mut kplus = 0;
    // The last position where `umin = lambda`.
    let mut kminus = 0;

    loop {
        if current_input_index == (width - 1) {
            // Reached the end of the input. Now process the last
            // set of jumps.
            if umin < num::zero() {
                // Negative jump is necessary as `segment_lower_bound`
                // is too high.
                output.extend(iter::repeat(segment_lower_bound).take(kminus - segment_start + 1));
                segment_start = kminus + 1;
                utils::sync_values(segment_start, &mut [&mut current_input_index, &mut kminus]);
                segment_lower_bound = input[kminus];
                umin = lambda;
                umax = segment_lower_bound + umin - segment_upper_bound;
            } else if umax > num::zero() {
                // If `segment_upper_bound` is too low, jump up.
                output.extend(iter::repeat(segment_upper_bound).take(kplus - segment_start + 1));
                segment_start = kplus + 1;
                utils::sync_values(segment_start, &mut [&mut current_input_index, &mut kplus]);
                segment_upper_bound = input[kplus];
                umax = minlambda;
                umin = segment_upper_bound + umax - segment_lower_bound;
            } else {
                // `segment_lower_bound` and `segment_upper_bound` are
                // not too high or not too low. Adjust the
                // `segment_lower_bound` to reflect the difference
                // between the current input value and value at the
                // beginning of the segment, and extend the output.
                segment_lower_bound +=
                    umin /
                    num::FromPrimitive::from_usize(current_input_index - segment_start + 1)
                        .expect("Unable to convert usize to num::FromPrimitive.");
                output.extend(iter::repeat(segment_lower_bound)
                    .take(current_input_index - segment_start + 1));
                return output;
            }
        } else {
            umin += input[current_input_index + 1] - segment_lower_bound;
            umax += input[current_input_index + 1] - segment_upper_bound;
            if umin < minlambda {
                // If next value (`input[current_input_index + 1]`is
                // much smaller than `segment_lower_bound`, make a
                // negative jump. Next value becomes the
                // `segment_lower_bound`, and `segment_upper_bound` is
                // adjusted accordingly.
                output.extend(iter::repeat(segment_lower_bound).take(kminus - segment_start + 1));
                segment_start = kminus + 1;
                utils::sync_values(segment_start,
                            &mut [&mut current_input_index, &mut kminus, &mut kplus]);
                segment_lower_bound = input[kplus];
                segment_upper_bound = segment_lower_bound + twolambda;
                umin = lambda;
                umax = minlambda;
            } else if umax > lambda {
                // If next value (`input[current_input_index + 1]`is
                // much larger than `segment_upper_bound`, make a
                // negative jump. Next value becomes the
                // `segment_upper_bound`, and `segment_lower_bound` is
                // adjusted accordingly.
                output.extend(iter::repeat(segment_upper_bound).take(kplus - segment_start + 1));
                segment_start = kplus + 1;
                utils::sync_values(segment_start,
                            &mut [&mut current_input_index, &mut kminus, &mut kplus]);
                segment_upper_bound = input[kplus];
                segment_lower_bound = segment_upper_bound - twolambda;
                umin = lambda;
                umax = minlambda;
            } else {
                // `segment_upper_bound` and `segment_lower_bound` are
                // appropriate, and therefore no jump is necessary.
                current_input_index += 1;
                if umin >= lambda {
                    // If `umin` is greater than lambda (threshold),
                    // adjust `segment_lower_bound` to be a little
                    // higher.
                    kminus = current_input_index;
                    segment_lower_bound += (umin - lambda) /
                                           num::FromPrimitive::from_usize(kminus - segment_start +
                                                                          1)
                        .expect("Unable to convert usize to num::FromPrimitive.");
                    umin = lambda;
                }
                if umax <= minlambda {
                    // If `umax` is smaller than -lambda (threshold),
                    // adjust `segment_upper_bound` to be a little
                    // lower.
                    kplus = current_input_index;
                    segment_upper_bound += (umax + lambda) /
                                           num::FromPrimitive::from_usize(kplus - segment_start +
                                                                          1)
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
        let output = tautstring(&input, 0.0);
        assert_eq!(input.len(), output.len());
    }

    #[test]
    fn tautstring_test_zero_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = tautstring(&input, 0.0);
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
        let output = tautstring(&input, 100.0);
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
        let output = tautstring(&input, 5.0);
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
    #[should_panic]
    fn tautstring_test_empty_input() {
        let input = vec![];
        condat(&input, 1.0);
    }

    #[test]
    fn condat_test_input_output_length() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat(&input, 0.0);
        assert_eq!(input.len(), output.len());
    }

    #[test]
    fn condat_test_zero_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat(&input, 0.0);
        let output_expected = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        assert_eq!(output, output_expected);
    }

    #[test]
    fn condat_test_large_lambda() {
        let input = vec![111.0, 422.1, 145.2, 248.2, 871.4, 675.2, 436.2, 310.1];
        let output = condat(&input, 700.0);
        // The expected output is taken from the Laurent Condat's C
        // implementation.
        let output_expected = vec![402.425049, 402.425049, 402.425049, 402.425049, 402.425049,
                                   402.425049, 402.425049, 402.425049];
        println!("{:?}", output);
        for i in 0..input.len() {
            let output_data = output[i] as f64;
            let expected_data = output_expected[i] as f64;
            assert!((output_data - expected_data).abs() <= 0.0001);
        }
    }

    #[test]
    fn condat_test_moderate_lambda() {
        let input = vec![1.0, 2.1, 5.2, 8.2, 1.4, 5.2, 6.2, 10.1];
        let output = condat(&input, 3.0);
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
        condat(&input, -1.0);
    }

    #[test]
    #[should_panic]
    fn condat_test_empty_input() {
        let input = vec![];
        condat(&input, 1.0);
    }
}
