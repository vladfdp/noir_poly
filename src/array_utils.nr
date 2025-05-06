// Define MAX_DEGREE or import it if it's defined elsewhere (e.g., in a constants file)
// Assuming it's defined in src/constants.nr or similar
use crate::MAX_DEGREE::{MAX_DEGREE, MAX_DEGREE_LOG, NORMALIZER};
use crate::domain::Domain;

pub fn add_arrays(a: [Field; MAX_DEGREE + 1], b: [Field; MAX_DEGREE + 1]) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    for i in 0..MAX_DEGREE + 1 {
        result[i] = a[i] + b[i];
    }
    result
}

pub fn sub_arrays(a: [Field; MAX_DEGREE + 1], b: [Field; MAX_DEGREE + 1]) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    for i in 0..MAX_DEGREE + 1 {
        result[i] = a[i] - b[i];
    }
    result
}

pub fn neg_array(a: [Field; MAX_DEGREE + 1]) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    for i in 0..MAX_DEGREE + 1 {
        result[i] = -a[i];
    }
    result
}

pub fn scalar_multiply_array(a: [Field; MAX_DEGREE + 1], scalar: Field) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    for i in 0..MAX_DEGREE + 1 {
        result[i] = a[i] * scalar;
    }
    result
}

pub fn pointwise_multiply_arrays(a: [Field; MAX_DEGREE + 1], b: [Field; MAX_DEGREE + 1]) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    for i in 0..MAX_DEGREE + 1 {
        result[i] = a[i] * b[i];
    }
    result
}

// Helper function to reverse the bits of the index
fn reverse_bits(index: u32) -> u32 {
    let mut reversed_index: u32 = 0;
    let mut temp_index = index;
    for _i in 0..MAX_DEGREE_LOG {
        // Prepend the LSB to the reversed_index (left shift)
        reversed_index <<= 1;
        // Check if the LSB of temp_index is 1
        if (temp_index & 1) == 1 {
            // If LSB is 1, set the LSB of reversed_index to 1 (using OR)
            reversed_index |= 1;
        }
        // Drop the LSB from temp_index (right shift)
        temp_index >>= 1;
    }
    reversed_index
}

// Reorders the elements of an array according to the bit-reversal permutation of indices.
pub fn bit_reverse_array(a: [Field; MAX_DEGREE + 1]) -> [Field; MAX_DEGREE + 1] {
    let mut result = [0; MAX_DEGREE + 1];
    let n = MAX_DEGREE + 1;

    for i in 0..n {
        let reversed_index = reverse_bits(i);
        result[reversed_index] = a[i];
    }
    result
}

// Cooley-Tukey FFT butterfly operation
pub fn CT_butterfly(coeffs: [Field; MAX_DEGREE + 1], domain: Domain) -> [Field; MAX_DEGREE + 1] {
    let N = MAX_DEGREE + 1;
    let mut new_coeffs = coeffs;
    for i in 1..(MAX_DEGREE_LOG + 1) {
        let block_size = 1 << i;
        let num_blocks = N / block_size;
        let twiddle_factor = domain.get_root(num_blocks);

        for j in 0..num_blocks {
            let mut wp = 1;
            for k in 0..(block_size / 2) {
                let index = j * block_size + k;
                let index_plus_half = index + block_size / 2;
                let u = new_coeffs[index];
                let v = new_coeffs[index_plus_half] * wp;
                new_coeffs[index] = u + v;
                new_coeffs[index_plus_half] = u - v;
                wp = wp * twiddle_factor;                
            }   
        }
    }
    new_coeffs
}

pub fn inverse_CT_butterfly(coeffs: [Field; MAX_DEGREE + 1], domain: Domain) -> [Field; MAX_DEGREE + 1] {
    let N = MAX_DEGREE + 1;
    let mut new_coeffs = coeffs;
    for i in 1..(MAX_DEGREE_LOG + 1) {
        let block_size = 1 << i;
        let num_blocks = N / block_size;
        let twiddle_factor = 1 / domain.get_root(num_blocks);

        for j in 0..num_blocks {
            let mut wp = 1;
            for k in 0..(block_size / 2) {
                let index = j * block_size + k;
                let index_plus_half = index + block_size / 2;
                let u = new_coeffs[index];
                let v = new_coeffs[index_plus_half] * wp;
                new_coeffs[index] = u + v;
                new_coeffs[index_plus_half] = u - v;
                wp = wp * twiddle_factor;                
            }   
        }
    }
    scalar_multiply_array(new_coeffs, NORMALIZER)
}
