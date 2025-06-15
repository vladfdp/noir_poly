use crate::domain::Domain;

pub fn add_arrays<let N: u32, let M: u32>(a: [Field; N], b: [Field; M]) -> [Field; N] {
    assert(N >= M, "first polynomial should be bigger than second polynomial");
    let mut result = a;
    for i in 0..M {
        result[i] += b[i];
    }
    result
}

pub fn sub_arrays<let N: u32, let M: u32>(a: [Field; N], b: [Field; M]) -> [Field; N] {
    assert(N >= M, "first polynomial should be bigger than second polynomial");
    let mut result = a;
    for i in 0..M {
        result[i] -= b[i];
    }
    result
}

pub fn neg_array<let N: u32>(a: [Field; N]) -> [Field; N] {
    let mut result = [0; N];
    for i in 0..N {
        result[i] = -a[i];
    }
    result
}

pub fn scalar_multiply_array<let N: u32>(a: [Field; N], scalar: Field) -> [Field; N] {
    let mut result = [0; N];
    for i in 0..N {
        result[i] = a[i] * scalar;
    }
    result
}

pub fn pointwise_multiply_arrays<let N: u32>(a: [Field; N], b: [Field; N]) -> [Field; N] {
    let mut result = [0; N];
    for i in 0..N {
        result[i] = a[i] * b[i];
    }
    result
}

pub fn get_last_non_zero_index<let N: u32>(coefficients: [Field; N]) -> u32 {
    let mut index = 0;
    for i in 0..N {
        if coefficients[i] != 0 {
            index = i;
        }
    }
    index
}

// Helper function to reverse the bits of the index
fn reverse_bits<let LOG_N: u32>(index: u32) -> u32 {
    let mut reversed_index: u32 = 0;
    let mut temp_index = index;
    for _i in 0..LOG_N {
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

comptime fn log2(N: u32) -> u32 {
    let mut count = 0;
    let mut n = N;
    while n > 1 {
        n >>= 1;
        count += 1;
    }
    count
}

// Reorders the elements of an array according to the bit-reversal permutation of indices.
pub fn bit_reverse_array<let N: u32, let LOG_N: u32>(a: [Field; N]) -> [Field; N] {
    let mut result = [0; N];

    for i in 0..N {
        let reversed_index = reverse_bits::<LOG_N>(i);
        result[reversed_index] = a[i];
    }
    result
}

// Cooley-Tukey FFT butterfly operation
pub fn CT_butterfly<let SIZE: u32, let LOG_SIZE: u32, let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(coeffs: [Field; SIZE], domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> [Field; SIZE] {
    let step = DOMAIN_SIZE/SIZE;

    let mut new_coeffs = coeffs;
    for i in 1..(LOG_SIZE + 1) {
        let block_size = 1 << (i as u8);
        let num_blocks = SIZE / block_size;
        let twiddle_factor = domain.get_root(num_blocks * step);

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

pub fn inverse_CT_butterfly<let SIZE: u32, let LOG_SIZE: u32, let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(coeffs: [Field; SIZE], domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> [Field; SIZE] {
    let step = DOMAIN_SIZE/SIZE;

    let mut new_coeffs = coeffs;
    for i in 1..(LOG_SIZE + 1) {
        let block_size = 1 << (i as u8);
        let num_blocks = SIZE / block_size;
        let twiddle_factor = 1 / domain.get_root(num_blocks * step);

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
    scalar_multiply_array::<SIZE>(new_coeffs, 1 / (SIZE as Field))
}

// ==========================================
// Tests for Array Utilities
// ==========================================

#[test]
fn test_add_arrays() {
    let a = [1, 2, 3, 4];
    let b = [5, 6];
    let result = add_arrays::<4, 2>(a, b);
    assert(result[0] == 6);
    assert(result[1] == 8);
    assert(result[2] == 3);
    assert(result[3] == 4);

    // Edge case: same size arrays
    let c = [1, 2];
    let d = [3, 4];
    let result2 = add_arrays::<2, 2>(c, d);
    assert(result2[0] == 4);
    assert(result2[1] == 6);
}

#[test]
fn test_sub_arrays() {
    let a = [10, 20, 30, 40];
    let b = [5, 6];
    let result = sub_arrays::<4, 2>(a, b);
    assert(result[0] == 5);
    assert(result[1] == 14);
    assert(result[2] == 30);
    assert(result[3] == 40);

    // Edge case: same size arrays
    let c = [10, 15];
    let d = [3, 5];
    let result2 = sub_arrays::<2, 2>(c, d);
    assert(result2[0] == 7);
    assert(result2[1] == 10);
}

#[test]
fn test_neg_array() {
    let a = [1, -2, 0, 5];
    let result = neg_array::<4>(a);
    assert(result[0] == -1);
    assert(result[1] == 2);
    assert(result[2] == 0);
    assert(result[3] == -5);
}

#[test]
fn test_scalar_multiply_array() {
    let a = [1, 2, 3];
    let scalar = 5;
    let result = scalar_multiply_array::<3>(a, scalar);
    assert(result[0] == 5);
    assert(result[1] == 10);
    assert(result[2] == 15);

    // Edge case: multiply by zero
    let result_zero = scalar_multiply_array::<3>(a, 0);
    assert(result_zero[0] == 0);
    assert(result_zero[1] == 0);
    assert(result_zero[2] == 0);
}

#[test]
fn test_pointwise_multiply_arrays() {
    let a = [1, 2, 3, 4];
    let b = [2, 3, 0, 5];
    let result = pointwise_multiply_arrays::<4>(a, b);
    assert(result[0] == 2);
    assert(result[1] == 6);
    assert(result[2] == 0);
    assert(result[3] == 20);
}

#[test]
fn test_get_last_non_zero_index() {
    let a = [1, 2, 3, 0];
    let index = get_last_non_zero_index::<4>(a);
    assert(index == 2);

    // Edge case: all zeros
    let b = [0, 0, 0, 0];
    let index_zero = get_last_non_zero_index::<4>(b);
    assert(index_zero == 0);

    // Edge case: last element non-zero
    let c = [0, 0, 0, 5];
    let index_last = get_last_non_zero_index::<4>(c);
    assert(index_last == 3);
}

#[test]
fn test_bit_reverse_array() {
    // Test with size 4 (LOG_N = 2)
    let a = [1, 2, 3, 4];
    let result = bit_reverse_array::<4, 2>(a);
    // Bit reversal: 0->0, 1->2, 2->1, 3->3
    assert(result[0] == 1); // a[0] -> result[0]
    assert(result[1] == 3); // a[2] -> result[1]
    assert(result[2] == 2); // a[1] -> result[2]
    assert(result[3] == 4); // a[3] -> result[3]

    // Test with size 2 (LOG_N = 1)
    let b = [5, 6];
    let result2 = bit_reverse_array::<2, 1>(b);
    // Bit reversal: 0->0, 1->1
    assert(result2[0] == 5);
    assert(result2[1] == 6);
}
