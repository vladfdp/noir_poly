// Define MAX_DEGREE or import it if it's defined elsewhere (e.g., in a constants file)
// Assuming it's defined in src/constants.nr or similar
use crate::MAX_DEGREE::{MAX_DEGREE, MAX_DEGREE_LOG};
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
pub fn bit_reverse_array<let N: u32>(a: [Field; N]) -> [Field; N] {
    comptime let LOG_N = 3; //comptime {log2(N)};        // TODO: remove this________________________________________________________

    let mut result = [0; N];

    for i in 0..N {
        let reversed_index = reverse_bits::<LOG_N>(i);
        result[reversed_index] = a[i];
    }
    result
}

// Cooley-Tukey FFT butterfly operation
pub fn CT_butterfly<let N: u32, let M: u32>(coeffs: [Field; N], domain: Domain<M>) -> [Field; N] {
    let LOG_N = comptime {log2(N)};
    let LOG_M = comptime {log2(M)};
    let step = comptime {1 << (LOG_M - LOG_N)};

    let mut new_coeffs = coeffs;
    for i in 1..(LOG_N + 1) {
        let block_size = 1 << i;
        let num_blocks = N / block_size;
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

pub fn inverse_CT_butterfly<let N: u32, let M: u32>(coeffs: [Field; N], domain: Domain<M>) -> [Field; N] {
    let LOG_N = comptime {log2(N)};
    let LOG_M = comptime {log2(M)};
    let step = comptime {1 << (LOG_M - LOG_N)};

    let mut new_coeffs = coeffs;
    for i in 1..(LOG_N + 1) {
        let block_size = 1 << i;
        let num_blocks = N / block_size;
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
    scalar_multiply_array::<N>(new_coeffs, 1 / (N as Field))
}
