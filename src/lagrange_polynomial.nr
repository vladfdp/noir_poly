use crate::domain::Domain;
use crate::polynomial::Polynomial;
use crate::array_utils::{add_arrays, sub_arrays, neg_array, scalar_multiply_array, pointwise_multiply_arrays};
use crate::array_utils::{bit_reverse_array, inverse_CT_butterfly};

pub struct LagrangePolynomial<let SIZE: u32, let LOG_SIZE: u32> {
    evals: [Field; SIZE],
}


impl<let SIZE: u32, let LOG_SIZE: u32> Eq for LagrangePolynomial<SIZE, LOG_SIZE> {
    fn eq(self, other: Self) -> bool {
        (self.evals == other.evals)
    }
}


impl<let SIZE: u32, let LOG_SIZE: u32> LagrangePolynomial<SIZE, LOG_SIZE> {
    pub fn new(evals: [Field; SIZE]) -> Self {
        assert(SIZE == 1 << (LOG_SIZE as u8));
        Self { evals }
    }

    pub fn new_from_slice(evals: [Field]) -> Self {
        assert(SIZE == 1 << (LOG_SIZE as u8));
        
        let mut evals_array = [0; SIZE];

        for i in 0..SIZE {
            if i < evals.len() {
                evals_array[i] = evals[i];
            }
        }

        Self::new(evals_array)
    }

    pub fn get_evals(self) -> [Field; SIZE] {
        self.evals
    }

    pub fn one() -> Self {
        let mut evals = [1; SIZE];
        Self::new(evals)
    }

    pub fn zero() -> Self {
        let mut evals = [0; SIZE];
        Self::new(evals)
    }

    pub fn is_not_zero(self) -> bool {
        self.evals.any(|eval| eval != 0)
    }


    pub fn add(self, other: Self) -> Self {
        let mut evals = add_arrays::<SIZE, SIZE>(self.evals, other.evals);
        Self::new(evals)
    }

    pub fn sub(self, other: Self) -> Self {
        let mut evals = sub_arrays::<SIZE, SIZE>(self.evals, other.evals);
        Self::new(evals)
    }

    pub fn neg(self) -> Self {
        Self { evals: neg_array::<SIZE>(self.evals) }
    }

    pub fn mul_scalar(self, scalar: Field) -> Self {
        Self { evals: scalar_multiply_array::<SIZE>(self.evals, scalar) }
    }

    pub fn mul(self, other: Self) -> Self {
        let mut evals = pointwise_multiply_arrays::<SIZE>(self.evals, other.evals);
        Self { evals }
    }

    pub fn intt<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(self, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> Polynomial<SIZE, LOG_SIZE> {
        let br_evals = bit_reverse_array::<SIZE, LOG_SIZE>(self.evals);
        let coeffs = inverse_CT_butterfly::<SIZE, LOG_SIZE, DOMAIN_SIZE, DOMAIN_LOG_SIZE>(br_evals, domain);
        Polynomial::new(coeffs)
    }
    
}


// =========================================
// Tests for LagrangePolynomial Implementation
// =========================================



#[test]
fn test_lagrange_add_sub_neg() {
    // Create LagrangePolynomials with capacity 4 (next power of 2 after length 3)
    let poly1: LagrangePolynomial<4, 2> = LagrangePolynomial::<4, 2>::new_from_slice([1, 2, 3]);
    let poly2: LagrangePolynomial<4, 2> = LagrangePolynomial::<4, 2>::new_from_slice([4, 5, 6]);

    // Test add
    let add_poly = poly1.add(poly2);
    let expected_add_evals_slice = [5, 7, 9];
    let add_evals = add_poly.get_evals();
    for i in 0..expected_add_evals_slice.len() {
        assert(add_evals[i] == expected_add_evals_slice[i]);
    }

    // Test sub
    let sub_poly = poly1.sub(poly2);
    let expected_sub_evals_slice = [-3, -3, -3];
    let sub_evals = sub_poly.get_evals();
    for i in 0..expected_sub_evals_slice.len() {
        assert(sub_evals[i] == expected_sub_evals_slice[i]);
    }
    
    // Test neg
    let neg_poly = poly1.neg();
    let expected_neg_evals_slice = [-1, -2, -3];
    let neg_evals = neg_poly.get_evals();
    for i in 0..expected_neg_evals_slice.len() {
        assert(neg_evals[i] == expected_neg_evals_slice[i]);
    }
}


#[test]
fn test_lagrange_mul_scalar() {
    // Create LagrangePolynomial with capacity 4 (next power of 2 after length 3)
    let poly: LagrangePolynomial<4, 2> = LagrangePolynomial::<4, 2>::new_from_slice([1, 2, 3]);
    let scalar = 5;

    let mul_poly = poly.mul_scalar(scalar);
    let expected_evals_slice = [5, 10, 15];
    let mul_evals = mul_poly.get_evals();
     for i in 0..expected_evals_slice.len() {
        assert(mul_evals[i] == expected_evals_slice[i]);
    }
}

#[test]
fn test_lagrange_mul() {
    // Create LagrangePolynomials with capacity 2 (power of 2 equal to length 2)
    let poly1: LagrangePolynomial<2, 1> = LagrangePolynomial::<2, 1>::new_from_slice([1, 2]); 
    let poly2: LagrangePolynomial<2, 1> = LagrangePolynomial::<2, 1>::new_from_slice([3, 4]);

    // Test mul
    let mul_poly = poly1.mul(poly2); // Represents P(x) = L(x) * M(x)
    let expected_evals_slice = [3, 8]; // P(w^0)=1*3=3, P(w^1)=2*4=8
    let mul_evals = mul_poly.get_evals();
    
    for i in 0..expected_evals_slice.len() {
        assert(mul_evals[i] == expected_evals_slice[i]);
    }
}

#[test]
fn test_lagrange_intt() {
    // Define a simple polynomial: P(x) = 1 + 2x + 3x^2
    // Capacity 4 (next power of 2 after length 3)
    let coefficients_slice = [1, 2, 3];
    let poly: Polynomial<4, 2> = Polynomial::<4, 2>::new_from_slice(coefficients_slice);
    assert(poly.get_degree() == 2);

    // Domain size must be >= poly capacity and a power of 2
    // We need a domain of size 4
    let domain: Domain<4, 2> = comptime { Domain::<4, 2>::new() }; 

    // 1. Perform NTT to get Lagrange evaluations
    let lagrange_poly: LagrangePolynomial<4, 2> = poly.ntt(domain);
    
    // 2. Perform INTT to get back the coefficients
    let recovered_poly: Polynomial<4, 2> = lagrange_poly.intt(domain);
    
    // Ensure the degree is correct
    assert(recovered_poly.get_degree() == 2);

    // 3. Check if the recovered polynomial matches the original
    let recovered_coeffs = recovered_poly.get_coefficients();
    let original_coeffs = poly.get_coefficients();

    // Compare all coefficients
    for i in 0..4 { 
         assert(original_coeffs[i] == recovered_coeffs[i]);
    }

    // Test with zero polynomial
    let zero_p: Polynomial<4, 2> = Polynomial::<4, 2>::zero();
    let zero_l: LagrangePolynomial<4, 2> = zero_p.ntt(domain);
    let recovered_zero: Polynomial<4, 2> = zero_l.intt(domain);
    assert(recovered_zero.get_degree() == 0); // Degree should be 0
    let zero_coeffs = zero_p.get_coefficients();
    let recovered_zero_coeffs = recovered_zero.get_coefficients();
    for i in 0..4 {
        assert(zero_coeffs[i] == recovered_zero_coeffs[i]);
    }

    // Test with one polynomial P(x) = 1
    let one_p: Polynomial<4, 2> = Polynomial::<4, 2>::one();
    let one_l: LagrangePolynomial<4, 2> = one_p.ntt(domain);
    let recovered_one: Polynomial<4, 2> = one_l.intt(domain);
    assert(recovered_one.get_degree() == 0); // Degree should be 0
    let one_coeffs = one_p.get_coefficients();
    let recovered_one_coeffs = recovered_one.get_coefficients();
    for i in 0..4 {
        assert(one_coeffs[i] == recovered_one_coeffs[i]);
    }
}
