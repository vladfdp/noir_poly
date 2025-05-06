use crate::MAX_DEGREE::MAX_DEGREE;
use crate::MAX_DEGREE::MAX_DEGREE_LOG;
use crate::MAX_DEGREE::NORMALIZER;
use crate::domain::Domain;
use crate::polynomial::Polynomial;
use crate::array_utils::{add_arrays, sub_arrays, neg_array, scalar_multiply_array, pointwise_multiply_arrays, bit_reverse_array, inverse_CT_butterfly};

pub struct LagrangePolynomial {
    evals: [Field; MAX_DEGREE + 1],
    degree: u32,
}


impl Eq for LagrangePolynomial {
    fn eq(self, other: Self) -> bool {
        (self.evals == other.evals) & (self.degree == other.degree)
    }
}


impl LagrangePolynomial {
    pub fn new(evals: [Field; MAX_DEGREE + 1], degree: u32) -> Self {
        assert(degree <= MAX_DEGREE);
        Self { evals, degree }
    }

    pub fn new_from_slice(evals: [Field]) -> Self {

        let mut evals_array = [0; MAX_DEGREE + 1];
        let degree = evals.len() - 1;
        assert(degree <= MAX_DEGREE, "Degree is too large");

        for i in 0..MAX_DEGREE + 1 {
            if i < degree + 1 {
                evals_array[i] = evals[i];
            }
        }

        Self::new(evals_array, degree)
    }

    pub fn get_degree(self) -> u32 {
        self.degree
    }

    pub fn get_evals(self) -> [Field; MAX_DEGREE + 1] {
        self.evals
    }

    pub fn one() -> Self {
        let mut evals = [1; MAX_DEGREE + 1];
        Self::new(evals, 0)
    }

    pub fn zero() -> Self {
        let mut evals = [0; MAX_DEGREE + 1];
        Self::new(evals, 0)
    }

    pub fn is_not_zero(self) -> bool {
        self.evals.any(|eval| eval != 0)
    }


    pub fn add(self, other: Self) -> Self {
        let mut evals = add_arrays(self.evals, other.evals);
        Self::new(evals, self.degree)
    }

    pub fn sub(self, other: Self) -> Self {
        let mut evals = sub_arrays(self.evals, other.evals);
        Self::new(evals, self.degree)
    }

    pub fn neg(self) -> Self {
        Self { evals: neg_array(self.evals), degree: self.degree }
    }

    pub fn mul_scalar(self, scalar: Field) -> Self {
        Self { evals: scalar_multiply_array(self.evals, scalar), degree: self.degree }
    }

    pub fn mul(self, other: Self) -> Self {
        let mut evals = pointwise_multiply_arrays(self.evals, other.evals);
        let degree = self.degree + other.degree;
        Self { evals, degree }
    }

    pub fn intt(self, domain: Domain) -> Polynomial {
        let mut evals = bit_reverse_array(self.evals);
        evals = inverse_CT_butterfly(evals, domain);
        let degree = self.degree;
        Polynomial::new(evals, degree)
    }
    
}


// =========================================
// Tests for LagrangePolynomial Implementation
// =========================================



#[test]
fn test_lagrange_add_sub_neg() {
    let poly1 = LagrangePolynomial::new_from_slice([1, 2, 3]);
    let poly2 = LagrangePolynomial::new_from_slice([4, 5, 6]);
    let expected_degree = 2; // Based on slices [1,2,3] and [4,5,6]

    assert(poly1.get_degree() == expected_degree);
    assert(poly2.get_degree() == expected_degree);

    // Test add
    let add_poly = poly1.add(poly2);
    let expected_add_evals_slice = [5, 7, 9];
    let add_evals = add_poly.get_evals();
    // Note: Lagrange add/sub don't change the degree field currently.
    // The degree field likely represents the max degree seen, not recalculated.
    assert(add_poly.get_degree() == expected_degree);
    for i in 0..expected_add_evals_slice.len() {
        assert(add_evals[i] == expected_add_evals_slice[i]);
    }
    // Check padding is zero
    for i in expected_add_evals_slice.len()..MAX_DEGREE+1 {
         assert(add_evals[i] == 0);
    }

    // Test sub
    let sub_poly = poly1.sub(poly2);
    let expected_sub_evals_slice = [-3, -3, -3];
    let sub_evals = sub_poly.get_evals();
    assert(sub_poly.get_degree() == expected_degree);
    for i in 0..expected_sub_evals_slice.len() {
        assert(sub_evals[i] == expected_sub_evals_slice[i]);
    }
     // Check padding is zero
    for i in expected_sub_evals_slice.len()..MAX_DEGREE+1 {
         assert(sub_evals[i] == 0);
    }

    // Test neg
    let neg_poly = poly1.neg();
    let expected_neg_evals_slice = [-1, -2, -3];
    let neg_evals = neg_poly.get_evals();
    assert(neg_poly.get_degree() == expected_degree);
    for i in 0..expected_neg_evals_slice.len() {
        assert(neg_evals[i] == expected_neg_evals_slice[i]);
    }
    // Check padding is zero
    for i in expected_neg_evals_slice.len()..MAX_DEGREE+1 {
         assert(neg_evals[i] == 0);
    }
}


#[test]
fn test_lagrange_mul_scalar() {
    let poly = LagrangePolynomial::new_from_slice([1, 2, 3]);
    let expected_degree = 2;
    let scalar = 5;

    assert(poly.get_degree() == expected_degree);

    let mul_poly = poly.mul_scalar(scalar);
    let expected_evals_slice = [5, 10, 15];
    let mul_evals = mul_poly.get_evals();
    assert(mul_poly.get_degree() == expected_degree);
     for i in 0..expected_evals_slice.len() {
        assert(mul_evals[i] == expected_evals_slice[i]);
    }
    // Check padding is zero
    for i in expected_evals_slice.len()..MAX_DEGREE+1 {
         assert(mul_evals[i] == 0);
    }
}

#[test]
fn test_lagrange_mul() {
    let poly1 = LagrangePolynomial::new_from_slice([1, 2]); // Degree 1
    let poly2 = LagrangePolynomial::new_from_slice([3, 4]); // Degree 1

    assert(poly1.get_degree() == 1);
    assert(poly2.get_degree() == 1);

    // Test mul
    let mul_poly = poly1.mul(poly2); // Represents P(x) = L(x) * M(x)
    let expected_evals_slice = [3, 8]; // P(w^0)=1*3=3, P(w^1)=2*4=8
    let expected_degree = 2; // Degree should be poly1.degree + poly2.degree
    let mul_evals = mul_poly.get_evals();
    
    assert(mul_poly.get_degree() == expected_degree);
    for i in 0..expected_evals_slice.len() {
        assert(mul_evals[i] == expected_evals_slice[i]);
    }
     // Check padding is zero
    for i in expected_evals_slice.len()..MAX_DEGREE+1 {
         assert(mul_evals[i] == 0);
    }
}

#[test]
fn test_lagrange_intt() {
    // Define a simple polynomial: P(x) = 1 + 2x + 3x^2
    let coefficients_slice = [1, 2, 3];
    let poly = Polynomial::new_from_slice(coefficients_slice);
    let degree = 2; // Degree from slice
    assert(poly.get_degree() == degree);

    // Domain size must be >= poly.degree + 1 and a power of 2.
    let domain_size = 1 << MAX_DEGREE_LOG; 
    assert(degree < domain_size); // Ensure domain is large enough (degree strictly less than size)

    // Domain::new must be called at compile time.
    let domain = comptime { Domain::new() }; 

    // 1. Perform NTT to get Lagrange evaluations
    let lagrange_poly = poly.ntt(domain);
    
    // Ensure the degree passed to LagrangePolynomial is correct
    assert(lagrange_poly.get_degree() == degree); 

    // 2. Perform INTT to get back the coefficients
    let recovered_poly = lagrange_poly.intt(domain);
    
    // Ensure the degree passed to Polynomial is correct
    assert(recovered_poly.get_degree() == degree);

    // 3. Check if the recovered polynomial matches the original
    let recovered_coeffs = recovered_poly.get_coefficients();
    let original_coeffs = poly.get_coefficients(); // This already has padding from new_from_slice

    // Assuming inverse_CT_butterfly includes scaling by 1/N.
    // Compare all coefficients in the padded array.
    for i in 0..MAX_DEGREE + 1 { 
         assert(original_coeffs[i] == recovered_coeffs[i]);
    }

     // Test with zero polynomial
     let zero_p = Polynomial::zero();
     let zero_l = zero_p.ntt(domain);
     assert(zero_l.get_degree() == 0);
     let recovered_zero = zero_l.intt(domain);
     assert(recovered_zero.get_degree() == 0); // Degree should be 0
     let zero_coeffs = zero_p.get_coefficients();
     let recovered_zero_coeffs = recovered_zero.get_coefficients();
     for i in 0..MAX_DEGREE + 1 {
         assert(zero_coeffs[i] == recovered_zero_coeffs[i]);
     }

     // Test with one polynomial P(x) = 1
     let one_p = Polynomial::one(); // degree 0, coeffs [1, 0, ...]
     let one_l = one_p.ntt(domain);
     assert(one_l.get_degree() == 0);
     let recovered_one = one_l.intt(domain);
     assert(recovered_one.get_degree() == 0); // Degree should be 0
     let one_coeffs = one_p.get_coefficients();
     let recovered_one_coeffs = recovered_one.get_coefficients();
     for i in 0..MAX_DEGREE + 1 {
         assert(one_coeffs[i] == recovered_one_coeffs[i]);
     }
}
