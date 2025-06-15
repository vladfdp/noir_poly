use crate::domain::Domain;
use crate::array_utils::{add_arrays, sub_arrays, neg_array, scalar_multiply_array, pointwise_multiply_arrays, get_last_non_zero_index};
use crate::array_utils::{bit_reverse_array, CT_butterfly};
use crate::lagrange_polynomial::LagrangePolynomial;

pub struct Polynomial<let SIZE: u32, let LOG_SIZE: u32> {
    coefficients: [Field; SIZE],
}


impl<let SIZE: u32, let LOG_SIZE: u32> Eq for Polynomial<SIZE, LOG_SIZE> {
    fn eq(self, other: Self) -> bool {
        (self.coefficients == other.coefficients)
    }
}


impl<let SIZE: u32, let LOG_SIZE: u32> Polynomial<SIZE, LOG_SIZE> {

    pub fn new(coefficients: [Field; SIZE]) -> Self {
        assert(SIZE == 1 << (LOG_SIZE as u8));
        Self { coefficients }
    }

    // Creates a polynomial from a variable-length slice of coefficients
    // Coefficients are ordered from constant term to highest degree term
    pub fn new_from_slice(coefficients: [Field]) -> Self {
        assert(SIZE == 1 << (LOG_SIZE as u8));

        let mut coeff_array = [0; SIZE];

        for i in 0..SIZE {
            if i < coefficients.len() {
                coeff_array[i] = coefficients[i];
            }
        }

        Self::new(coeff_array)
    }


    pub fn get_degree(self) -> u32 {
        get_last_non_zero_index(self.coefficients)
    }

    pub fn get_coefficients(self) -> [Field; SIZE] {
        self.coefficients
    }

    pub fn one() -> Self {
        let mut coeff_array = [0; SIZE];
        coeff_array[0] = 1;
        Self::new(coeff_array)
    }

    pub fn zero() -> Self {
        let mut coeff_array = [0; SIZE];
        Self::new(coeff_array)
    }

    pub fn is_not_zero(self) -> bool {
        self.coefficients.any(|coeff| coeff != 0)
    }

    // Expands a polynomial to a larger capacity by adding zeros
    pub fn expand<let M: u32, let LOG_M: u32>(self) -> Polynomial<M, LOG_M> {
        assert(SIZE <= M);
        let mut coeff_array = [0; M];
        for i in 0..SIZE {
            coeff_array[i] = self.coefficients[i];
        }
        Polynomial::<M, LOG_M>::new(coeff_array)
    }

    // Reduces a polynomial to a smaller capacity by removing trailing zeros
    pub fn reduce<let M: u32, let LOG_M: u32>(self) -> Polynomial<M, LOG_M> {
        assert(SIZE >= M);
        let mut coeff_array = [0; M];
        for i in 0..M {
            coeff_array[i] = self.coefficients[i];
        }
        for i in M..SIZE {
            assert(self.coefficients[i] == 0);
        }
        Polynomial::<M, LOG_M>::new(coeff_array)
    }
    

    pub fn add<let OTHER_SIZE: u32, let OTHER_LOG_SIZE: u32>(self, other: Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>) -> Self {
        let mut coeff_array = add_arrays::<SIZE, OTHER_SIZE>(self.coefficients, other.coefficients);
        Self { coefficients: coeff_array }
    }

    pub fn sub<let OTHER_SIZE: u32, let OTHER_LOG_SIZE: u32>(self, other: Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>) -> Self {
        let mut coeff_array = sub_arrays::<SIZE, OTHER_SIZE>(self.coefficients, other.coefficients);
        Self { coefficients: coeff_array }
    }

    pub fn neg(self) -> Self {
        Self { coefficients: neg_array(self.coefficients) }
    }

    // Creates a polynomial (x - a) that evaluates to 0 at point x = a
    pub fn vanishing_polynomial(a: Field) -> Self {
        let mut coeff_array = [0; SIZE];
        coeff_array[0] = -a;
        coeff_array[1] = 1;
        Self::new(coeff_array)
    }

    pub fn mul_scalar(self, scalar: Field) -> Self {
        Self { coefficients: scalar_multiply_array(self.coefficients, scalar) }
    }


    pub fn evaluate(self, x: Field) -> Field {
        let mut result: Field = 0;
        let mut x_power = 1;
        for i in 0..SIZE {
            result += self.coefficients[i] * x_power;
            x_power *= x;
        }
        result
    }

    pub fn mul_by_x(self) -> Self {
        assert(self.coefficients[SIZE - 1] == 0);
        let mut coeff_array = [0; SIZE];
        for i in 1..SIZE {
            coeff_array[i] = self.coefficients[i - 1];
        }
        Self { coefficients: coeff_array }
    }

    pub fn ntt<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(self, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> LagrangePolynomial<SIZE, LOG_SIZE> {
        let mut new_coeffs = bit_reverse_array::<SIZE, LOG_SIZE>(self.coefficients);
        new_coeffs = CT_butterfly::<SIZE, LOG_SIZE, DOMAIN_SIZE, DOMAIN_LOG_SIZE>(new_coeffs, domain);
        LagrangePolynomial::<SIZE, LOG_SIZE>::new(new_coeffs)
    }

    // Polynomial multiplication using ntt and intt
    // Consider using unconstrained operations with verification, if hash is fast, we can check if f(x) * g(x) = h(x) with x pseudo-randomly chosen
    pub fn mul<let OTHER_SIZE: u32, let OTHER_LOG_SIZE: u32, let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(self, other: Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> Polynomial<SIZE + SIZE, LOG_SIZE + 1> {
        assert(SIZE <= DOMAIN_SIZE);
        assert(SIZE >= OTHER_SIZE);
        
        let expanded_self: Polynomial<SIZE + SIZE, LOG_SIZE + 1> = self.expand();
        let expanded_other: Polynomial<SIZE + SIZE, LOG_SIZE + 1> = other.expand();
        let lagrange_self = expanded_self.ntt(domain);
        let lagrange_other = expanded_other.ntt(domain);
        let lagrange_result = lagrange_self.mul(lagrange_other);
        lagrange_result.intt(domain)
    }


    // Polynomial multiplication mod x^N - 1
    pub fn mul_mod<let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(self, other: Self, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> Self {
        assert(SIZE <= DOMAIN_SIZE);
        let lagrange_self = self.ntt(domain);
        let lagrange_other = other.ntt(domain);
        let lagrange_result = lagrange_self.mul(lagrange_other);
        lagrange_result.intt(domain)
    }



    // Naive polynomial multiplication using coefficient-wise operations
    // pub fn naive_mul(self, other: Self) -> Self {
    //     assert(self.degree + other.degree <= MAX_DEGREE);
    //     let degree = self.degree + other.degree;
    //     let mut coeff_array = [0; MAX_DEGREE + 1];
    //     let mut other_times_x = other;
    //     for i in 0..MAX_DEGREE + 1 {
    //         if i <= self.degree {
    //             coeff_array = add_arrays(coeff_array, other_times_x.mul_scalar(self.coefficients[i]).get_coefficients());
    //             other_times_x = other_times_x.mul_by_x();
    //         }
    //     }
    //     Self { coefficients: coeff_array, degree: degree }
    // }



    // we compute the quotient and remainder unconstrained then check the result with mul and add
    pub fn div<let OTHER_SIZE: u32, let OTHER_LOG_SIZE: u32, let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(self, other: Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> (Polynomial<SIZE, LOG_SIZE>, Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>) {
        assert(SIZE >= OTHER_SIZE);
        let other_degree = other.get_degree();
        if other_degree == 0 {
            assert(other.coefficients[0] != 0, "Division by zero polynomial");
        }

        //We divide the polynomials unconstrained then check the result with mul and add
        let (quotient, remainder) = unsafe {div_polynomials::<SIZE, LOG_SIZE, OTHER_SIZE, OTHER_LOG_SIZE, DOMAIN_SIZE, DOMAIN_LOG_SIZE>(self, other, domain)};
        assert((remainder.get_degree() < other_degree) | (other_degree == 0));

        let mul_result = quotient.mul(other, domain);
        let add_result = mul_result.add(remainder);
        let recovered: Polynomial<SIZE, LOG_SIZE> = add_result.reduce();

        assert(self == recovered);
        (quotient, remainder)
    }



}

unconstrained pub fn div_polynomials<let SIZE: u32, let LOG_SIZE: u32, let OTHER_SIZE: u32, let OTHER_LOG_SIZE: u32, let DOMAIN_SIZE: u32, let DOMAIN_LOG_SIZE: u32>(a: Polynomial<SIZE, LOG_SIZE>, b: Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>, domain: Domain<DOMAIN_SIZE, DOMAIN_LOG_SIZE>) -> (Polynomial<SIZE, LOG_SIZE>, Polynomial<OTHER_SIZE, OTHER_LOG_SIZE>) {
    let mut quotient = Polynomial::<SIZE, LOG_SIZE>::zero();
    let mut remainder = a;
    //Having so many get_degree() calls is not super efficient, but it's ok since we're doing it in an unconstrained context
    while remainder.get_degree() >= b.get_degree() & remainder.is_not_zero(){
        let leading_coeff = remainder.coefficients[remainder.get_degree()] / b.coefficients[b.get_degree()];
        let degree_diff = remainder.get_degree() - b.get_degree();
        
        let mut temp = Polynomial::<SIZE, LOG_SIZE>::zero();
        temp.coefficients[degree_diff] = leading_coeff;
        
        quotient = quotient.add(temp);
        let product: Polynomial<SIZE, LOG_SIZE> = temp.mul(b, domain).reduce();
        remainder = remainder.sub(product);
    }
    
    (quotient, remainder.reduce::<OTHER_SIZE, OTHER_LOG_SIZE>())
}





// ==========================================
// Tests for Polynomial Implementation
// ==========================================


#[test]
fn test_evaluate() {
    // Test polynomial: 2x^2 + 3x + 1
    let coefficients = [1, 3, 2];
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice(coefficients);
    
    // Evaluate at x = 2: 2(4) + 3(2) + 1 = 15
    let result = poly.evaluate(2);
    assert(result == 15);

    // Evaluate at x = 0: should return constant term
    let result = poly.evaluate(0);
    assert(result == 1);

    // Test polynomial: x^3 - 2x + 4
    let coefficients2 = [4, -2, 0, 1];
    let poly2: Polynomial<4, 2> = Polynomial::new_from_slice(coefficients2);
    
    // Evaluate at x = 3: 27 - 6 + 4 = 25
    let result = poly2.evaluate(3);
    assert(result == 25);

    // Test constant polynomial
    let coefficients3 = [5];
    let poly3: Polynomial<1, 0> = Polynomial::new_from_slice(coefficients3);
    
    // Should return 5 for any input
    let result = poly3.evaluate(100);
    assert(result == 5);

    let poly4: Polynomial<2, 1> = Polynomial::new_from_slice([7]);
    assert(poly4.evaluate(10) == 7);
}



#[test]
fn test_new_from_slice() {
    // Test normal case (degree 2), capacity 4
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    assert(poly.get_degree() == 2);
    let coeffs = poly.get_coefficients();
    assert(coeffs[0] == 1);
    assert(coeffs[1] == 2);
    assert(coeffs[2] == 3);
    assert(coeffs[3] == 0);
    
    // Test empty polynomial (degree 0), capacity 1
    let poly_zero: Polynomial<1, 0> = Polynomial::new_from_slice([0]);
    assert(poly_zero.get_degree() == 0);
    let coeffs_zero = poly_zero.get_coefficients();
    assert(coeffs_zero[0] == 0);

    // Test with specific capacity matching slice length (power of 2)
    let poly_deg3_cap4: Polynomial<4, 2> = Polynomial::new_from_slice([1,2,3,4]);
    assert(poly_deg3_cap4.get_degree() == 3);
    assert(poly_deg3_cap4.get_coefficients()[3] == 4);

    // Test with slice smaller than capacity
    let poly_deg0_cap2: Polynomial<2, 1> = Polynomial::new_from_slice([5]);
    assert(poly_deg0_cap2.get_degree() == 0);
    assert(poly_deg0_cap2.get_coefficients()[0] == 5);
    assert(poly_deg0_cap2.get_coefficients()[1] == 0);
}

#[test]
fn test_add() {
    // Basic addition test
    let poly1: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]); // 1 + 2x + 3x^2
    let poly2: Polynomial<4, 2> = Polynomial::new_from_slice([4, 5, 6]); // 4 + 5x + 6x^2
    let result = poly1.add(poly2); // 5 + 7x + 9x^2
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 5);
    assert(coeffs[1] == 7);
    assert(coeffs[2] == 9);
    assert(coeffs[3] == 0);
    assert(result.get_degree() == 2);

    // Test with different capacities
    let poly_large: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3, 0]);
    let poly_small: Polynomial<2, 1> = Polynomial::new_from_slice([1, 1]);
    let result_ls = poly_large.add(poly_small);
    let coeffs_ls = result_ls.get_coefficients();
    assert(coeffs_ls[0] == 2);
    assert(coeffs_ls[1] == 3);
    assert(coeffs_ls[2] == 3);
    assert(coeffs_ls[3] == 0);
    assert(result_ls.get_degree() == 2);

    // Test symmetry
    let poly_small_2: Polynomial<2, 1> = Polynomial::new_from_slice([1,1]);
    let poly_large_2: Polynomial<4, 2> = Polynomial::new_from_slice([1,2,3,0]);
    let result_sl = poly_large_2.add(poly_small_2);
    let coeffs_sl = result_sl.get_coefficients();
    assert(coeffs_sl[0] == 2);
    assert(coeffs_sl[1] == 3);
    assert(coeffs_sl[2] == 3);
    assert(coeffs_sl[3] == 0);

    // Test degree preservation
    let p1_deg1_cap4: Polynomial<4, 2> = Polynomial::new_from_slice([0,1]);
    let p2_deg1_cap4: Polynomial<4, 2> = Polynomial::new_from_slice([1,0]);
    let res_deg1_sum: Polynomial<4, 2> = p1_deg1_cap4.add(p2_deg1_cap4);
    assert(res_deg1_sum.get_coefficients()[0] == 1);
    assert(res_deg1_sum.get_coefficients()[1] == 1);
    assert(res_deg1_sum.get_degree() == 1);
}

#[test]
fn test_sub() {
    // Test normal case: poly1 (deg 2, cap 4), poly2 (deg 2, cap 4)
    let poly1: Polynomial<4, 2> = Polynomial::new_from_slice([10, 20, 30]); // 10 + 20x + 30x^2
    let poly2: Polynomial<4, 2> = Polynomial::new_from_slice([4, 5, 6]);   // 4 + 5x + 6x^2
    let result = poly1.sub(poly2); // Result cap 4: 6 + 15x + 24x^2
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 6);
    assert(coeffs[1] == 15);
    assert(coeffs[2] == 24);
    assert(coeffs[3] == 0);
    assert(result.get_degree() == 2);
    
    // Test case where leading coefficients are equal (should reduce degree)
    // poly3 (deg 2, cap 4), poly4 (deg 2, cap 4)
    let poly3: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]); // 1 + 2x + 3x^2
    let poly4: Polynomial<4, 2> = Polynomial::new_from_slice([1, 1, 3]); // 1 + x + 3x^2
    let result2 = poly3.sub(poly4); // Result cap 4: 0 + x + 0x^2 => x
    
    let coeffs2 = result2.get_coefficients();
    assert(coeffs2[0] == 0);
    assert(coeffs2[1] == 1);
    assert(coeffs2[2] == 0);
    assert(coeffs2[3] == 0);
    assert(result2.get_degree() == 1); // Degree should be 1

    // Test sub with different capacities: poly_large (cap 4), poly_small (cap 2)
    let poly_large: Polynomial<4, 2> = Polynomial::new_from_slice([5, 5, 5, 0]); // 5 + 5x + 5x^2
    let poly_small: Polynomial<2, 1> = Polynomial::new_from_slice([1, 1]);       // 1 + x
    // poly_large.sub(poly_small) -> result cap 4
    let result_ls = poly_large.sub(poly_small); // (5-1) + (5-1)x + 5x^2 = 4 + 4x + 5x^2
    let coeffs_ls = result_ls.get_coefficients();
    assert(coeffs_ls[0] == 4);
    assert(coeffs_ls[1] == 4);
    assert(coeffs_ls[2] == 5);
    assert(coeffs_ls[3] == 0);
    assert(result_ls.get_degree() == 2);
}

#[test]
fn test_neg() {
    // poly (deg 2, cap 4)
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    let result = poly.neg(); // Result cap 4
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == -1);
    assert(coeffs[1] == -2);
    assert(coeffs[2] == -3);
    assert(coeffs[3] == 0);
}

#[test]
fn test_vanishing_polynomial() {
    let x = 5;
    // Vanishing poly is (X - x), degree 1. Needs capacity 2.
    let vanishing: Polynomial<2, 1> = Polynomial::<2, 1>::vanishing_polynomial(x);
    
    let coeffs = vanishing.get_coefficients();
    assert(coeffs[0] == -5);
    assert(coeffs[1] == 1);
    assert(vanishing.get_degree() == 1);
    assert(vanishing.evaluate(x) == 0);

    // Test with another value
    let y = 0;
    let vanishing_at_0: Polynomial<2, 1> = Polynomial::<2, 1>::vanishing_polynomial(y); // X - 0 = X
    let coeffs_y = vanishing_at_0.get_coefficients();
    assert(coeffs_y[0] == 0);
    assert(coeffs_y[1] == 1);
    assert(vanishing_at_0.evaluate(0) == 0);
    assert(vanishing_at_0.evaluate(1) == 1);
}

#[test]
fn test_mul_scalar() {
    // poly (deg 2, cap 4)
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    let scalar = 5;
    let result = poly.mul_scalar(scalar); // Result cap 4
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 5);
    assert(coeffs[1] == 10);
    assert(coeffs[2] == 15);
    assert(coeffs[3] == 0);

    // Test with zero scalar
    let result_zero_scalar = poly.mul_scalar(0);
    let coeffs_zero_scalar = result_zero_scalar.get_coefficients();
    assert(coeffs_zero_scalar[0] == 0);
    assert(coeffs_zero_scalar[1] == 0);
    assert(coeffs_zero_scalar[2] == 0);
    assert(coeffs_zero_scalar[3] == 0);
    assert(result_zero_scalar.get_degree() == 0); // Or handle as special zero poly case
}

#[test]
fn test_mul_by_x() {
    // Test polynomial: 1 + 2x
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2]);
    let result = poly.mul_by_x(); // x + 2x^2
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 0);
    assert(coeffs[1] == 1);
    assert(coeffs[2] == 2);
    assert(coeffs[3] == 0);
    assert(result.get_degree() == 2);

    // Test with constant polynomial: 3
    let poly_const: Polynomial<2, 1> = Polynomial::new_from_slice([3]);
    let result_const = poly_const.mul_by_x(); // 3x
    let coeffs_const = result_const.get_coefficients();
    assert(coeffs_const[0] == 0);
    assert(coeffs_const[1] == 3);
    assert(result_const.get_degree() == 1);

    // Test with polynomial: 3x
    let poly_3x: Polynomial<4, 2> = Polynomial::new_from_slice([0,3]);
    let result_3x_sq = poly_3x.mul_by_x(); // 3x^2
    assert(result_3x_sq.get_coefficients()[0]==0);
    assert(result_3x_sq.get_coefficients()[1]==0);
    assert(result_3x_sq.get_coefficients()[2]==3);
    assert(result_3x_sq.get_degree()==2);
}

#[test]
fn test_mul() {
    let domain: Domain<8, 3> = comptime { Domain::<8, 3>::new() };
    
    // Test (1 + 2x)(2 + 3x) = 2 + 7x + 6x^2
    let poly1: Polynomial<2, 1> = Polynomial::new_from_slice([1, 2]);
    let poly2: Polynomial<2, 1> = Polynomial::new_from_slice([2, 3]);
    let result: Polynomial<4, 2> = poly1.mul(poly2, domain);
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 2);  // constant term
    assert(coeffs[1] == 7);  // x term
    assert(coeffs[2] == 6);  // x^2 term
    assert(coeffs[3] == 0);  // x^3 term
    assert(result.get_degree() == 2);

    // Test (x + 1)(x - 1) = x^2 - 1
    let poly3: Polynomial<2, 1> = Polynomial::new_from_slice([1, 1]);
    let poly4: Polynomial<2, 1> = Polynomial::new_from_slice([-1, 1]);
    let result2: Polynomial<4, 2> = poly3.mul(poly4, domain);
    
    let coeffs2 = result2.get_coefficients();
    assert(coeffs2[0] == -1); // constant term
    assert(coeffs2[1] == 0);  // x term
    assert(coeffs2[2] == 1);  // x^2 term
    assert(coeffs2[3] == 0);
    assert(result2.get_degree() == 2);

    // Test multiplication by constant polynomial
    let poly5: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]); 
    let poly6: Polynomial<1, 0> = Polynomial::new_from_slice([2]);
    let result3: Polynomial<8, 3> = poly5.mul(poly6, domain);
    
    let coeffs3 = result3.get_coefficients();
    assert(coeffs3[0] == 2);   // constant term
    assert(coeffs3[1] == 4);   // x term
    assert(coeffs3[2] == 6);   // x^2 term
    assert(result3.get_degree() == 2);

    // Test constant multiplication
    let p_const1: Polynomial<1, 0> = Polynomial::new_from_slice([3]);
    let p_const2: Polynomial<1, 0> = Polynomial::new_from_slice([5]);
    let res_const_mul: Polynomial<2, 1> = p_const1.mul(p_const2, domain);
    let coeffs_rcm = res_const_mul.get_coefficients();
    assert(coeffs_rcm[0] == 15);
}


#[test]
fn test_div() {
    let domain: Domain<8, 3> = comptime { Domain::<8, 3>::new() };

    // Test case 1: Simple division with no remainder
    // (x^2 + 2x + 1) / (x + 1) = x + 1
    let dividend1: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 1]); // x^2 + 2x + 1
    let divisor1: Polynomial<2, 1> = Polynomial::new_from_slice([1, 1]);     // x + 1
    let (quotient1, remainder1) = dividend1.div(divisor1, domain);
    
    assert(quotient1.get_degree() == 1);
    assert(remainder1.get_degree() == 0);
    let q1_coeffs = quotient1.get_coefficients();
    assert(q1_coeffs[0] == 1); // constant term
    assert(q1_coeffs[1] == 1); // x term
    assert(remainder1.get_coefficients()[0] == 0); // remainder should be 0

    // Test case 2: Division with remainder
    // (x^2 + 3x + 2) / (x + 1) = x + 2 with remainder 0
    let dividend2: Polynomial<4, 2> = Polynomial::new_from_slice([2, 3, 1]); // x^2 + 3x + 2
    let divisor2: Polynomial<2, 1> = Polynomial::new_from_slice([1, 1]);     // x + 1
    let (quotient2, remainder2) = dividend2.div(divisor2, domain);
    
    assert(quotient2.get_degree() == 1);
    assert(remainder2.get_degree() == 0);
    let q2_coeffs = quotient2.get_coefficients();
    assert(q2_coeffs[0] == 2); // constant term
    assert(q2_coeffs[1] == 1); // x term
    assert(remainder2.get_coefficients()[0] == 0); // remainder should be 0

    // Test case 3: Division by constant
    // (2x^2 + 4x + 6) / 2 = x^2 + 2x + 3
    let dividend3: Polynomial<4, 2> = Polynomial::new_from_slice([6, 4, 2]); // 2x^2 + 4x + 6
    let divisor3: Polynomial<1, 0> = Polynomial::new_from_slice([2]);        // 2
    let (quotient3, remainder3) = dividend3.div(divisor3, domain);
    
    assert(quotient3.get_degree() == 2);
    assert(remainder3.get_degree() == 0);
    let q3_coeffs = quotient3.get_coefficients();
    assert(q3_coeffs[0] == 3); // constant term
    assert(q3_coeffs[1] == 2); // x term
    assert(q3_coeffs[2] == 1); // x^2 term
    assert(remainder3.get_coefficients()[0] == 0); // remainder should be 0

    // Test case 4: Division with non-zero remainder
    // (x^3 + 2x^2 + 3x + 4) / (x^2 + 1) = x + 2 with remainder 2x + 2
    let dividend4: Polynomial<4, 2> = Polynomial::new_from_slice([4, 3, 2, 1]); // x^3 + 2x^2 + 3x + 4
    let divisor4: Polynomial<4, 2> = Polynomial::new_from_slice([1, 0, 1]);     // x^2 + 1
    let (quotient4, remainder4) = dividend4.div(divisor4, domain);
    
    assert(quotient4.get_degree() == 1);
    assert(remainder4.get_degree() == 1);
    let q4_coeffs = quotient4.get_coefficients();
    let r4_coeffs = remainder4.get_coefficients();
    assert(q4_coeffs[0] == 2); // constant term
    assert(q4_coeffs[1] == 1); // x term
    assert(r4_coeffs[0] == 2); // remainder constant term
    assert(r4_coeffs[1] == 2); // remainder x term
}


#[test]
fn test_expand_reduce() {
    // Test expand
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    let expanded: Polynomial<8, 3> = poly.expand();
    let expanded_coeffs = expanded.get_coefficients();
    assert(expanded_coeffs[0] == 1);
    assert(expanded_coeffs[1] == 2);
    assert(expanded_coeffs[2] == 3);
    assert(expanded_coeffs[3] == 0);
    assert(expanded_coeffs[4] == 0);

    // Test reduce
    let large_poly: Polynomial<8, 3> = Polynomial::new_from_slice([5, 6, 7, 0, 0]);
    let reduced: Polynomial<4, 2> = large_poly.reduce();
    let reduced_coeffs = reduced.get_coefficients();
    assert(reduced_coeffs[0] == 5);
    assert(reduced_coeffs[1] == 6);
    assert(reduced_coeffs[2] == 7);
    assert(reduced_coeffs[3] == 0);
}

#[test(should_fail)]
fn test_reduce_without_trailing_zeros() {
    // Create a polynomial with non-zero coefficients that we want to "cut off"
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3, 4]);
    
    // This should fail because we're trying to reduce to size 2, 
    // but coefficients at indices 2 and 3 are non-zero (3 and 4)
    let _reduced: Polynomial<2, 1> = poly.reduce();
}

#[test]
fn test_get_degree_edge_cases() {
    // Zero polynomial
    let zero: Polynomial<4, 2> = Polynomial::zero();
    assert(zero.get_degree() == 0);

    // Polynomial with trailing zeros
    let poly_trailing: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 0, 0]);
    assert(poly_trailing.get_degree() == 1);

    // Polynomial with all non-zero coefficients
    let poly_full: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3, 4]);
    assert(poly_full.get_degree() == 3);
}

#[test]
fn test_ntt() {
    let domain = comptime { Domain::<4, 2>::new() };
    let poly: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    
    // Forward NTT
    let lagrange_poly = poly.ntt(domain);
    
    // Verify that we can convert back
    let recovered = lagrange_poly.intt(domain);
    let recovered_coeffs = recovered.get_coefficients();
    let original_coeffs = poly.get_coefficients();
    
    for i in 0..4 {
        assert(recovered_coeffs[i] == original_coeffs[i]);
    }
}

#[test]
fn test_mul_mod() {
    let domain = comptime { Domain::<4, 2>::new() };
    let poly1: Polynomial<4, 2> = Polynomial::new_from_slice([1, 2, 3]);
    let poly2: Polynomial<4, 2> = Polynomial::new_from_slice([4, 5, 6]);
    let result: Polynomial<4, 2> = poly1.mul_mod(poly2, domain);
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 22);
    assert(coeffs[1] == 13);
    assert(coeffs[2] == 28);
    assert(coeffs[3] == 27);
}
