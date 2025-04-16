use crate::MAX_DEGREE::MAX_DEGREE;
use crate::domain::Domain;
pub struct Polynomial {
    coefficients: [Field; MAX_DEGREE + 1],
    degree: u32,
}


impl Eq for Polynomial {
    fn eq(self, other: Self) -> bool {
        (self.coefficients == other.coefficients) & (self.degree == other.degree)
    }
}


impl Polynomial {

    pub fn new(coefficients: [Field; MAX_DEGREE + 1], degree: u32) -> Self {
        assert(degree <= MAX_DEGREE);
        Self { coefficients, degree }
    }

    // Creates a polynomial from a variable-length array of coefficients
    // Coefficients are ordered from constant term to highest degree term
    pub fn new_from_slice(coefficients: [Field]) -> Self {

        let mut coeff_array = [0; MAX_DEGREE + 1];
        let degree = coefficients.len() - 1;
        assert(degree <= MAX_DEGREE, "Degree is too large");

        for i in 0..MAX_DEGREE + 1 {
            if i < degree + 1 {
                coeff_array[i] = coefficients[i];
            }
        }

        Self::new(coeff_array, degree)
    }


    pub fn get_degree(self) -> u32 {
        self.degree
    }

    pub fn get_coefficients(self) -> [Field; MAX_DEGREE + 1] {
        self.coefficients
    }

    pub fn one() -> Self {
        let mut coeff_array = [0; MAX_DEGREE + 1];
        coeff_array[0] = 1;
        Self::new(coeff_array, 0)
    }

    pub fn zero() -> Self {
        let mut coeff_array = [0; MAX_DEGREE + 1];
        Self::new(coeff_array, 0)
    }

    pub fn is_not_zero(self) -> bool {
        self.coefficients.any(|coeff| coeff != 0)
    }

    pub fn eq(self, other: Self) -> bool {
        self.coefficients == other.coefficients
    }
    

    pub fn add(self, other: Self) -> Self {
        let mut coeff_array = add_arrays(self.coefficients, other.coefficients);
        let degree = Self::get_degree_from_coefficients(coeff_array);
        Self { coefficients: coeff_array, degree: degree }
    }

    pub fn sub(self, other: Self) -> Self {
        let mut coeff_array = sub_arrays(self.coefficients, other.coefficients);
        let degree = Self::get_degree_from_coefficients(coeff_array);
        Self { coefficients: coeff_array, degree: degree }
    }

    pub fn neg(self) -> Self {
        let mut coeff_array = [0; MAX_DEGREE + 1];
        for i in 0..MAX_DEGREE + 1 {
            coeff_array[i] = -self.coefficients[i];
        }
        Self { coefficients: coeff_array, degree: self.degree }
    }

    // Creates a polynomial (x - a) that evaluates to 0 at point x = a
    pub fn vanishing_polynomial(x: Field) -> Self {
        let mut coeff_array = [0; MAX_DEGREE + 1];
        coeff_array[0] = -x;
        coeff_array[1] = 1;
        Self::new(coeff_array, 1)
    }

    pub fn mul_scalar(self, scalar: Field) -> Self {
        let mut coeff_array = [0; MAX_DEGREE + 1];
        for i in 0..MAX_DEGREE + 1 {
            coeff_array[i] = self.coefficients[i] * scalar;
        }
        Self { coefficients: coeff_array, degree: self.degree }
    }

    // Determines the actual degree by finding the highest non-zero coefficient
    pub fn get_degree_from_coefficients(coefficients: [Field; MAX_DEGREE + 1]) -> u32 {
        let mut degree = 0;
        for i in 0..MAX_DEGREE + 1 {
            if coefficients[i] != 0 {
                degree = i;
            }
        }
        degree
    }

    pub fn evaluate(self, x: Field) -> Field {
        let mut result: Field = 0;
        let mut x_power = 1;
        for i in 0..MAX_DEGREE + 1 {
            result += self.coefficients[i] * x_power;
            x_power *= x;
        }
        result
    }

    pub fn mul_by_x(self) -> Self {
        assert(self.degree < MAX_DEGREE);
        let degree = self.degree + 1;
        let mut coeff_array = [0; MAX_DEGREE + 1];
        for i in 1..MAX_DEGREE + 1 {
            coeff_array[i] = self.coefficients[i - 1];
        }
        Self { coefficients: coeff_array, degree: degree }
    }

    // Basic polynomial multiplication using coefficient-wise operations
    // TODO: Improve this function:
    // - Use Karatsuba algorithm for O(n^1.58) complexity
    // - Use FFT for O(n log n) complexity when applicable
    // - Consider using unconstrained operations with verification, if hash is fast, we can check if f(x) * g(x) = h(x) with x pseudo-randomly chosen
    pub fn mul(self, other: Self) -> Self {
        assert(self.degree + other.degree <= MAX_DEGREE);
        let degree = self.degree + other.degree;
        let mut coeff_array = [0; MAX_DEGREE + 1];
        let mut other_times_x = other;
        for i in 0..MAX_DEGREE + 1 {
            if i <= self.degree {
                coeff_array = add_arrays(coeff_array, other_times_x.mul_scalar(self.coefficients[i]).get_coefficients());
                other_times_x = other_times_x.mul_by_x();
            }
        }
        Self { coefficients: coeff_array, degree: degree }
    }

    // Polynomial division returning quotient and remainder
    // Requires that self.degree >= other.degree
    // we compute the quotient and remainder unconstrained then check the result with mul and add
    // TODO: Improve this function:
    pub fn div(self, other: Self) -> (Self, Self) {
        assert(self.degree >= other.degree);
        
        let result = unsafe {div_polynomials(self, other)};
        assert(result.0.degree == self.degree - other.degree);
        assert((result.1.degree < other.degree) | (result.1.degree == 0));
        assert(self == other.mul(result.0).add(result.1));
        result
    }

}

pub fn multiplier<let N: u32,let M: u32>(a: [Field; N], b: [Field; M]) -> [Field; N + M - 1] {
    let mut result = [0; N + M - 1];
    for i in 0..N {
        for j in 0..M {
            result[i + j] += a[i] * b[j];
        }
    }
    result
}

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

unconstrained pub fn div_polynomials(a: Polynomial, b: Polynomial) -> (Polynomial, Polynomial) {
    let mut quotient = Polynomial::zero();
    let mut remainder = a;
    
    while remainder.get_degree() >= b.get_degree() & remainder.is_not_zero(){
        let leading_coeff = remainder.coefficients[remainder.get_degree()] / b.coefficients[b.get_degree()];
        let degree_diff = remainder.get_degree() - b.get_degree();
        
        let mut temp = Polynomial::zero();
        temp.coefficients[degree_diff] = leading_coeff;
        temp.degree = degree_diff;
        
        quotient = quotient.add(temp);
        remainder = remainder.sub(b.mul(temp));
    }
    
    (quotient, remainder)
}


// ==========================================
// Tests for Polynomial Implementation
// ==========================================


#[test]
fn test_evaluate() {
    // Test polynomial: 2x^2 + 3x + 1
    let coefficients = [1, 3, 2];
    let poly = Polynomial::new_from_slice(coefficients);
    
    // Evaluate at x = 2: 2(2^2) + 3(2) + 1 = 8 + 6 + 1 = 15
    let result = poly.evaluate(2);
    assert(result == 15);

    // Evaluate at x = 0: should return constant term
    let result = poly.evaluate(0);
    assert(result == 1);

    // Test polynomial: x^3 - 2x + 4
    let coefficients2 = [4, -2, 0, 1];
    let poly2 = Polynomial::new_from_slice(coefficients2);
    
    // Evaluate at x = 3: 3^3 - 2(3) + 4 = 27 - 6 + 4 = 25
    let result = poly2.evaluate(3);
    assert(result == 25);

    // Test constant polynomial: 5
    let coefficients3 = [5];
    let poly3 = Polynomial::new_from_slice(coefficients3);
    
    // Should return 5 for any input
    let result = poly3.evaluate(100);
    assert(result == 5);
}



#[test]
fn test_new_from_slice() {
    // Test normal case
    let poly = Polynomial::new_from_slice([1, 2, 3]);
    assert(poly.get_degree() == 2);
    let coeffs = poly.get_coefficients();
    assert(coeffs[0] == 1);
    assert(coeffs[1] == 2);
    assert(coeffs[2] == 3);
    
    // Test empty polynomial
    let poly_zero = Polynomial::new_from_slice([0]);
    assert(poly_zero.get_degree() == 0);
}

#[test]
fn test_add() {
    let poly1 = Polynomial::new_from_slice([1, 2, 3]);
    let poly2 = Polynomial::new_from_slice([4, 5, 6]);
    let result = poly1.add(poly2);
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 5);
    assert(coeffs[1] == 7);
    assert(coeffs[2] == 9);
}

#[test]
fn test_sub() {
    // Test normal case
    let poly1 = Polynomial::new_from_slice([10, 20, 30]);
    let poly2 = Polynomial::new_from_slice([4, 5, 6]);
    let result = poly1.sub(poly2);
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 6);
    assert(coeffs[1] == 15);
    assert(coeffs[2] == 24);
    assert(result.get_degree() == 2);
    
    // Test case where leading coefficients are equal (should reduce degree)
    let poly3 = Polynomial::new_from_slice([1, 2, 3]);
    let poly4 = Polynomial::new_from_slice([1, 1, 3]);
    let result2 = poly3.sub(poly4);
    
    let coeffs2 = result2.get_coefficients();
    assert(coeffs2[0] == 0);
    assert(coeffs2[1] == 1);
    assert(result2.get_degree() == 1); // Degree should be 1, not 2
}

#[test]
fn test_neg() {
    let poly = Polynomial::new_from_slice([1, 2, 3]);
    let result = poly.neg();
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == -1);
    assert(coeffs[1] == -2);
    assert(coeffs[2] == -3);
}

#[test]
fn test_vanishing_polynomial() {
    let x = 5;
    let vanishing = Polynomial::vanishing_polynomial(x);
    
    let coeffs = vanishing.get_coefficients();
    assert(coeffs[0] == -5);
    assert(coeffs[1] == 1);
    assert(vanishing.get_degree() == 1);
    assert(vanishing.evaluate(x) == 0);
}

#[test]
fn test_mul_scalar() {
    let poly = Polynomial::new_from_slice([1, 2, 3]);
    let scalar = 5;
    let result = poly.mul_scalar(scalar);
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 5);
    assert(coeffs[1] == 10);
    assert(coeffs[2] == 15);
}

#[test]
fn test_mul_by_x() {
    // Test polynomial: 2x + 1 -> x(2x + 1) = 2x^2 + x
    let poly = Polynomial::new_from_slice([1, 2]);
    let result = poly.mul_by_x();
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 0);
    assert(coeffs[1] == 1);
    assert(coeffs[2] == 2);
    assert(result.get_degree() == 2);
}

#[test]
fn test_mul() {
    // Test case 1: (2x + 1)(3x + 2) = 6x^2 + 7x + 2
    let poly1 = Polynomial::new_from_slice([1, 2]);
    let poly2 = Polynomial::new_from_slice([2, 3]);
    let result = poly1.mul(poly2);
    
    let coeffs = result.get_coefficients();
    assert(coeffs[0] == 2);  // constant term
    assert(coeffs[1] == 7);  // x term
    assert(coeffs[2] == 6);  // x^2 term
    assert(result.get_degree() == 2);

    // Test case 2: (x + 1)(x - 1) = x^2 - 1
    let poly3 = Polynomial::new_from_slice([1, 1]);
    let poly4 = Polynomial::new_from_slice([-1, 1]);
    let result2 = poly3.mul(poly4);
    
    let coeffs2 = result2.get_coefficients();
    assert(coeffs2[0] == -1);  // constant term
    assert(coeffs2[1] == 0);   // x term
    assert(coeffs2[2] == 1);   // x^2 term
    assert(result2.get_degree() == 2);

    // Test case 3: multiplication by constant polynomial
    let poly5 = Polynomial::new_from_slice([1, 2, 3]);  // 3x^2 + 2x + 1
    let poly6 = Polynomial::new_from_slice([2]);        // 2
    let result3 = poly5.mul(poly6);
    
    let coeffs3 = result3.get_coefficients();
    assert(coeffs3[0] == 2);   // constant term
    assert(coeffs3[1] == 4);   // x term
    assert(coeffs3[2] == 6);   // x^2 term
    assert(result3.get_degree() == 2);


     let poly7 = Polynomial::new_from_slice([42, 87, 15, 63, 29, 91, 8]);  // 3x^2 + 2x + 1
    let poly8 = Polynomial::new_from_slice([77, 34, 5, 68, 92, 13, 56, 70, 24]);        // 2
    let result4 = poly7.mul(poly8);
    let poly9 = Polynomial::new_from_slice([3234, 8127, 4323, 8652, 14230, 17878, 13002, 16502, 17653, 15959, 8313, 8742, 7514, 2744, 192]);
    
    let coeffs4 = result4.get_coefficients();
    for i in 0..MAX_DEGREE + 1 {
        assert(coeffs4[i] == poly9.coefficients[i]);
    }
    assert(result4.get_degree() == 14);


}

#[test]
fn test_div() {
    // Test case 1: Exact division (x^2 + 2x + 1) / (x + 1) = x + 1
    let dividend = Polynomial::new_from_slice([1, 2, 1]);  // x^2 + 2x + 1
    let divisor = Polynomial::new_from_slice([1, 1]);      // x + 1
    let (quotient, remainder) = dividend.div(divisor);
    
    let expected_quotient = Polynomial::new_from_slice([1, 1]);  // x + 1
    let expected_remainder = Polynomial::zero();
    
    assert(quotient == expected_quotient);
    assert(remainder == expected_remainder);

    //Test case 2: Division with remainder (x^3 + 2x^2 + 3x + 4) / (x^2 + 1) = x + 2, remainder 2x + 2
    let dividend2 = Polynomial::new_from_slice([4, 3, 2, 1]);  // x^3 + 2x^2 + 3x + 4
    let divisor2 = Polynomial::new_from_slice([1, 0, 1]);      // x^2 + 1
    let (quotient2, remainder2) = dividend2.div(divisor2);
    
    let expected_quotient2 = Polynomial::new_from_slice([2, 1]);  // x + 2
    let expected_remainder2 = Polynomial::new_from_slice([2, 2]); // 2x + 2
    
    assert(quotient2 == expected_quotient2);
    assert(remainder2 == expected_remainder2);

    // // Test case 3: Division by constant (2x^2 + 4x + 6) / 2 = x^2 + 2x + 3
    let dividend3 = Polynomial::new_from_slice([6, 4, 2]);  // 2x^2 + 4x + 6
    let divisor3 = Polynomial::new_from_slice([2]);         // 2
    let (quotient3, remainder3) = dividend3.div(divisor3);
    
    let expected_quotient3 = Polynomial::new_from_slice([3, 2, 1]);  // x^2 + 2x + 3
    let expected_remainder3 = Polynomial::zero();
    
    assert(quotient3 == expected_quotient3);
    assert(remainder3 == expected_remainder3);
}

