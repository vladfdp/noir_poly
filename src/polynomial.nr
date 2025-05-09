use crate::MAX_DEGREE::MAX_DEGREE;
use crate::MAX_DEGREE::MAX_DEGREE_LOG;
use crate::MAX_DEGREE::NORMALIZER;
use crate::domain::Domain;
use crate::array_utils::{add_arrays, sub_arrays, neg_array, scalar_multiply_array, pointwise_multiply_arrays, CT_butterfly, bit_reverse_array};
use crate::lagrange_polynomial::LagrangePolynomial;

pub struct Polynomial<let N: u32> {
    coefficients: [Field; N],
}


impl<let N: u32> Eq for Polynomial<N> {
    fn eq(self, other: Self) -> bool {
        (self.coefficients == other.coefficients)
    }
}


impl<let N: u32> Polynomial<N> {

    pub fn new(coefficients: [Field; N]) -> Self {
        Self { coefficients }
    }

    // Creates a polynomial from a variable-length slice of coefficients
    // Coefficients are ordered from constant term to highest degree term
    pub fn new_from_slice(coefficients: [Field]) -> Self {

        let mut coeff_array = [0; N];

        for i in 0..N {
            if i < coefficients.len() {
                coeff_array[i] = coefficients[i];
            }
        }

        Self::new(coeff_array)
    }


    pub fn get_degree(self) -> u32 {
        get_last_non_zero_index(self.coefficients)
    }

    pub fn get_coefficients(self) -> [Field; N] {
        self.coefficients
    }

    pub fn one() -> Self {
        let mut coeff_array = [0; N];
        coeff_array[0] = 1;
        Self::new(coeff_array, 0)
    }

    pub fn zero() -> Self {
        let mut coeff_array = [0; N];
        Self::new(coeff_array, 0)
    }

    pub fn is_not_zero(self) -> bool {
        self.coefficients.any(|coeff| coeff != 0)
    }

    // Expands a polynomial to a larger capacity by adding zeros
    pub fn expand<let M: u32>(self) -> Polynomial<M> {
        assert(N <= M);
        let mut coeff_array = [0; M];
        for i in 0..N {
            coeff_array[i] = self.coefficients[i];
        }
        Self::new(coeff_array)
    }

    // Reduces a polynomial to a smaller capacity by removing trailing zeros
    pub fn reduce<let M: u32>(self) -> Polynomial<M> {
        assert(N >= M);
        let mut coeff_array = [0; M];
        for i in 0..M {
            coeff_array[i] = self.coefficients[i];
        }
        for i in M..N {
            assert(self.coefficients[i] == 0);
        }
        Self::new(coeff_array)
    }
    

    pub fn add<let M: u32>(self, other: Polynomial<M>) -> Self {
        let mut coeff_array = add_arrays::<N, M>(self.coefficients, other.coefficients);
        Self { coefficients: coeff_array }
    }

    pub fn sub<let M: u32>(self, other: Polynomial<M>) -> Self {
        let mut coeff_array = sub_arrays::<N, M>(self.coefficients, other.coefficients);
        Self { coefficients: coeff_array }
    }

    pub fn neg(self) -> Self {
        Self { coefficients: neg_array(self.coefficients) }
    }

    // Creates a polynomial (x - a) that evaluates to 0 at point x = a
    pub fn vanishing_polynomial(a: Field) -> Self {
        let mut coeff_array = [0; N];
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
        for i in 0..N {
            result += self.coefficients[i] * x_power;
            x_power *= x;
        }
        result
    }

    pub fn mul_by_x(self) -> Self {
        assert(self.coefficients[N - 1] == 0);
        let mut coeff_array = [0; N];
        for i in 1..N {
            coeff_array[i] = self.coefficients[i - 1];
        }
        Self { coefficients: coeff_array }
    }

    pub fn ntt<let L: u32>(self, domain: Domain<L>) -> LagrangePolynomial<N> {
        let mut new_coeffs = bit_reverse_array::<N>(self.coefficients);
        new_coeffs = CT_butterfly::<N, L>(new_coeffs, domain);
        LagrangePolynomial::<N>::new(new_coeffs)
    }

    // Polynomial multiplication using ntt and intt
    // - Consider using unconstrained operations with verification, if hash is fast, we can check if f(x) * g(x) = h(x) with x pseudo-randomly chosen
    pub fn mul<let M: u32, let L: u32>(self, other: Polynomial<M>, domain: Domain<L>) -> Polynomial<N + N> {
        assert(N >= M);
        let expanded_self = self.expand::<N + N>();
        let expanded_other = other.expand::<N + N>();
        let lagrange_self = expanded_self.ntt::<N + N, L>(domain);
        let lagrange_other = expanded_other.ntt::<N + N, L>(domain);
        let lagrange_result = lagrange_self.mul(lagrange_other);
        lagrange_result.intt::<N + N, L>(domain)
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
    pub fn div<let M: u32, let L: u32>(self, other: Polynomial<M>, domain: Domain<L>) -> (Polynomial<N>, Polynomial<M>) {
        assert(N >= M);
        let other_degree = other.get_degree();
        if other_degree == 0 {
            assert(other.coefficients[0] != 0, "Division by zero polynomial");
        }

        //We divide the polynomials unconstrained then check the result with mul and add
        let (quotient, remainder) = unsafe {div_polynomials::<N, M, L>(self, other, domain)};
        assert(remainder.get_degree() < other_degree | other_degree == 0);

        let mul_result = quotient.mul::<M, L>(other, domain);
        let add_result = mul_result.add::<N+N, M>(remainder);
        let recovered: Polynomial<N> = add_result.reduce::<N>();

        assert(self == recovered);
        (quotient, remainder)
    }



}

unconstrained pub fn div_polynomials<let N: u32, let M: u32, let L: u32>(a: Polynomial<N>, b: Polynomial<M>, domain: Domain<L>) -> (Polynomial<N>, Polynomial<M>) {
    let mut quotient = Polynomial::<N>::zero();
    let mut remainder = a;
    //Having so many get_degree() calls is not super efficient, but it's ok since we're doing it in an unconstrained context
    while remainder.get_degree() >= b.get_degree() & remainder.is_not_zero(){
        let leading_coeff = remainder.coefficients[remainder.get_degree()] / b.coefficients[b.get_degree()];
        let degree_diff = remainder.get_degree() - b.get_degree();
        
        let mut temp = Polynomial::<N>::zero();
        temp.coefficients[degree_diff] = leading_coeff;
        
        quotient = quotient.add::<N>(temp);
        let product = temp.mul::<M, L>(b, domain);
        remainder = remainder.sub::<N+N>(product);
    }
    
    (quotient, remainder.reduce::<M>())
}





// ==========================================
// Tests for Polynomial Implementation
// ==========================================


// #[test]
// fn test_evaluate() {
//     // Test polynomial: 2x^2 + 3x + 1 (degree 2), capacity 4
//     let coefficients = [1, 3, 2];
//     let poly: Polynomial<4> = Polynomial::new_from_slice(coefficients);
    
//     // Evaluate at x = 2: 2(2^2) + 3(2) + 1 = 8 + 6 + 1 = 15
//     let result = poly.evaluate(2);
//     assert(result == 15);

//     // Evaluate at x = 0: should return constant term
//     let result = poly.evaluate(0);
//     assert(result == 1);

//     // Test polynomial: x^3 - 2x + 4 (degree 3), capacity 4
//     let coefficients2 = [4, -2, 0, 1];
//     let poly2: Polynomial<4> = Polynomial::new_from_slice(coefficients2);
    
//     // Evaluate at x = 3: 3^3 - 2(3) + 4 = 27 - 6 + 4 = 25
//     let result = poly2.evaluate(3);
//     assert(result == 25);

//     // Test constant polynomial: 5 (degree 0), capacity 1
//     let coefficients3 = [5];
//     let poly3: Polynomial<1> = Polynomial::new_from_slice(coefficients3);
    
//     // Should return 5 for any input
//     let result = poly3.evaluate(100);
//     assert(result == 5);

//     // Test with capacity 2 for a degree 0 polynomial
//     let poly4: Polynomial<2> = Polynomial::new_from_slice([7]);
//     assert(poly4.evaluate(10) == 7);
// }



// #[test]
// fn test_new_from_slice() {
//     // Test normal case (degree 2), capacity 4
//     let poly: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]);
//     assert(poly.get_degree() == 2);
//     let coeffs = poly.get_coefficients();
//     assert(coeffs[0] == 1);
//     assert(coeffs[1] == 2);
//     assert(coeffs[2] == 3);
//     assert(coeffs[3] == 0);
    
//     // Test empty polynomial (degree 0), capacity 1
//     let poly_zero: Polynomial<1> = Polynomial::new_from_slice([0]);
//     assert(poly_zero.get_degree() == 0);
//     let coeffs_zero = poly_zero.get_coefficients();
//     assert(coeffs_zero[0] == 0);

//     // Test with specific capacity matching slice length (power of 2)
//     let poly_deg3_cap4: Polynomial<4> = Polynomial::new_from_slice([1,2,3,4]);
//     assert(poly_deg3_cap4.get_degree() == 3);
//     assert(poly_deg3_cap4.get_coefficients()[3] == 4);

//     // Test with slice smaller than capacity
//     let poly_deg0_cap2: Polynomial<2> = Polynomial::new_from_slice([5]);
//     assert(poly_deg0_cap2.get_degree() == 0);
//     assert(poly_deg0_cap2.get_coefficients()[0] == 5);
//     assert(poly_deg0_cap2.get_coefficients()[1] == 0);
// }

// #[test]
// fn test_add() {
//     // Both polynomials degree 2, capacity 4
//     let poly1: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]); // 1 + 2x + 3x^2
//     let poly2: Polynomial<4> = Polynomial::new_from_slice([4, 5, 6]); // 4 + 5x + 6x^2
//     let result = poly1.add(poly2); // Result capacity 4: 5 + 7x + 9x^2
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == 5);
//     assert(coeffs[1] == 7);
//     assert(coeffs[2] == 9);
//     assert(coeffs[3] == 0);
//     assert(result.get_degree() == 2);

//     // Test add with different capacities: poly_large (cap 4), poly_small (cap 2)
//     let poly_large: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3, 0]); // 1 + 2x + 3x^2
//     let poly_small: Polynomial<2> = Polynomial::new_from_slice([1, 1]);       // 1 + x
//     // poly_large.add(poly_small) -> result cap 4
//     let result_ls = poly_large.add(poly_small);
//     let coeffs_ls = result_ls.get_coefficients();
//     assert(coeffs_ls[0] == 2); // 1+1
//     assert(coeffs_ls[1] == 3); // 2+1
//     assert(coeffs_ls[2] == 3); // 3+0
//     assert(coeffs_ls[3] == 0);
//     assert(result_ls.get_degree() == 2);

//     // Test add where smaller capacity is first (should ensure self.capacity is large enough for result)
//     // If poly_small.add(poly_large) -> result cap 2, this would truncate.
//     // This depends on add_arrays behavior; assuming self.capacity must be >= other.capacity
//     // To be safe and follow the "bigger first" rule:
//     let poly_small_2: Polynomial<2> = Polynomial::new_from_slice([1,1]);
//     let poly_large_2: Polynomial<4> = Polynomial::new_from_slice([1,2,3,0]);
//     // Call poly_large_2.add(poly_small_2)
//     let result_sl = poly_large_2.add(poly_small_2); // Should be same as result_ls
//     let coeffs_sl = result_sl.get_coefficients();
//     assert(coeffs_sl[0] == 2);
//     assert(coeffs_sl[1] == 3);
//     assert(coeffs_sl[2] == 3);
//     assert(coeffs_sl[3] == 0);

//     // Test add resulting in degree increase (still within capacity)
//     let p1_deg1_cap4: Polynomial<4> = Polynomial::new_from_slice([0,1]); // x
//     let p2_deg1_cap4: Polynomial<4> = Polynomial::new_from_slice([1,0]); // 1
//     let res_deg1_sum: Polynomial<4> = p1_deg1_cap4.add(p2_deg1_cap4); // 1+x
//     assert(res_deg1_sum.get_coefficients()[0] == 1);
//     assert(res_deg1_sum.get_coefficients()[1] == 1);
//     assert(res_deg1_sum.get_degree() == 1);
// }

// #[test]
// fn test_sub() {
//     // Test normal case: poly1 (deg 2, cap 4), poly2 (deg 2, cap 4)
//     let poly1: Polynomial<4> = Polynomial::new_from_slice([10, 20, 30]); // 10 + 20x + 30x^2
//     let poly2: Polynomial<4> = Polynomial::new_from_slice([4, 5, 6]);   // 4 + 5x + 6x^2
//     let result = poly1.sub(poly2); // Result cap 4: 6 + 15x + 24x^2
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == 6);
//     assert(coeffs[1] == 15);
//     assert(coeffs[2] == 24);
//     assert(coeffs[3] == 0);
//     assert(result.get_degree() == 2);
    
//     // Test case where leading coefficients are equal (should reduce degree)
//     // poly3 (deg 2, cap 4), poly4 (deg 2, cap 4)
//     let poly3: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]); // 1 + 2x + 3x^2
//     let poly4: Polynomial<4> = Polynomial::new_from_slice([1, 1, 3]); // 1 + x + 3x^2
//     let result2 = poly3.sub(poly4); // Result cap 4: 0 + x + 0x^2 => x
    
//     let coeffs2 = result2.get_coefficients();
//     assert(coeffs2[0] == 0);
//     assert(coeffs2[1] == 1);
//     assert(coeffs2[2] == 0);
//     assert(coeffs2[3] == 0);
//     assert(result2.get_degree() == 1); // Degree should be 1

//     // Test sub with different capacities: poly_large (cap 4), poly_small (cap 2)
//     let poly_large: Polynomial<4> = Polynomial::new_from_slice([5, 5, 5, 0]); // 5 + 5x + 5x^2
//     let poly_small: Polynomial<2> = Polynomial::new_from_slice([1, 1]);       // 1 + x
//     // poly_large.sub(poly_small) -> result cap 4
//     let result_ls = poly_large.sub(poly_small); // (5-1) + (5-1)x + 5x^2 = 4 + 4x + 5x^2
//     let coeffs_ls = result_ls.get_coefficients();
//     assert(coeffs_ls[0] == 4);
//     assert(coeffs_ls[1] == 4);
//     assert(coeffs_ls[2] == 5);
//     assert(coeffs_ls[3] == 0);
//     assert(result_ls.get_degree() == 2);
// }

// #[test]
// fn test_neg() {
//     // poly (deg 2, cap 4)
//     let poly: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]);
//     let result = poly.neg(); // Result cap 4
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == -1);
//     assert(coeffs[1] == -2);
//     assert(coeffs[2] == -3);
//     assert(coeffs[3] == 0);
// }

// #[test]
// fn test_vanishing_polynomial() {
//     let x = 5;
//     // Vanishing poly is (X - x), degree 1. Needs capacity 2.
//     let vanishing: Polynomial<2> = Polynomial::<2>::vanishing_polynomial(x);
    
//     let coeffs = vanishing.get_coefficients();
//     assert(coeffs[0] == -5);
//     assert(coeffs[1] == 1);
//     assert(vanishing.get_degree() == 1);
//     assert(vanishing.evaluate(x) == 0);

//     // Test with another value
//     let y = 0;
//     let vanishing_at_0: Polynomial<2> = Polynomial::<2>::vanishing_polynomial(y); // X - 0 = X
//     let coeffs_y = vanishing_at_0.get_coefficients();
//     assert(coeffs_y[0] == 0);
//     assert(coeffs_y[1] == 1);
//     assert(vanishing_at_0.evaluate(0) == 0);
//     assert(vanishing_at_0.evaluate(1) == 1);
// }

// #[test]
// fn test_mul_scalar() {
//     // poly (deg 2, cap 4)
//     let poly: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]);
//     let scalar = 5;
//     let result = poly.mul_scalar(scalar); // Result cap 4
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == 5);
//     assert(coeffs[1] == 10);
//     assert(coeffs[2] == 15);
//     assert(coeffs[3] == 0);

//     // Test with zero scalar
//     let result_zero_scalar = poly.mul_scalar(0);
//     let coeffs_zero_scalar = result_zero_scalar.get_coefficients();
//     assert(coeffs_zero_scalar[0] == 0);
//     assert(coeffs_zero_scalar[1] == 0);
//     assert(coeffs_zero_scalar[2] == 0);
//     assert(coeffs_zero_scalar[3] == 0);
//     assert(result_zero_scalar.get_degree() == 0); // Or handle as special zero poly case
// }

// #[test]
// fn test_mul_by_x() {
//     // Test polynomial: 2x + 1 (deg 1). Choose capacity 4 to allow for degree increase.
//     // P(x) = 1 + 2x. Coeffs [1, 2].
//     // x * P(x) = x + 2x^2. Coeffs [0, 1, 2]. Degree 2.
//     let poly: Polynomial<4> = Polynomial::new_from_slice([1, 2]);
//     let result = poly.mul_by_x(); // Result cap 4
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == 0);
//     assert(coeffs[1] == 1);
//     assert(coeffs[2] == 2);
//     assert(coeffs[3] == 0);
//     assert(result.get_degree() == 2);

//     // Test with polynomial of degree 0: P(x) = 3. Cap 2.
//     // x * P(x) = 3x. Degree 1.
//     let poly_const: Polynomial<2> = Polynomial::new_from_slice([3]);
//     let result_const = poly_const.mul_by_x(); // Result cap 2
//     let coeffs_const = result_const.get_coefficients();
//     assert(coeffs_const[0] == 0);
//     assert(coeffs_const[1] == 3);
//     assert(result_const.get_degree() == 1);

//     // Test when highest coefficient is already at N-2, so x*P(x) fills N-1
//     // P(x) = x. Coeffs [0,1]. Cap 2.
//     // x * P(x) = x^2. Coeffs [0,0,1]. This requires capacity 4 for mul_by_x if result must have space for x^2.
//     // The current mul_by_x returns Self, so N must be large enough for x * P(x).
//     // If P(x) = x (coeffs [0,1]) with poly_deg1_cap2: Polynomial<2>
//     // mul_by_x shifts coefficients: [0,0]. This is wrong.
//     // The assert `self.coefficients[N-1] == 0` is key.
//     // For P(x) = 1+2x, let N=4. Coeffs [1,2,0,0]. N-1=3. coeff[3]=0. OK. Result [0,1,2,0].
//     // For P(x) = 3x (deg 1), let's use N=4. Coeffs [0,3,0,0]. N-1=3. coeff[3]=0. OK. Result [0,0,3,0] (3x^2).
//     let poly_3x: Polynomial<4> = Polynomial::new_from_slice([0,3]);
//     let result_3x_sq = poly_3x.mul_by_x();
//     assert(result_3x_sq.get_coefficients()[0]==0);
//     assert(result_3x_sq.get_coefficients()[1]==0);
//     assert(result_3x_sq.get_coefficients()[2]==3);
//     assert(result_3x_sq.get_degree()==2);
// }

// #[test]
// fn test_mul() {
//     // Test case 1: (2x + 1)(3x + 2) = 6x^2 + 7x + 2
//     // P1(x) = 1 + 2x (deg 1), N1=2. P2(x) = 2 + 3x (deg 1), N2=2.
//     // Result deg 2. mul result capacity N1+N1 = 2+2=4. Domain capacity L = N1+N1 = 4.
//     let poly1: Polynomial<2> = Polynomial::new_from_slice([1, 2]);
//     let poly2: Polynomial<2> = Polynomial::new_from_slice([2, 3]);
//     let domain: Domain<4> = comptime { Domain::<4>::new() };
//     let result: Polynomial<4> = poly1.mul::<2, 4>(poly2, domain);
    
//     let coeffs = result.get_coefficients();
//     assert(coeffs[0] == 2);  // constant term
//     assert(coeffs[1] == 7);  // x term
//     assert(coeffs[2] == 6);  // x^2 term
//     assert(coeffs[3] == 0);  // x^3 term (zero)
//     assert(result.get_degree() == 2);

//     // Test case 2: (x + 1)(x - 1) = x^2 - 1
//     // P3(x) = 1 + x (deg 1), N3=2. P4(x) = -1 + x (deg 1), N4=2.
//     // Result deg 2. mul result capacity N3+N3 = 2+2=4. Domain capacity L = N3+N3 = 4.
//     let poly3: Polynomial<2> = Polynomial::new_from_slice([1, 1]);
//     let poly4: Polynomial<2> = Polynomial::new_from_slice([-1, 1]);
//     // Domain can be reused if size matches
//     let result2: Polynomial<4> = poly3.mul::<2, 4>(poly4, domain); // Using domain from above
    
//     let coeffs2 = result2.get_coefficients();
//     assert(coeffs2[0] == -1); // constant term
//     assert(coeffs2[1] == 0);  // x term
//     assert(coeffs2[2] == 1);  // x^2 term
//     assert(coeffs2[3] == 0);
//     assert(result2.get_degree() == 2);

//     // Test case 3: multiplication by constant polynomial
//     // P5(x) = 1 + 2x + 3x^2 (deg 2), N5=4. P6(x) = 2 (deg 0), N6=1.
//     // N5 >= N6 (4 >= 1).
//     // Result deg 2. mul result capacity N5+N5 = 4+4=8. Domain capacity L = N5+N5 = 8.
//     let poly5: Polynomial<4> = Polynomial::new_from_slice([1, 2, 3]); 
//     let poly6: Polynomial<1> = Polynomial::new_from_slice([2]);       
//     let domain_8: Domain<8> = comptime { Domain::<8>::new() };
//     let result3: Polynomial<8> = poly5.mul::<1, 8>(poly6, domain_8);
    
//     let coeffs3 = result3.get_coefficients();
//     assert(coeffs3[0] == 2);   // constant term: 1*2=2
//     assert(coeffs3[1] == 4);   // x term: 2*2=4
//     assert(coeffs3[2] == 6);   // x^2 term: 3*2=6
//     for i in 3..8 { assert(coeffs3[i] == 0); }
//     assert(result3.get_degree() == 2);

//     // Test case 4: P1 (deg 0, N1=1), P2 (deg 0, N2=1)
//     // P1(x)=3, P2(x)=5. Result P(x)=15 (deg 0).
//     // N1=1, N2=1. N1+N1=2. Domain L=2.
//     let p_const1: Polynomial<1> = Polynomial::new_from_slice([3]);
//     let p_const2: Polynomial<1> = Polynomial::new_from_slice([5]);
//     let domain_2: Domain<2> = comptime { Domain::<2>::new() };
//     let res_const_mul: Polynomial<2> = p_const1.mul::<1, 2>(p_const2, domain_2);
//     let coeffs_rcm = res_const_mul.get_coefficients();
//     assert(coeffs_rcm[0] == 15);
//     assert(coeffs_rcm[1] == 0);
//     assert(res_const_mul.get_degree() == 0);


//     //  let poly7 = Polynomial::new_from_slice([42, 87, 15, 63, 29, 91, 8]);  // 3x^2 + 2x + 1
//     // assert(result4.get_degree() == 14);


// }

// #[test]
// fn test_div() {
//     // Test case 1: Exact division (x^2 + 2x + 1) / (x + 1) = x + 1
//     // Dividend: 1+2x+x^2 (deg 2). Capacity N1=4.
//     // Divisor: 1+x (deg 1). Capacity N2=2.
//     // N1 >= N2 (4 >= 2).
//     // Domain for mul check: L = N1+N1 = 4+4=8.
//     let dividend: Polynomial<4> = Polynomial::new_from_slice([1, 2, 1]); 
//     let divisor: Polynomial<2> = Polynomial::new_from_slice([1, 1]);     
//     let domain_8: Domain<8> = comptime { Domain::<8>::new() };
//     let (quotient, remainder) = dividend.div::<2, 8>(divisor, domain_8);
    
//     // Expected Q: 1+x (deg 1). Capacity N1=4 (from div signature).
//     // Expected R: 0 (deg 0). Capacity N2=2 (from div signature).
//     let expected_quotient: Polynomial<4> = Polynomial::new_from_slice([1, 1]); 
//     let expected_remainder: Polynomial<2> = Polynomial::<2>::zero();
    
//     assert(quotient == expected_quotient);
//     assert(remainder == expected_remainder);

//     //Test case 2: Division with remainder (x^3 + 2x^2 + 3x + 4) / (x^2 + 1) = x + 2, remainder 2x + 2
//     // Dividend: 4+3x+2x^2+x^3 (deg 3). Capacity N1=4.
//     // Divisor: 1+0x+x^2 (deg 2). Capacity N2=4.
//     // N1 >= N2 (4 >= 4).
//     // Domain for mul check: L = N1+N1 = 4+4=8.
//     let dividend2: Polynomial<4> = Polynomial::new_from_slice([4, 3, 2, 1]); 
//     let divisor2: Polynomial<4> = Polynomial::new_from_slice([1, 0, 1]);     
//     // Can reuse domain_8
//     let (quotient2, remainder2) = dividend2.div::<4, 8>(divisor2, domain_8);
    
//     // Expected Q: 2+x (deg 1). Capacity N1=4.
//     // Expected R: 2+2x (deg 1). Capacity N2=4.
//     let expected_quotient2: Polynomial<4> = Polynomial::new_from_slice([2, 1]); 
//     let expected_remainder2: Polynomial<4> = Polynomial::new_from_slice([2, 2]);
    
//     assert(quotient2 == expected_quotient2);
//     assert(remainder2 == expected_remainder2);

//     // Test case 3: Division by constant (2x^2 + 4x + 6) / 2 = x^2 + 2x + 3
//     // Dividend: 6+4x+2x^2 (deg 2). Capacity N1=4.
//     // Divisor: 2 (deg 0). Capacity N2=1.
//     // N1 >= N2 (4 >= 1).
//     // Domain for mul check: L = N1+N1 = 4+4=8.
//     let dividend3: Polynomial<4> = Polynomial::new_from_slice([6, 4, 2]); 
//     let divisor3: Polynomial<1> = Polynomial::new_from_slice([2]);        
//     // Can reuse domain_8
//     let (quotient3, remainder3) = dividend3.div::<1, 8>(divisor3, domain_8);
    
//     // Expected Q: 3+2x+x^2 (deg 2). Capacity N1=4.
//     // Expected R: 0 (deg 0). Capacity N2=1.
//     let expected_quotient3: Polynomial<4> = Polynomial::new_from_slice([3, 2, 1]); 
//     let expected_remainder3: Polynomial<1> = Polynomial::<1>::zero();
    
//     assert(quotient3 == expected_quotient3);
//     assert(remainder3 == expected_remainder3);

//     // Test case 4: Dividend smaller degree than divisor
//     // (x+1) / (x^2+1). Q=0, R=x+1
//     // Dividend: 1+x (deg 1). N1=2
//     // Divisor: 1+x^2 (deg 2). N2=4
//     // Need N1 >= N2 for current `div`. So we must expand dividend or use a different div.
//     // Let's make dividend N1=4 by expansion conceptually: 1+x+0x^2+0x^3
//     let dividend4_small: Polynomial<2> = Polynomial::new_from_slice([1,1]);
//     // let dividend4: Polynomial<4> = dividend4_small.expand(); // Assuming expand exists and is used appropriately.
//     // For the test, ensure N1 >= N2.
//     // Let dividend: Polynomial<4> = [1,1,0,0]
//     // Let divisor: Polynomial<4> = [1,0,1,0]
//     let dividend4: Polynomial<4> = Polynomial::new_from_slice([1,1]); // Padded to [1,1,0,0] by new_from_slice
//     let divisor4: Polynomial<4> = Polynomial::new_from_slice([1,0,1]); // Padded to [1,0,1,0]
//     // N1=4, N2=4. N1>=N2 holds. Domain L = 4+4=8.
//     let (quotient4, remainder4) = dividend4.div::<4, 8>(divisor4, domain_8);

//     let expected_quotient4: Polynomial<4> = Polynomial::<4>::zero();
//     let expected_remainder4: Polynomial<4> = Polynomial::new_from_slice([1,1]); // Remainder is 1+x

//     assert(quotient4 == expected_quotient4);
//     // The div implementation returns remainder.reduce<M> where M is original divisor capacity.
//     // Here M (N2) is 4. So remainder should be Polynomial<4>.
//     assert(remainder4 == expected_remainder4);
// }
