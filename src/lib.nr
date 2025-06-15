pub mod polynomial;
pub mod lagrange_polynomial;
pub mod domain;
pub mod array_utils;
pub mod ROOT_OF_UNITY;

#[test]
fn test_basic_usage_example() {
    // Create polynomials with different sizes
    // Polynomial<8, 3> can hold up to degree 7 polynomials (8 coefficients)
    // new_from_slice allows you to pass any size slice as long as it's smaller than the capacity
    let poly1 = polynomial::Polynomial::<8, 3>::new_from_slice([1, 2, 3]); // 3x^2 + 2x + 1
    let poly2 = polynomial::Polynomial::<4, 2>::new_from_slice([4, 5]);    // 5x + 4
    
    // Basic operations
    let sum = poly1.add(poly2);                    // 3x^2 + 7x + 5 // You can use different sizes polys as long
    let diff = poly1.sub(poly2);                   // 3x^2 - 3x - 3 // as the bigger capacity poly is first
    let scaled = poly1.mul_scalar(2);              // 6x^2 + 4x + 2
    let negated = poly1.neg();                     // -3x^2 - 2x - 1
    
    // Evaluate at a point
    let value = poly1.evaluate(2);                 // Returns 17
    
    // Create a domain for NTT-based operations
    let domain = comptime { domain::Domain::<16, 4>::new() }; //you only need to create it once with the max size
    
    // Multiplication (returns larger polynomial)
    let product = poly1.mul(poly2, domain);        // Returns Polynomial<16, 4>
    
    // Division with quotient and remainder
    let dividend = polynomial::Polynomial::<8, 3>::new_from_slice([1, 2, 1]);  // x^2 + 2x + 1
    let divisor = polynomial::Polynomial::<8, 3>::new_from_slice([1, 1]);      // x + 1
    let (quotient, remainder) = dividend.div(divisor, domain);
    
    // Convert between different sizes
    let larger_poly = poly1.expand::<16, 4>();     // Convert to larger capacity
    let smaller_poly = larger_poly.reduce::<8, 3>(); // Convert back (must have trailing zeros)
    
    // Lagrange basis polynomials (NTT form)
    let lagrange_poly1 = poly1.ntt(domain);        // Convert to Lagrange basis
    let poly2_expanded = poly2.expand::<8, 3>();
    let lagrange_poly2 = poly2_expanded.ntt(domain); // Expand poly2 first, then convert
    
    // Operations in Lagrange basis (often faster)
    let lagrange_sum = lagrange_poly1.add(lagrange_poly2);
    let lagrange_product = lagrange_poly1.mul(lagrange_poly2);
    
    // Convert back to coefficient form
    let recovered_sum = lagrange_sum.intt(domain);
    let recovered_product = lagrange_product.intt(domain);
    
    // Modular multiplication (mod x^N - 1, where N is the first polynomial's capacity)
    let mod_product = poly1.mul_mod(poly2_expanded, domain); // Result has same size as poly1
    
    // Verify some results
    assert(poly1.evaluate(2) == 17);
    assert(sum.get_coefficients()[0] == 5); // constant term should be 5
    assert(sum.get_coefficients()[1] == 7); // x term should be 7
    assert(sum.get_coefficients()[2] == 3); // x^2 term should be 3
    
    // Verify basic arithmetic properties
    assert(value == 17); // 3*4 + 2*2 + 1 = 12 + 4 + 1 = 17
    
    // Verify negation
    let neg_coeffs = negated.get_coefficients();
    assert(neg_coeffs[0] == -1);
    assert(neg_coeffs[1] == -2);
    assert(neg_coeffs[2] == -3);
    
    // Verify scaling
    let scaled_coeffs = scaled.get_coefficients();
    assert(scaled_coeffs[0] == 2);
    assert(scaled_coeffs[1] == 4);
    assert(scaled_coeffs[2] == 6);
}