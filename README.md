# Noir Polynomial Library

A polynomial library implemented in Noir.

## Features

- Core polynomial coefficient representation with constant factor first
- Basic arithmetic operations: addition, subtraction, multiplication, division
- Polynomial evaluation at field points
- Lagrange Basis Polynomials, basically just arrays of a polynomials evaluation on the roots of unity
- NTT and INTT to go from coefficient basis to lagrange basis and back

Due to the way Noir works, every polynomial is a fixed length array (because loops must have known size at compile time). The length of that array, which gives you the maximum degree is a global variable that can be found in `src/MAX_DEGREE.nr`. You can modify it for your need


## Usage

```noir
// Create polynomials - coefficients are ordered from constant term to highest degree
let poly1 = Polynomial::new_from_slice([1, 2, 3]); // 3x^2 + 2x + 1
let poly2 = Polynomial::new_from_slice([4, 5]);    // 5x + 4

// Special polynomials
let one = Polynomial::one();                    // 1
let zero = Polynomial::zero();                  // 0
let vanishing = Polynomial::vanishing_polynomial(3); // (x - 3)

// Basic operations
let sum = poly1.add(poly2);                     // 3x^2 + 7x + 5
let diff = poly1.sub(poly2);                    // 3x^2 - 3x - 3
let scaled = poly1.mul_scalar(2);               // 6x^2 + 4x + 2

// Evaluate at a point
let value = poly1.evaluate(2);                  // Returns 17

// You need to create a domain for multiplication and division as these operations use ntt/intt
// The domain is created at compile time.
let domain = comptime { Domain::new() };

// Multiplication
let product = poly1.mul(poly2, domain);                 // 15x^3 + 23x^2 + 13x + 4

// Division - returns (quotient, remainder)
let dividend = Polynomial::new_from_slice([1, 2, 1]);  // x^2 + 2x + 1
let divisor = Polynomial::new_from_slice([1, 1]);      // x + 1
let (quotient, remainder) = dividend.div(divisor, domain);     // quotient = (x + 1), remainder = 0

// NTT give the evaluations of the polynomial on the domain elements.
let lb_poly1 = poly1.ntt(domain);
let lb_poly2 = poly2.ntt(domain);

// Note that this still represents the same polynomial, just in a different form.
// you can operate on them similarly to how you would with the coefficient form.
let lb_sum = lb_poly1.add(lb_poly2);
let lb_product = lb_poly1.mul(lb_poly2);

// INTT will give you back the poly in coefficient form.
let recovered_sum = lb_sum.intt(domain);
let recovered_product = lb_product.intt(domain);

assert(recovered_sum.coefficients == sum.coefficients);
assert(recovered_product.coefficients == product.coefficients);
``` 