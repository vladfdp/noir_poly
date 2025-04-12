# Noir Polynomial Library

A polynomial library implemented in Noir.

## Features

- Core polynomial representation with coefficient-based storage
- Basic arithmetic operations: addition, subtraction, multiplication, division
- Polynomial evaluation at field points

Due to the way Noir works, every polynomial is a fixed length array (because loops must have known size at compile time). The length of that array, which gives you the maximum degree is a global variable that can be found in `src/MAX_DEGREE.nr`. You can modify it for your need


## Usage

```noir
// Create polynomials - coefficients are ordered from constant term to highest degree
let poly1 = Polynomial::new_from_slice([1, 2, 3]); // 3x^2 + 2x + 1
let poly2 = Polynomial::new_from_slice([4, 5]);    // 5x + 4

// Basic operations
let sum = poly1.add(poly2);                     // 3x^2 + 7x + 5
let diff = poly1.sub(poly2);                    // 3x^2 - 3x - 3
let product = poly1.mul(poly2);                 // 15x^3 + 23x^2 + 13x + 4
let scaled = poly1.mul_scalar(2);               // 6x^2 + 4x + 2

// Evaluate at a point
let value = poly1.evaluate(2);                  // Returns 17

// Special polynomials
let one = Polynomial::one();                    // 1
let zero = Polynomial::zero();                  // 0
let vanishing = Polynomial::vanishing_polynomial(3); // (x - 3)

// Division - returns (quotient, remainder)
let dividend = Polynomial::new_from_slice([1, 2, 1]);  // x^2 + 2x + 1
let divisor = Polynomial::new_from_slice([1, 1]);      // x + 1
let (quotient, remainder) = dividend.div(divisor);     // quotient = (x + 1), remainder = 0
``` 