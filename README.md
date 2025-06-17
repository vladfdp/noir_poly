# Noir Polynomial Library

A polynomial library implemented in Noir with support for polynomials of different sizes.

## Features

- **Parameterized polynomial sizes**: Create polynomials with different maximum degrees using compile-time parameters
- Core polynomial coefficient representation with constant factor first
- Polynomials in Lagrange basis
- Basic arithmetic operations: addition, subtraction, multiplication, division

- NTT and INTT transformations for fast polynomial operations
- Polynomial division with quotient and remainder
- Size conversion operations (expand/reduce)

The library uses parameterized types `Polynomial<SIZE, LOG_SIZE>` where `SIZE` is the maximum number of coefficients and `LOG_SIZE = log₂(SIZE)`. This allows you to work with polynomials of different capacities in the same program.

## Basic Usage

```noir

    // Create polynomials with different sizes
    // Polynomial<8, 3> can hold up to degree 7 polynomials (8 coefficients)
    // new_from_slice allows you to pass any size slice as long as it's smaller than the capacity
    let poly1 = Polynomial::<8, 3>::new_from_slice([1, 2, 3]); // 3x² + 2x + 1
    let poly2 = Polynomial::<4, 2>::new_from_slice([4, 5]);    // 5x + 4
    
    // Basic operations
    let sum = poly1.add(poly2);                    // 3x² + 7x + 5 // You can use different sizes polys as long
    let diff = poly1.sub(poly2);                   // 3x² - 3x - 3 // as the bigger capacity poly is first
    let scaled = poly1.mul_scalar(2);              // 6x² + 4x + 2
    let negated = poly1.neg();                     // -3x² - 2x - 1
    
    // Evaluate at a point
    let value = poly1.evaluate(2);                 // Returns 17
    
    // Create a domain for NTT-based operations
    let domain = comptime { Domain::<16, 4>::new() }; //you only need to create it once with the max size, it's required for mul and div
    
    // Multiplication (returns larger polynomial)
    let product = poly1.mul(poly2, domain);        // Returns Polynomial<16, 4>
    
    // Division with quotient and remainder
    let dividend = Polynomial::<8, 3>::new_from_slice([1, 2, 1]);  // x² + 2x + 1
    let divisor = Polynomial::<8, 3>::new_from_slice([1, 1]);      // x + 1
    let (quotient, remainder) = dividend.div(divisor, domain);
    
    // Convert between different sizes
    let larger_poly = poly1.expand::<16, 4>();     // Convert to larger capacity
    let smaller_poly = larger_poly.reduce::<8, 3>(); // Convert back (must have trailing zeros)
    
    // Lagrange basis polynomials
    let lagrange_poly1 = poly1.ntt(domain);        // Convert to Lagrange basis
    let lagrange_poly2 = poly2.expand::<8, 3>().ntt(domain); // Expand poly2 first, then convert
    
    // Operations in Lagrange basis
    let lagrange_sum = lagrange_poly1.add(lagrange_poly2);
    let lagrange_product = lagrange_poly1.mul(lagrange_poly2);
    
    // Convert back to coefficient form
    let recovered_sum = lagrange_sum.intt(domain);
    let recovered_product = lagrange_product.intt(domain);
    
    // Modular multiplication (mod x^N - 1, where N is the first polynomial's capacity)
    let mod_product = poly1.mul_mod(poly2, domain); // Result has same size as poly1

```

## Library Capabilities

- **Arithmetic**: `add()`, `sub()`, `mul()`, `mul_scalar()`, `neg()`
- **Evaluation**: `evaluate(x)` - evaluate polynomial at a given point
- **Special polynomials**: `one()`, `zero()`, `vanishing_polynomial(a)`
- **Polynomial division**: `div()` - returns quotient and remainder
- **NTT transformations**: `ntt()` and `intt()` for fast operations
- **Modular multiplication**: `mul_mod()` for multiplication mod (x^N - 1 with N the capacity of the first poly)
- **Size conversion**: `expand()` and `reduce()` to change polynomial capacity

### Requirements
- Domain size must be a power of 2 and specified at compile time
- Polynomial sizes must be a power of 2 and specified at compile time
- LOG_SIZE must equal log₂(SIZE) for a Polynomial<SIZE,LOG_SIZE>, same for lagrange polys and domains
- Operations requiring NTT need a domain with sufficient capacity 
