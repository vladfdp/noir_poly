# Noir Polynomial Library

A polynomial library implemented in Noir.

## Features

- Core polynomial representation with coefficient-based storage
- Basic arithmetic operations: addition, subtraction, multiplication, division
- Polynomial evaluation at field points

Due to the way Noir works, every polynomial is a fixed length array (because loops must have known size at compile time). The length of that array, which gives you the maximum degree is a global variable that can be found in `src/MAX_DEGREE.nr`. You can modify it for your need


## Usage

```noir
let poly1 = Polynomial::new_from_slice([1, 2, 3]); // 3x^2 + 2x + 1
let poly2 = Polynomial::new_from_slice([4, 5]);    // 5x + 4
let result = poly1.add(poly2);                     // 3x^2 + 7x + 5

// Evaluate at a point
let value = poly1.evaluate(2); // Returns 17
``` 