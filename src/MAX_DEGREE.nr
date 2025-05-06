// This file is used to define the maximum degree of the polynomial

pub global MAX_DEGREE_LOG: u8 = 7;
pub global MAX_DEGREE: u32 = (1 << MAX_DEGREE_LOG) - 1;
pub global NORMALIZER: Field = 1 / (MAX_DEGREE as Field + 1);

