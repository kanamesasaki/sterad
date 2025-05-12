#[cfg(feature = "double_precision")]
pub type Float = f64;
#[cfg(not(feature = "double_precision"))]
pub type Float = f32;

pub mod float {
    use super::Float;

    pub const FLOAT_INFINITY: Float = Float::INFINITY;
    pub const FLOAT_MACHINE_EPSILON: Float = Float::EPSILON * 0.5;

    #[cfg(feature = "double_precision")]
    pub const FLOAT_ONE_MINUS_EPSILON: Float = Float::from_bits(0x3fefffffffffffff);
    #[cfg(not(feature = "double_precision"))]
    pub const FLOAT_ONE_MINUS_EPSILON: Float = Float::from_bits(0x3f7fffff);

    #[cfg(feature = "double_precision")]
    pub const PI: Float = std::f64::consts::PI;
    #[cfg(not(feature = "double_precision"))]
    pub const PI: Float = std::f32::consts::PI;

    /// Returns the next representable floating-point value greater than `v`.
    ///
    /// This function manipulates the binary representation of the floating-point value
    /// to find the next greater value in the IEEE-754 representation.
    ///
    /// # Special Cases
    ///
    /// - If `v` is positive infinity, returns positive infinity
    /// - If `v` is negative zero (-0.0), returns positive zero (0.0)
    ///
    #[inline]
    pub fn next_float_up(v: Float) -> Float {
        if v.is_infinite() && v > 0.0 {
            return v;
        }
        if v == -0.0 {
            return 0.0;
        }
        let mut bits = v.to_bits();
        if v > 0.0 {
            bits += 1;
        } else {
            bits -= 1;
        }
        Float::from_bits(bits)
    }

    /// Returns the next representable floating-point value less than `v`.
    ///
    /// This function manipulates the binary representation of the floating-point value
    /// to find the next smaller value in the IEEE-754 representation.
    ///
    /// # Special Cases
    ///
    /// - If `v` is negative infinity, returns negative infinity
    /// - If `v` is positive zero (0.0), returns negative zero (-0.0)
    ///
    #[inline]
    pub fn next_float_down(v: Float) -> Float {
        if v.is_infinite() && v < 0.0 {
            return v;
        }
        if v == 0.0 {
            return -0.0;
        }
        let mut bits = v.to_bits();
        if v > 0.0 {
            bits -= 1;
        } else {
            bits += 1;
        }
        Float::from_bits(bits)
    }

    #[inline]
    pub fn exponent(v: Float) -> i32 {
        #[cfg(feature = "double_precision")]
        {
            // Convert Float to bit representation
            let bits: u64 = v.to_bits();
            // Shift the exponent part to the least significant bits
            let shifted: u64 = bits >> 52;
            // bitmask to extract the exponent part
            let exponent_bits: i32 = (shifted & 0x7ff) as i32;
            // Subtract the bias
            exponent_bits - 1023
        }
        #[cfg(not(feature = "double_precision"))]
        {
            // Convert Float to bit representation
            let bits: u32 = v.to_bits();
            // Shift the exponent part to the least significant bits
            let shifted: u32 = bits >> 23;
            // bitmask to extract the exponent part
            let exponent_bits: i32 = (shifted & 0xff) as i32;
            // Subtract the bias
            exponent_bits - 127
        }
    }

    #[inline]
    pub fn significand(v: Float) -> u64 {
        #[cfg(feature = "double_precision")]
        {
            v.to_bits() & ((1u64 << 52) - 1)
        }
        #[cfg(not(feature = "double_precision"))]
        {
            v.to_bits() as u64 & ((1u64 << 23) - 1)
        }
    }

    #[inline]
    pub fn sign_bit(v: Float) -> u64 {
        #[cfg(feature = "double_precision")]
        {
            v.to_bits() & 0x8000000000000000
        }
        #[cfg(not(feature = "double_precision"))]
        {
            v.to_bits() as u64 & 0x80000000
        }
    }

    #[inline]
    pub fn add_round_up(a: Float, b: Float) -> Float {
        next_float_up(a + b)
    }

    #[inline]
    pub fn add_round_down(a: Float, b: Float) -> Float {
        next_float_down(a + b)
    }

    #[inline]
    pub fn mul_round_up(a: Float, b: Float) -> Float {
        next_float_up(a * b)
    }

    #[inline]
    pub fn mul_round_down(a: Float, b: Float) -> Float {
        next_float_down(a * b)
    }

    #[inline]
    pub fn div_round_up(a: Float, b: Float) -> Float {
        next_float_up(a / b)
    }

    #[inline]
    pub fn div_round_down(a: Float, b: Float) -> Float {
        next_float_down(a / b)
    }

    #[inline]
    pub fn sqrt_round_up(a: Float) -> Float {
        next_float_up(a.sqrt())
    }

    #[inline]
    pub fn sqrt_round_down(a: Float) -> Float {
        next_float_down(a.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next_float_up() {
        let value = 1.0;
        assert!(float::next_float_up(value) > value);
    }

    #[test]
    fn test_next_float_down() {
        let value = 1.0;
        assert!(float::next_float_down(value) < value);
        assert!(float::next_float_down(value) == float::FLOAT_ONE_MINUS_EPSILON);
    }

    #[test]
    fn test_exponent() {
        let value: Float = 1.0;
        assert_eq!(float::exponent(value), 0);

        let value: Float = 2.0;
        assert_eq!(float::exponent(value), 1);

        let value: Float = 0.5;
        assert_eq!(float::exponent(value), -1);

        let value: Float = 4.0;
        assert_eq!(float::exponent(value), 2);

        let value: Float = 0.25;
        assert_eq!(float::exponent(value), -2);
    }
}
