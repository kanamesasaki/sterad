use crate::error::SolidAngleError;

/// Calculates Carlson's elliptic integral of the first kind, rf(x, y, z).
/// rf(x, y, z) = 1/2 * int_{0}^{infty} dt / sqrt((t+x) * (t+y) * (t+z))
///
/// # Arguments
///
/// * `x` - The first input parameter, a non-negative floating-point number.
/// * `y` - The second input parameter, a non-negative floating-point number.
/// * `z` - The third input parameter, a non-negative floating-point number.
///
/// # Returns
///
/// The value of the Carlson's elliptic integral of the first kind, RF(x, y, z),
/// as a floating-point number.
///
/// # Errors
///
/// * If any of the input values is negative.
/// * If the sum of any two input values is less than TINY.
/// * If any of the input values is greater than BIG.
///
/// # Example
///
/// ```
/// use sterad_solid_angle::elliptic_integral::rf;
///
/// let x = 1.0;
/// let y = 2.0;
/// let z = 3.0;
/// let result = rf(x, y, z).unwrap();
/// assert!((result - 0.726946).abs() < 1e-6);
/// ```
pub fn rf(x: f64, y: f64, z: f64) -> Result<f64, SolidAngleError> {
    let error_tolerance = 0.0025;
    let third = 1.0 / 3.0;
    let c1 = 1.0 / 24.0;
    let c2 = 0.1;
    let c3 = 3.0 / 44.0;
    let c4 = 1.0 / 14.0;
    let tiny = 5.0 * f64::MIN_POSITIVE;
    let big = 0.2 * f64::MAX;

    if x < 0.0 || y < 0.0 || z < 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "x, y, and z must be 0 or greater".to_string(),
        });
    }
    if x + y < tiny || x + z < tiny || y + z < tiny {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "Two inputs are too close to zero".to_string(),
        });
    }
    if x > big || y > big || z > big {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "One of the inputs is too large".to_string(),
        });
    }

    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut delx: f64;
    let mut dely: f64;
    let mut delz: f64;
    let mut ave: f64;
    let mut counter: i32 = 0;

    loop {
        let sqrtx = xt.sqrt();
        let sqrty = yt.sqrt();
        let sqrtz = zt.sqrt();
        let alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        ave = third * (xt + yt + zt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;

        counter += 1;
        if counter > 1000 {
            return Err(SolidAngleError::RuntimeError {
                message: "Calculation of carlson_1st_kind did not converge".to_string(),
            });
        }

        if delx.abs().max(dely.abs()).max(delz.abs()) <= error_tolerance {
            break;
        }
    }

    let e2 = delx * dely - delz * delz;
    let e3 = delx * dely * delz;

    let result = (1.0 + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3) / ave.sqrt();
    Ok(result)
}

/// Calculates the Legendre form elliptic integral of the first kind, ellf(phi, k).
/// ellf(phi, k) = int_{0}^{phi} dtheta / sqrt(1 - k^2 sin^2(theta))
///
/// # Arguments
///
/// * `phi` - The amplitude, a floating-point number in radians.
/// * `k` - The modulus, a floating-point number between -1 and 1 (inclusive).
///
/// # Returns
///
/// * `Ok(f64)` - The value of the elliptic integral of the first kind, ellf(phi, k),
///   as a floating-point number, if the input is valid.
/// * `Err(SolidAngleError)` - An error if the input `k` is not within the valid range.
///
/// # Example
///
/// ```
/// use sterad_solid_angle::elliptic_integral::ellf;
///
/// let phi: f64 = 0.3;
/// let k2: f64 = 0.8;
/// let result = ellf(phi, k2.sqrt()).unwrap();
/// assert!((result - 0.303652).abs() < 1e-6);
/// ```
pub fn ellf(phi: f64, k: f64) -> Result<f64, SolidAngleError> {
    if !(-1.0..=1.0).contains(&k) {
        return Err(SolidAngleError::InvalidInput {
            param_name: "k",
            message: "k must be between -1 and 1".to_string(),
        });
    }

    let x = phi.cos().powi(2);
    let y = 1.0 - k.powi(2) * phi.sin().powi(2);
    let result = phi.sin() * rf(x, y, 1.0)?;
    Ok(result)
}

/// Calculates the complete elliptic integral of the first kind, ellk(k).
/// ellik(k) = int_{0}^{pi/2} dtheta / sqrt(1 - k^2 sin^2(theta))
/// This is equivalent to rf(0, 1-k^2, 1).
///
/// # Arguments
///
/// * `k` - The modulus, a floating-point number between -1 and 1 (inclusive).
///
/// # Returns
///
/// * `Ok(f64)` - The value of the complete elliptic integral of the first kind, K(k),
///   as a floating-point number, if the input is valid.
/// * `Err(SolidAngleError)` - An error if the input `k` is not within the valid range.
///
/// # Example
///
/// ```
/// use sterad_solid_angle::elliptic_integral::ellk;
///
/// let k2: f64 = 0.5;
/// let result = ellk(k2.sqrt()).unwrap();
/// assert!((result - 1.85407).abs() < 1e-5);
/// ```
pub fn ellk(k: f64) -> Result<f64, SolidAngleError> {
    if !(-1.0..=1.0).contains(&k) {
        return Err(SolidAngleError::InvalidInput {
            param_name: "k",
            message: "k must be between -1 and 1".to_string(),
        });
    }

    rf(0.0, 1.0 - k.powi(2), 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rf() {
        let x = 1.0;
        let y = 2.0;
        let z = 3.0;
        let expected = 0.726946;
        let result = rf(x, y, z).unwrap();
        assert!(
            (result - expected).abs() < 1e-6,
            "result: {}, expected: {}",
            result,
            expected
        );
    }

    #[test]
    fn test_ellf() {
        let phi: f64 = 0.3;
        let k2: f64 = 0.8;
        let expected = 0.303652;
        let result = ellf(phi, k2.sqrt()).unwrap();
        assert!(
            (result - expected).abs() < 1e-6,
            "result: {}, expected: {}",
            result,
            expected
        );
    }

    #[test]
    fn test_ellk() {
        let k2: f64 = 0.5;
        let expected = 1.85407;
        let result = ellk(k2.sqrt()).unwrap();
        assert!(
            (result - expected).abs() < 1e-5,
            "result: {}, expected: {}",
            result,
            expected
        );
    }
}
