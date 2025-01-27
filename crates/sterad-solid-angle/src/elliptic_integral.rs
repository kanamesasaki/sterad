use crate::error::SolidAngleError;
use std::f64::consts::PI;

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

/// Computes Carlson's elliptic integral of the third kind rj(x,y,z,p)
///
/// The function calculates the value of the integral defined as:
/// rj(x, y, z, p) = 3/2 * int_{0}^{infty} dt / ((t+p) sqrt((t+x) * (t+y) * (t+z)))
/// This funciton is currently implemented only for p > 0 case
///
/// # Arguments
/// * `x` - First parameter, must be non-negative
/// * `y` - Second parameter, must be non-negative
/// * `z` - Third parameter, must be non-negative
/// * `p` - Fourth parameter, must be positive
///
/// # Returns
/// * `Ok(f64)` - The computed value of the elliptic integral
/// * `Err(SolidAngleError)` - If any input parameter is invalid or computation fails
///
/// # Errors
/// Returns error if:
/// * Any of x, y, z, p is negative
/// * If p equals 0
/// * If numerical computation fails to converge
///
/// # Example
/// ```
/// use sterad_solid_angle::elliptic_integral::rj;
/// let x = 2.0;
/// let y = 3.0;
/// let z = 5.0;
/// let p = 7.2;
/// let result = rj(x, y, z, p).unwrap();
/// assert!((result - 0.104459).abs() < 1e-6);
/// ```
pub fn rj(x: f64, y: f64, z: f64, p: f64) -> Result<f64, SolidAngleError> {
    let tiny: f64 = (5.0 * f64::MIN_POSITIVE).cbrt();
    let big: f64 = 0.3 * (0.2 * f64::MAX).cbrt();
    if x < 0.0 || y < 0.0 || z < 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "x, y, and z must be 0 or greater".to_string(),
        });
    }
    if p <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "p",
            message: "p must be greater than zero".to_string(),
        });
    }
    if x + y < tiny || x + z < tiny || y + z < tiny {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "Two inputs are zero or too close to zero".to_string(),
        });
    }
    if x > big || y > big || z > big {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x, y, z",
            message: "One of the inputs is too large".to_string(),
        });
    }

    const ERRTOL: f64 = 0.0015;
    const C1: f64 = 3.0 / 14.0;
    const C2: f64 = 1.0 / 3.0;
    const C3: f64 = 3.0 / 22.0;
    const C4: f64 = 3.0 / 26.0;
    const C5: f64 = 0.75 * C3;
    const C6: f64 = 1.5 * C4;
    const C7: f64 = 0.5 * C2;
    const C8: f64 = C3 + C3;

    let mut fac = 1.0;
    let mut sum = 0.0;

    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut pt = p;

    let mut delx: f64;
    let mut dely: f64;
    let mut delz: f64;
    let mut delp: f64;
    let mut ave: f64;
    let mut counter: i32 = 0;

    loop {
        let sqrtx = xt.sqrt();
        let sqrty = yt.sqrt();
        let sqrtz = zt.sqrt();
        let alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
        let alpha = (pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz).powi(2);
        let beta = pt * (pt + alamb).powi(2);
        sum += fac * rf(alpha, beta, beta)?;
        fac *= 0.25;
        xt = 0.25 * (xt + alamb);
        yt = 0.25 * (yt + alamb);
        zt = 0.25 * (zt + alamb);
        pt = 0.25 * (pt + alamb);
        ave = 0.2 * (xt + yt + zt + pt + pt);
        delx = (ave - xt) / ave;
        dely = (ave - yt) / ave;
        delz = (ave - zt) / ave;
        delp = (ave - pt) / ave;

        counter += 1;
        if counter > 1000 {
            return Err(SolidAngleError::RuntimeError {
                message: "Calculation of carlson_3rd_kind did not converge".to_string(),
            });
        }

        if delx.abs() <= ERRTOL
            || dely.abs() <= ERRTOL
            || delz.abs() <= ERRTOL
            || delp.abs() <= ERRTOL
        {
            break;
        }
    }

    let ea = delx * (dely + delz) + dely * delz;
    let eb = delx * dely * delz;
    let ec = delp * delp;
    let ed = ea - 3.0 * ec;
    let ee = eb + 2.0 * delp * (ea - ec);
    let ans = 3.0 * sum
        + fac
            * (1.0
                + ed * (-C1 + C5 * ed - C6 * ee)
                + eb * (C7 + delp * (-C8 + delp * C4))
                + delp * ea * (C2 - delp * C3)
                - C2 * delp * ec)
            / (ave * ave.sqrt());

    Ok(ans)
}

/// Computes the incomplete elliptic integral of the third kind ellpii(phi, n, k)
///
/// The function calculates the elliptic integral defined as:
/// ellpii(phi, n, k) = int_{0}^{phi} dθ / ((1 - n sin²θ) sqrt(1 - k² sin²θ))
///
/// # Arguments
/// * `phi` - Amplitude angle in radians, must be in range [0, pi/2]
/// * `n` - Characteristic (parameter), must be less than 1
/// * `k` - Modulus, must be in range [-1, 1]
///
/// # Returns
/// * `Ok(f64)` - Computed value of the elliptic integral
/// * `Err(SolidAngleError)` - If parameters are out of valid ranges
///
/// # Errors
/// Returns error if:
/// * phi is not in [0, π/2]
/// * n is greater than or equal to 1
/// * k is not in [-1, 1]
///
/// # Example
/// ```
/// use sterad_solid_angle::elliptic_integral::ellpi;
/// use std::f64::consts::PI;
///
/// let phi = PI / 2.0;
/// let n: f64 = 0.4;
/// let k2: f64 = 0.6;
/// let result = ellpi(phi, n, k2.sqrt()).unwrap();
/// assert!((result - 2.59092).abs() < 1e-5);
/// ```
pub fn ellpii(phi: f64, n: f64, k: f64) -> Result<f64, SolidAngleError> {
    if !(0.0..=PI / 2.0).contains(&phi) {
        return Err(SolidAngleError::InvalidInput {
            param_name: "phi",
            message: "phi must be between 0 and PI/2".to_string(),
        });
    }
    if !(-1.0..=1.0).contains(&k) {
        return Err(SolidAngleError::InvalidInput {
            param_name: "k",
            message: "k must be between -1 and 1".to_string(),
        });
    }
    if n >= 1.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "n",
            message: "n must be less than 1".to_string(),
        });
    }

    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let q = 1.0 - k.powi(2) * sin_phi.powi(2);
    let result = sin_phi * rf(cos_phi.powi(2), q, 1.0)?
        + n * sin_phi.powi(3) / 3.0 * rj(cos_phi.powi(2), q, 1.0, 1.0 - n * sin_phi.powi(2))?;
    Ok(result)
}

/// Computes the incomplete elliptic integral of the third kind ellpic(n, k)
///
/// The function calculates the elliptic integral defined as:
/// ellpic(n, k) = ellipii(π/2, n, k)
///
/// # Arguments
/// * `n` - Characteristic (parameter), must be less than 1
/// * `k` - Modulus, must be in range [-1, 1]
///
/// # Returns
/// * `Ok(f64)` - Computed value of the elliptic integral
/// * `Err(SolidAngleError)` - If parameters are out of valid ranges
///
/// # Errors
/// Returns error if:
/// * n is greater than or equal to 1
/// * k is not in [-1, 1]
///
/// # Example
/// ```
/// use sterad_solid_angle::elliptic_integral::ellpi;
///
/// let n: f64 = 0.4;
/// let k2: f64 = 0.6;
/// let result = ellpi(n, k2.sqrt()).unwrap();
/// assert!((result - 2.59092).abs() < 1e-5);
/// ```
pub fn ellpic(n: f64, k: f64) -> Result<f64, SolidAngleError> {
    if !(-1.0..=1.0).contains(&k) {
        return Err(SolidAngleError::InvalidInput {
            param_name: "k",
            message: "k must be between -1 and 1".to_string(),
        });
    }
    if n >= 1.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "n",
            message: "n must be less than 1".to_string(),
        });
    }

    let q = 1.0 - k.powi(2);
    let result = rf(0.0, q, 1.0)? + n / 3.0 * rj(0.0, q, 1.0, 1.0 - n)?;
    Ok(result)
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
