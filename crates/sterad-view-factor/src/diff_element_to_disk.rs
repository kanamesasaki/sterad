use crate::error::ViewFactorError;
use std::f64::consts::PI;

/// B-12: Planar element dA1 to circular disk A2 in parallel plane. Normal to element passes through center of disk.
///
/// # Arguments
///
/// * `h` - A floating point number representing the distance between the planar element and the circular disk.
/// * `r` - A floating point number representing the radius of the circular disk.
///
/// # Returns
///
/// The View Factor value.
///
/// # Errors
///
/// Returns an error if `h` or `r` is less than or equal to zero.
///
/// # Example
///
/// ```
/// use sterad_view_factor::diff_element_to_disk;
/// let h: f64 = 2.0;
/// let r: f64 = 1.0;
/// let vf: f64 = diff_element_to_disk::parallel_center(h, r).unwrap();
/// assert!((vf - 0.2).abs() < 1e-8);
/// ```
pub fn parallel_center(h: f64, r: f64) -> Result<f64, ViewFactorError> {
    if h <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    let h_over_r = h / r;
    let f = 1.0 / (h_over_r.powi(2) + 1.0);
    Ok(f)
}

fn tilted_center_partial(h: f64, r: f64, theta: f64) -> f64 {
    let cos_alpha: f64 = -h / (r * theta.tan());
    let sin_alpha: f64 = (1.0 - cos_alpha.powi(2)).sqrt();
    let alpha: f64 = cos_alpha.acos();

    let vf1: f64 = -r * h * theta.sin() * sin_alpha / (PI * (r.powi(2) + h.powi(2)))
        + r.powi(2) * theta.cos() * alpha / (PI * (r.powi(2) + h.powi(2)));

    let vf2: f64 = h * theta.sin() / (PI * (r.powi(2) * cos_alpha.powi(2) + h.powi(2)).sqrt())
        * (r * sin_alpha / (r.powi(2) * cos_alpha.powi(2) + h.powi(2)).sqrt()).atan()
        - r * theta.cos() * cos_alpha / (PI * (r.powi(2) * cos_alpha.powi(2) + h.powi(2)).sqrt())
            * (r * sin_alpha / (r.powi(2) * cos_alpha.powi(2) + h.powi(2)).sqrt()).atan();

    vf1 + vf2
}

fn tilted_center_full(h: f64, r: f64, theta: f64) -> f64 {
    r.powi(2) / (r.powi(2) + h.powi(2)) * theta.cos()
}

/// B-13: Differential tilted planar element dA1 to disk A2. Element lies on normal to disk passing through disk center.
///
/// # Arguments
///
/// * `h` - A floating point number representing the distance between the planar element and the circular disk. Must be greater than 0.
/// * `r` - A floating point number representing the radius of the circular disk. Must be greater than 0.
/// * `theta` - A floating point number representing the tilt angle in radians. Must be between 0 and PI.
///
/// # Returns
///
/// The View Factor value as a `Result<f64, ViewFactorError>`.
///
/// # Errors
///
/// Returns an error if `h` or `r` is less than or equal to zero, or if `theta` is not between 0 and PI.
pub fn tilted_center(h: f64, r: f64, theta: f64) -> Result<f64, ViewFactorError> {
    if h <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let vf: f64 = if theta <= (h / r).atan() {
        tilted_center_full(h, r, theta)
    } else if theta < PI - (h / r).atan() {
        tilted_center_partial(h, r, theta)
    } else {
        0.0
    };

    Ok(vf)
}

/// B-14: Planar element dA1 to a circular disk A2 in a parallel plane. Element is offset from normal to disk center by distance a.
///
/// # Arguments
///
/// * `h` - A floating point number representing the distance between the planar element and the circular disk.
/// * `r` - A floating point number representing the radius of the circular disk.
/// * `offset` - A floating point number representing the distance between the planar element and the normal line passing the center of the circular disk.
///
/// # Returns
///
/// The View Factor value.
///
/// # Errors
///
/// Returns an error if `h` or `r` is less than or equal to zero.
pub fn parallel_offset(h: f64, r: f64, offset: f64) -> Result<f64, ViewFactorError> {
    if h <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    let h_over_a: f64 = h / offset;
    let r_over_a: f64 = r / offset;
    let z: f64 = 1.0 + h_over_a.powi(2) + r_over_a.powi(2);
    let f: f64 =
        0.5 * (1.0 - (z - 2.0 * r_over_a.powi(2)) / f64::sqrt(z.powi(2) - 4.0 * r_over_a.powi(2)));
    Ok(f)
}
