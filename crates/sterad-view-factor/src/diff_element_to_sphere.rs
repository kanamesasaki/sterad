use crate::error::ViewFactorError;
use std::f64::consts::PI;

/// B-39: Plane element to sphere; normal to center of element passes through center of sphere.
///
/// # Arguments
///
/// * `h` - Distance between the planar element and the center of the sphere.
/// * `r` - Radius of the sphere.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the sphere from the planar element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `r` is less than or equal to 0.
/// * `h` is less than or equal to `r`.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_sphere;
///
/// let h = 5.0;
/// let r = 2.0;
/// let result = diff_element_to_sphere::center(h, r).unwrap();
/// assert!((result - 0.16).abs() < 1e-10);
/// ```
pub fn center(h: f64, r: f64) -> Result<f64, ViewFactorError> {
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    if h <= r {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than r".to_string(),
        });
    }

    Ok((r / h).powi(2))
}

fn tilted_full(h: f64, r: f64, theta: f64) -> Result<f64, ViewFactorError> {
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    if h <= r {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than r".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let param_h = h / r;
    Ok(theta.cos() / param_h.powi(2))
}

fn tilted_partial(h: f64, r: f64, theta: f64) -> Result<f64, ViewFactorError> {
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    if h <= r {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than r".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let param_h = h / r;
    let sqrt_h2 = (param_h.powi(2) - 1.0).sqrt();
    let vf = 0.5 - (sqrt_h2 / (param_h * theta.sin())).asin() / PI
        + (theta.cos() * (-sqrt_h2 / theta.tan()).acos()
            - sqrt_h2 * (1.0 - param_h.powi(2) * theta.cos().powi(2)).sqrt())
            / (PI * param_h.powi(2));
    Ok(vf)
}

/// B-43: Arbitrarily oriented differential planar element to a sphere
///
/// # Arguments
///
/// * `h` - Distance between the planar element and the center of the sphere.
/// * `r` - Radius of the sphere.
/// * `theta` - Tilt angle of the planar element in radians.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the sphere from the planar element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `r` is less than or equal to 0.
/// * `h` is less than or equal to `r`.
/// * `theta` is not between 0 and PI.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_sphere;
///
/// let h = 5.0;
/// let r = 2.0;
/// let theta = 0.0;
/// let result = diff_element_to_sphere::tilted(h, r, theta).unwrap();
/// assert!((result - 0.16).abs() < 1e-10);
///
/// ```
pub fn tilted(h: f64, r: f64, theta: f64) -> Result<f64, ViewFactorError> {
    if r <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }
    if h <= r {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than r".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let phi = (r / h).asin();
    let vf: f64 = if theta <= PI / 2.0 - phi {
        tilted_full(h, r, theta)?
    } else if theta >= PI / 2.0 + phi {
        0.0
    } else {
        tilted_partial(h, r, theta)?
    };

    Ok(vf)
}
