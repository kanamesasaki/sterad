use crate::error::ViewFactorError;
use std::f64::consts::PI;

/// B-3: Differential planar element to finite parallel rectangle. Normal to element passes through corner of rectangle.
///
/// # Arguments
///
/// * `a` - Length of one side of the rectangle.
/// * `b` - Length of the other side of the rectangle.
/// * `c` - Distance between the planar element and the plane of the rectangle.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the rectangle from the planar element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `a` is less than or equal to 0.
/// * `b` is less than or equal to 0.
/// * `c` is less than or equal to 0.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_rectangle;
///
/// let a = 1.0;
/// let b = 1.0;
/// let c = 1.0;
/// let result = diff_element_to_rectangle::parallel_corner(a, b, c).unwrap();
/// assert!((result - 0.13853160599489298).abs() < 1e-10);
/// ```
pub fn parallel_corner(a: f64, b: f64, c: f64) -> Result<f64, ViewFactorError> {
    if a <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a",
            message: "a must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "b",
            message: "b must be greater than 0".to_string(),
        });
    }
    if c <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "c",
            message: "c must be greater than 0".to_string(),
        });
    }

    let a_over_c = a / c;
    let b_over_c = b / c;
    let vf = 1.0 / (2.0 * PI)
        * (a_over_c / (1.0 + a_over_c.powi(2)).sqrt()
            * (b_over_c / (1.0 + a_over_c.powi(2)).sqrt()).atan()
            + b_over_c / (1.0 + b_over_c.powi(2)).sqrt()
                * (a_over_c / (1.0 + b_over_c.powi(2)).sqrt()).atan());

    Ok(vf)
}

/// B-4: Differential planar element to rectangle in plane 90° to plane of element and perpendicular to corner of plane.
///
/// # Arguments
///
/// * `a` - Length of one side of the rectangle.
/// * `b` - Length of the other side of the rectangle.
/// * `c` - Distance between the planar element and the plane of the rectangle.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the rectangle from the planar element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `a` is less than or equal to 0.
/// * `b` is less than or equal to 0.
/// * `c` is less than or equal to 0.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_rectangle;
///
/// let a = 1.0;
/// let b = 1.0;
/// let c = 1.0;
/// let result = diff_element_to_rectangle::perpendicular_corner(a, b, c).unwrap();
/// assert!((result - 0.05573419700255351).abs() < 1e-10);
/// ```
pub fn perpendicular_corner(a: f64, b: f64, c: f64) -> Result<f64, ViewFactorError> {
    if a <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a",
            message: "a must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "b",
            message: "b must be greater than 0".to_string(),
        });
    }
    if c <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "c",
            message: "c must be greater than 0".to_string(),
        });
    }

    let a_over_b = a / c;
    let c_over_b = b / c;
    let y = (a_over_b.powi(2) + c_over_b.powi(2)).sqrt();
    let vf = 1.0 / (2.0 * PI) * ((1.0 / c_over_b).atan() - c_over_b / y * (1.0 / y).atan());
    Ok(vf)
}
