use crate::error::SolidAngleError;

/// Solid angle subtended by a rectangle at a point located on the plane normal from the corner of the rectangle.
///
/// # Arguments
///
/// * `h` - Distance between the point and the plane of the rectangle.
/// * `a` - Length of one side of the rectangle.
/// * `b` - Length of the other side of the rectangle.
///
/// # Returns
///
/// * `Result<f64, SolidAngleError>` - The solid angle subtended by the rectangle, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `SolidAngleError` if:
///
/// * `h` is less than or equal to 0.
/// * `a` is less than or equal to 0.
/// * `b` is less than or equal to 0.
///
/// # Examples
///
/// ```
/// use sterad_solid_angle::point_to_rectangle;
/// use std::f64::consts::PI;
///
/// let h = 1.0;
/// let a = 1.0;
/// let b = 1.0;
/// let result = point_to_rectangle::corner(h, a, b).unwrap();
/// assert!((result - PI / 6.0).abs() < 1e-10);
/// ```
pub fn corner(h: f64, a: f64, b: f64) -> Result<f64, SolidAngleError> {
    if h <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if a <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "a",
            message: "a must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "b",
            message: "b must be greater than 0".to_string(),
        });
    }

    let omega = (a * b / (h * (a.powi(2) + b.powi(2) + h.powi(2)).sqrt())).atan();
    Ok(omega)
}
