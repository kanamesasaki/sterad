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

/// Solid angle subtended by a rectangle at a point.
///
/// # Arguments
///
/// * `h` - Distance between the point and the plane of the rectangle.
/// * `x1` - Distance from the point to one side of the rectangle along the x-axis.
/// * `x2` - Distance from the point to the opposite side of the rectangle along the x-axis.
/// * `y1` - Distance from the point to one side of the rectangle along the y-axis.
/// * `y2` - Distance from the point to the opposite side of the rectangle along the y-axis.
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
/// * `x1` is less than 0.
/// * `x2` is less than or equal to `x1`.
/// * `y1` is less than 0.
/// * `y2` is less than or equal to `y1`.
///
/// # Examples
///
/// ```
/// use sterad_solid_angle::point_to_rectangle;
/// use std::f64::consts::PI;
///
/// let h = 1.0;
/// let x1 = 0.0;
/// let x2 = 1.0;
/// let y1 = 0.0;
/// let y2 = 1.0;
/// let result = point_to_rectangle::offset(h, x1, x2, y1, y2).unwrap();
/// assert!((result - PI / 6.0).abs() < 1e-10);
/// ```
pub fn offset(h: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> Result<f64, SolidAngleError> {
    if h <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if x1 < 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x1",
            message: "x1 must be 0 or greater".to_string(),
        });
    }
    if y1 < 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "y1",
            message: "y1 must be 0 or greater".to_string(),
        });
    }
    if x2 <= x1 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "x2",
            message: "x2 must be greater than x1".to_string(),
        });
    }
    if y2 <= y1 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "y2",
            message: "y2 must be greater than y1".to_string(),
        });
    }

    let omega = (x2 * y2 / (h * (x2.powi(2) + y2.powi(2) + h.powi(2)).sqrt())).atan()
        - (x1 * y2 / (h * (x1.powi(2) + y2.powi(2) + h.powi(2)).sqrt())).atan()
        - (x2 * y1 / (h * (x2.powi(2) + y1.powi(2) + h.powi(2)).sqrt())).atan()
        + (x1 * y1 / (h * (x1.powi(2) + y1.powi(2) + h.powi(2)).sqrt())).atan();
    Ok(omega)
}
