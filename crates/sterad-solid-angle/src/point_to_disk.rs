use crate::error::SolidAngleError;
use std::f64::consts::PI;

/// A point to a circular disk, where the point located on the normal passing through the center of the disk.
///
/// # Arguments
///
/// * `h` - Distance between the point and the center of the disk.
/// * `r` - Radius of the circular disk.
///
/// # Returns
///
/// * `Result<f64, SolidAngleError>` - The solid angle subtended by the disk, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `SolidAngleError` if:
///
/// * `h` is less than or equal to 0.
/// * `r` is less than or equal to 0.
///
/// # Examples
///
/// ```
/// use sterad_solid_angle::point_to_disk;
/// use std::f64::consts::PI;
///
/// let h = 4.0;
/// let r = 3.0;
/// let result = point_to_disk::center(h, r).unwrap();
/// assert!((result - 0.4 * PI).abs() < 1e-10);
/// ```
pub fn center(h: f64, r: f64) -> Result<f64, SolidAngleError> {
    if h <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if r <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "r",
            message: "r must be greater than 0".to_string(),
        });
    }

    let omega = 2.0 * PI * (1.0 - h / (h.powi(2) + r.powi(2)).sqrt());
    Ok(omega)
}
