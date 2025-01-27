use crate::elliptic_integral::{ellk, ellpic};
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

/// Calculates the solid angle subtended by a disk, where the source point is not on the central axis.
///
/// This function computes the solid angle of a disk with radius `rm`.
/// The source point is at height `h` above the disk and at a distance `ro` from the disk central axis.
///
/// # Arguments
/// * `h` - Height of the source point above the disk
/// * `rm` - Radius of the disk (must be positive)
/// * `ro` - Distance between the source point and the disk central axis
///
/// # Returns
/// * `Ok(f64)` - Computed solid angle in steradians
/// * `Err(SolidAngleError)` - If parameters are invalid
///
/// # Errors
/// Returns error if:
/// * h is not positive
/// * rm is not positive
/// * ro is not positive
///
/// # Examples
///
/// ```
/// use sterad_solid_angle::point_to_disk;
///
/// let h = 2.0;
/// let rm = 2.0;
/// let ro = 4.0;
/// let result = point_to_disk::offset(h, rm, ro).unwrap();
/// assert!((result - 0.32580).abs() < 1e-5);
/// ```
pub fn offset(h: f64, rm: f64, ro: f64) -> Result<f64, SolidAngleError> {
    if h <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if rm <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "rm",
            message: "rm must be greater than 0".to_string(),
        });
    }
    if ro <= 0.0 {
        return Err(SolidAngleError::InvalidInput {
            param_name: "ro",
            message: "ro must be greater than 0".to_string(),
        });
    }

    let omega: f64;
    let r_max = ((ro + rm).powi(2) + h.powi(2)).sqrt();
    let alpha2 = 4.0 * ro * rm / (ro + rm).powi(2);
    let k2 = 4.0 * ro * rm / (h.powi(2) + (ro + rm).powi(2));
    let k = k2.sqrt();
    if ro < rm {
        let beta_max = PI;
        let omega_half = beta_max - h / r_max * ellk(k)?
            + h / r_max * (ro - rm) / (ro + rm) * ellpic(alpha2, k)?;
        omega = 2.0 * omega_half;
    } else if ro == rm {
        omega = PI - 2.0 * h / r_max * ellk(k)?;
    } else {
        let omega_half =
            -h / r_max * ellk(k)? + h / r_max * (ro - rm) / (ro + rm) * ellpic(alpha2, k)?;
        omega = 2.0 * omega_half;
    }
    Ok(omega)
}
