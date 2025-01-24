use crate::error::SolidAngleError;

/// Calculate the solid angle subtended by a triangle at a given point using the vertices of the triangle.
///
/// # Arguments
///
/// * `r1` - A 3-element array representing the coordinates of the first vertex of the triangle.
/// * `r2` - A 3-element array representing the coordinates of the second vertex of the triangle.
/// * `r3` - A 3-element array representing the coordinates of the third vertex of the triangle.
///
/// # Returns
///
/// * `Result<f64, SolidAngleError>` - The solid angle subtended by the triangle, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `SolidAngleError` if:
///
/// * Any of the vertices (`r1`, `r2`, `r3`) are at the origin `[0.0, 0.0, 0.0]`.
///
/// # Examples
///
/// ```
/// use sterad_solid_angle::point_to_triangle;
/// use std::f64::consts::PI;
///
/// let r1 = [1.0, 0.0, 0.0];
/// let r2 = [0.0, 1.0, 0.0];
/// let r3 = [0.0, 0.0, 1.0];
/// let result = point_to_triangle::vertices(r1, r2, r3).unwrap();
/// assert!((result - 0.5 * PI).abs() < 1e-10);
/// ```
pub fn vertices(r1: [f64; 3], r2: [f64; 3], r3: [f64; 3]) -> Result<f64, SolidAngleError> {
    if r1 == [0.0, 0.0, 0.0] || r2 == [0.0, 0.0, 0.0] || r3 == [0.0, 0.0, 0.0] {
        return Err(SolidAngleError::InvalidInput {
            param_name: "r1, r2, r3",
            message: "r1, r2, r3 must not be zero".to_string(),
        });
    }

    let r1_norm = norm(&r1);
    let r2_norm = norm(&r2);
    let r3_norm = norm(&r3);
    let volume = dot(&r1, &cross(&r2, &r3));
    let dot_r1_r2 = dot(&r1, &r2);
    let dot_r2_r3 = dot(&r2, &r3);
    let dot_r3_r1 = dot(&r3, &r1);
    let omega = 2.0
        * (volume
            / (r1_norm * r2_norm * r3_norm
                + dot_r1_r2 * r3_norm
                + dot_r2_r3 * r1_norm
                + dot_r3_r1 * r2_norm))
            .atan();
    Ok(omega)
}

fn cross(v1: &[f64; 3], v2: &[f64; 3]) -> [f64; 3] {
    [
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ]
}

fn dot(v1: &[f64; 3], v2: &[f64; 3]) -> f64 {
    v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
}

fn norm(v: &[f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
