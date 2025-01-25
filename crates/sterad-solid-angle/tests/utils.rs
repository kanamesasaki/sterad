// Divide `start` and `end` into `num` equal parts and return the middle points as Vector.
fn linspace_middle(start: f64, end: f64, num: usize) -> Vec<f64> {
    let step = (end - start) / num as f64;
    let start_middle = start + step / 2.0;
    (0..num).map(|i| start_middle + i as f64 * step).collect()
}

/// Calculate the numerical view factor of an ellipse from a differential plate.
/// - the differential plate is located at (0, 0, 0)
/// - the ellipse is given by: (x-xc)^2/a^2 + (y-yc)^2/b^2 = 1, z = zc
///
/// # Arguments
///
/// * `a` - Semi-major axis of the ellipse.
/// * `b` - Semi-minor axis of the ellipse.
/// * `theta` - Polar angle of the differential planar element.
/// * `phi` - Azimuthal angle of the differential planar element.
/// * `xc` - X-coordinate of the center of the ellipse.
/// * `yc` - Y-coordinate of the center of the ellipse.
/// * `zc` - Z-coordinate of the center of the ellipse.
/// * `a_num` - Number of divisions along the semi-major axis.
/// * `b_num` - Number of divisions along the semi-minor axis.
///
/// # Returns
///
/// * `f64` - The numerical view factor of the ellipse from the differential planar element.
///
pub fn point_to_ellipse_numerical(
    a: f64,
    b: f64,
    xc: f64,
    yc: f64,
    zc: f64,
    a_num: usize,
    b_num: usize,
) -> f64 {
    let vec_a = linspace_middle(-a + xc, a + xc, a_num);
    let vec_b = linspace_middle(-b + yc, b + yc, b_num);
    let dsize = 4.0 * a * b / (a_num * b_num) as f64;
    let a2 = a.powi(2);
    let b2 = b.powi(2);

    // direction of the mesh element is always (0, 0, -1)

    let mut omega: f64 = 0.0;
    for x in &vec_a {
        for y in &vec_b {
            if (x - xc).powi(2) / a2 + (y - yc).powi(2) / b2 > 1.0 {
                continue;
            }
            let s2 = x.powi(2) + y.powi(2) + zc.powi(2);
            let s = s2.sqrt(); // distance from the diff_plate to the mesh
            let plate_to_mesh_z = zc / s;
            let cos_theta = plate_to_mesh_z;
            let omega_mesh = dsize * cos_theta / s2;
            omega += omega_mesh.max(0.0);
        }
    }
    omega
}
