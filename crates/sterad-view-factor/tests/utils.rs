use std::f64::consts::{FRAC_PI_2, PI};

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
#[allow(clippy::too_many_arguments)]
pub fn diff_element_to_ellipse_numerical(
    a: f64,
    b: f64,
    theta: f64,
    phi: f64,
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

    // direction of the differential plate
    let dir_x = theta.sin() * phi.cos();
    let dir_y = theta.sin() * phi.sin();
    let dir_z = theta.cos();

    // direction of the mesh element is always (0, 0, -1)

    let mut vf: f64 = 0.0;
    for x in &vec_a {
        for y in &vec_b {
            if (x - xc).powi(2) / a2 + (y - yc).powi(2) / b2 > 1.0 {
                continue;
            }
            let s2 = x.powi(2) + y.powi(2) + zc.powi(2);
            let s = s2.sqrt(); // distance from the diff_plate to the mesh
            let plate_to_mesh_x = x / s;
            let plate_to_mesh_y = y / s;
            let plate_to_mesh_z = zc / s;
            let cos_theta_1 =
                plate_to_mesh_x * dir_x + plate_to_mesh_y * dir_y + plate_to_mesh_z * dir_z;
            let cos_theta_2 = plate_to_mesh_z;
            let vf_mesh = dsize * cos_theta_1 * cos_theta_2 / (PI * s2);
            vf += vf_mesh.max(0.0);
        }
    }
    vf
}

#[allow(clippy::too_many_arguments)]
pub fn diff_element_to_ellipsoid_numerical(
    a: f64,
    b: f64,
    c: f64,
    theta: f64,
    phi: f64,
    xc: f64,
    yc: f64,
    zc: f64,
    alpha_num: usize,
    beta_num: usize,
) -> f64 {
    let vec_alpha = linspace_middle(0.0, PI, alpha_num);
    let vec_beta = linspace_middle(0.0, 2.0 * PI, beta_num);
    let d_alpha = PI / alpha_num as f64;
    let d_beta = 2.0 * PI / beta_num as f64;

    // direction of the differential plate
    let plate_x = theta.sin() * phi.cos();
    let plate_y = theta.sin() * phi.sin();
    let plate_z = theta.cos();

    let mut vf: f64 = 0.0;
    for alpha in &vec_alpha {
        for beta in &vec_beta {
            // coordinates of each mesh surface
            let mesh_x = xc + a * alpha.sin() * beta.cos();
            let mesh_y = yc + b * alpha.sin() * beta.sin();
            let mesh_z = zc + c * alpha.cos();
            // normal direction of the mesh surface
            let dir_x = b * c * alpha.sin() * beta.cos();
            let dir_y = a * c * alpha.sin() * beta.sin();
            let dir_z = a * b * alpha.cos();

            let dir_length = ((dir_x).powi(2) + (dir_y).powi(2) + (dir_z).powi(2)).sqrt();
            let mesh_size = alpha.sin() * dir_length * d_alpha * d_beta;
            // normalized mesh surface direction
            let norm_x = dir_x / dir_length;
            let norm_y = dir_y / dir_length;
            let norm_z = dir_z / dir_length;

            let mesh_s2 = mesh_x.powi(2) + mesh_y.powi(2) + mesh_z.powi(2);
            let mesh_s1 = mesh_s2.sqrt();
            let plate_to_mesh_x = mesh_x / mesh_s1;
            let plate_to_mesh_y = mesh_y / mesh_s1;
            let plate_to_mesh_z = mesh_z / mesh_s1;

            let cos_theta_1 =
                plate_to_mesh_x * plate_x + plate_to_mesh_y * plate_y + plate_to_mesh_z * plate_z;
            let cos_theta_2 =
                -plate_to_mesh_x * norm_x - plate_to_mesh_y * norm_y - plate_to_mesh_z * norm_z;

            if cos_theta_1 > 0.0 && cos_theta_2 > 0.0 {
                vf += mesh_size * cos_theta_1 * cos_theta_2 / (PI * mesh_s2);
            }
        }
    }

    vf
}

/// Numerical integration of the spheroid view factor from a plate element.
///
/// Spheroid: x²/a² + y²/a² + z²/b² = 1
/// Plate element position: (x0, y0, z0)
/// Plate element normal: (cos(theta)cos(phi), cos(theta)sin(phi), sin(theta))
///
/// # Arguments
///   * `a` - spheroid semi-major axis
///   * `b` - spheroid semi-minor axis
///   * `x0` - plate element x position
///   * `y0` - plate element y position
///   * `z0` - plate element z position
///   * `theta` - plate element elevation angle
///   * `phi` - plate element azimuthal angle
///   * `alpha_num` - number of divisions in the elevation direction
///   * `beta_num` - number of divisions in the azimuth direction
///
/// # Returns
///   * `Result<f64, ViewFactorError>` - view factor from the plate element
#[allow(clippy::too_many_arguments)]
pub fn diff_element_to_spheroid_numerical(
    a: f64,
    b: f64,
    x0: f64,
    y0: f64,
    z0: f64,
    theta: f64,
    phi: f64,
    alpha_num: usize,
    beta_num: usize,
) -> f64 {
    let d_alpha = PI / alpha_num as f64;
    let d_beta = 2.0 * PI / beta_num as f64;
    let vec_alpha = linspace_middle(-FRAC_PI_2, FRAC_PI_2, alpha_num);
    let vec_beta = linspace_middle(0.0, 2.0 * PI, beta_num);

    // create a mesh grid (alpha, beta)
    let mesh = vec_alpha
        .iter()
        .flat_map(|alpha| vec_beta.iter().map(move |beta| (alpha, beta)));

    // plate element normal vector
    let plate_dir_x = theta.cos() * phi.cos();
    let plate_dir_y = theta.cos() * phi.sin();
    let plate_dir_z = theta.sin();

    let vf: f64 = mesh
        .map(|(alpha, beta)| {
            // coordinates of each mesh surface
            let mesh_x = a * alpha.cos() * beta.cos();
            let mesh_y = a * alpha.cos() * beta.sin();
            let mesh_z = b * alpha.sin();

            // area of the mesh surface
            let coeff =
                f64::sqrt(a.powi(2) * alpha.sin().powi(2) + b.powi(2) * alpha.cos().powi(2));
            let mesh_area = a * alpha.cos() * coeff * d_alpha * d_beta;
            // normal direction of the mesh surface
            let dir_x = b * alpha.cos() * beta.cos() / coeff;
            let dir_y = b * alpha.cos() * beta.sin() / coeff;
            let dir_z = a * alpha.sin() / coeff;
            // distance from the plate element to the mesh surface
            let mesh_s2 = (mesh_x - x0).powi(2) + (mesh_y - y0).powi(2) + (mesh_z - z0).powi(2);
            let mesh_s1 = mesh_s2.sqrt();

            let cos_theta_1 = ((mesh_x - x0) * plate_dir_x
                + (mesh_y - y0) * plate_dir_y
                + (mesh_z - z0) * plate_dir_z)
                / mesh_s1;
            let cos_theta_2 =
                ((x0 - mesh_x) * dir_x + (y0 - mesh_y) * dir_y + (z0 - mesh_z) * dir_z) / mesh_s1;
            let df = mesh_area * cos_theta_1 * cos_theta_2 / (PI * mesh_s2);
            if cos_theta_1 > 0.0 && cos_theta_2 > 0.0 {
                df
            } else {
                0.0
            }
        })
        .sum();
    vf
}
