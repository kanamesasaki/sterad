use std::f64::consts::{FRAC_PI_2, PI};
use sterad_view_factor::diff_element_to_spheroid;
mod utils;
use utils::diff_element_to_spheroid_numerical;

#[test]
fn test_spheroid_plate_element_on_z() {
    // Test case: spheroid on x-y plane, differential element at origin
    let a = 2.0; // semi-major axis
    let b = 1.0; // semi-minor axis
    let x0 = 0.0; // element x position
    let y0 = 0.0; // element y position
    let z0 = 3.0; // element z position (above the spheroid)
    let theta = -FRAC_PI_2; // element points straight down
    let phi = 0.0;

    let vf_ana = diff_element_to_spheroid::tilted_offset(a, b, x0, y0, z0, theta, phi).unwrap();
    // let vf_ref = a.powi(2) / (a.powi(2) + z0.powi(2) - b.powi(2));
    let vf_num = diff_element_to_spheroid_numerical(a, b, x0, y0, z0, theta, phi, 1000, 1000);
    assert!(
        (vf_ana - vf_num).abs() < 1e-6,
        "vf_ana: {}, vf_num: {}",
        vf_ana,
        vf_num
    );
}
