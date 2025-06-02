use std::f64::consts::{FRAC_PI_2, PI};
use sterad_view_factor::diff_element_to_spheroid;
mod utils;
use utils::diff_element_to_spheroid_numerical;

#[derive(Debug)]
struct TestResult {
    theta: f64,
    phi: f64,
    vf_ana: f64,
    vf_num: f64,
    error: f64,
}

#[allow(clippy::too_many_arguments)]
fn test_diff_element_to_spheroid(
    a: f64,
    b: f64,
    theta_vec: &[f64],
    phi_vec: &[f64],
    x0: f64,
    y0: f64,
    z0: f64,
    a_num: usize,
    b_num: usize,
) {
    // create a direct product (theta, phi)
    let theta_phi = theta_vec
        .iter()
        .flat_map(|theta| phi_vec.iter().map(move |phi| (theta, phi)));
    let results: Vec<TestResult> = theta_phi
        .map(|(theta, phi)| {
            let vf_ana =
                diff_element_to_spheroid::tilted_offset(a, b, x0, y0, z0, *theta, *phi).unwrap();
            let vf_num =
                diff_element_to_spheroid_numerical(a, b, x0, y0, z0, *theta, *phi, a_num, b_num);
            let error = (vf_ana - vf_num).abs();

            TestResult {
                theta: *theta,
                phi: *phi,
                vf_ana,
                vf_num,
                error,
            }
        })
        .collect();
    let failed_cases: Vec<&TestResult> = results
        .iter()
        .filter(|result| result.error > 1e-4)
        .collect();
    if !failed_cases.is_empty() {
        let num_failed = failed_cases.len();
        for result in failed_cases {
            println!(
                "theta: {:.1}, phi: {:.1}, vf_ana: {:.8}, vf_num: {:.8}, error: {:.8}",
                result.theta * 180.0 / PI,
                result.phi * 180.0 / PI,
                result.vf_ana,
                result.vf_num,
                result.error
            );
        }
        panic!(
            "Failed {}/{} test cases with error > 1e-4",
            num_failed,
            results.len()
        );
    }
}

#[test]
fn test_diff_element_to_spheroid_a05_r11() {
    let a = 0.5;
    let b = 1.0 / a;
    let r = 1.1;
    let theta_vec: Vec<f64> = (-6..7).map(|i| i as f64 * FRAC_PI_2 / 6.0).collect();
    let phi_vec: Vec<f64> = (-9..10).map(|i| i as f64 * PI / 9.0).collect();
    let psi_vec: Vec<f64> = (-6..7).map(|i| i as f64 * FRAC_PI_2 / 6.0).collect();
    let a_num = 2000;
    let b_num = 2000;

    psi_vec.iter().for_each(|psi| {
        let x0 = r * a * psi.cos();
        let y0 = 1.0;
        let z0 = r * b * psi.sin();
        println!(
            "Testing with (a, r, psi): {:.1}, {:.1}, {:.1}",
            a,
            r,
            psi * 180.0 / PI
        );
        test_diff_element_to_spheroid(a, b, &theta_vec, &phi_vec, x0, y0, z0, a_num, b_num);
    });
}
