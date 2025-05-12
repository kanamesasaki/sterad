mod utils;

// cargo test --test diff_element_to_ellipse_test
#[cfg(test)]
mod tests {
    use crate::utils::diff_element_to_ellipse_numerical;
    use std::f64::consts::{FRAC_PI_2, PI};
    use sterad_view_factor::diff_element_to_ellipse;

    #[derive(Debug)]
    struct TestResult {
        theta: f64,
        phi: f64,
        vf_ana: f64,
        vf_num: f64,
        error: f64,
    }

    #[allow(clippy::too_many_arguments)]
    fn test_diff_element_to_ellipse_center(
        a: f64,
        b: f64,
        theta_vec: &[f64],
        phi_vec: &[f64],
        xc: f64,
        yc: f64,
        zc: f64,
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
                    diff_element_to_ellipse::tilted_center(zc, a, b, *theta, *phi).unwrap();
                let vf_num =
                    diff_element_to_ellipse_numerical(a, b, *theta, *phi, xc, yc, zc, a_num, b_num);
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
            .filter(|result| result.error > 1e-6)
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
                "Failed {}/{} test cases with error > 1e-6",
                num_failed,
                results.len()
            );
        }
    }

    /// Test case for ellipse with a = 2.0, b = 0.5
    /// Plate element at the central axis with zc = 1.0
    /// Last verified: 2025-05-12
    #[test]
    fn test_diff_element_to_ellipse_center_a20_h010() {
        let a = 2.0;
        let b = 1.0 / a;
        let theta_vec: Vec<f64> = (0..10).map(|i| i as f64 * PI / 9.0).collect();
        let phi_vec: Vec<f64> = (0..7).map(|i| i as f64 * FRAC_PI_2 / 6.0).collect();
        let xc = 0.0;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;

        test_diff_element_to_ellipse_center(a, b, &theta_vec, &phi_vec, xc, yc, zc, a_num, b_num);
    }

    /// Test case for ellipse with a = 2.0, b = 0.5
    /// Plate element at the central axis with zc = 0.1
    /// Last verified: 2025-05-12
    #[test]
    fn test_diff_element_to_ellipse_center_a20_h001() {
        let a = 2.0;
        let b = 1.0 / a;
        let theta_vec: Vec<f64> = (0..10).map(|i| i as f64 * PI / 9.0).collect();
        let phi_vec: Vec<f64> = (0..7).map(|i| i as f64 * FRAC_PI_2 / 6.0).collect();
        let xc = 0.0;
        let yc = 0.0;
        let zc = 0.1;
        let a_num = 10000;
        let b_num = 10000;

        test_diff_element_to_ellipse_center(a, b, &theta_vec, &phi_vec, xc, yc, zc, a_num, b_num);
    }

    /// Test case for ellipse with a = 2.0, b = 0.5
    /// Plate element at the central axis with zc = 10.0
    /// Last verified: 2025-05-12
    #[test]
    fn test_diff_element_to_ellipse_center_a20_h100() {
        let a = 2.0;
        let b = 1.0 / a;
        let theta_vec: Vec<f64> = (0..10).map(|i| i as f64 * PI / 9.0).collect();
        let phi_vec: Vec<f64> = (0..7).map(|i| i as f64 * FRAC_PI_2 / 6.0).collect();
        let xc = 0.0;
        let yc = 0.0;
        let zc = 10.0;
        let a_num = 10000;
        let b_num = 10000;

        test_diff_element_to_ellipse_center(a, b, &theta_vec, &phi_vec, xc, yc, zc, a_num, b_num);
    }

    #[test]
    fn test_diff_element_to_ellipse_offset() {
        let a = 2.0;
        let b = 1.0;
        let theta = PI / 6.0;
        let phi = PI / 6.0;
        let xc = 1.0;
        let yc = 1.0;
        let zc = 2.0;
        let a_num = 10000;
        let b_num = 10000;
        let vf_ana = diff_element_to_ellipse::tilted_offset(a, b, xc, yc, zc, theta, phi).unwrap();
        let vf_num = diff_element_to_ellipse_numerical(a, b, theta, phi, xc, yc, zc, a_num, b_num);
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
