mod utils;

// cargo test --test diff_element_to_ellipsoi_test
#[cfg(test)]
mod tests {
    use crate::utils::diff_element_to_ellipsoid_numerical;
    use std::f64::consts::PI;
    use sterad_view_factor::diff_element_to_ellipsoid;

    #[test]
    fn test_diff_element_to_ellipsoid() {
        let a = 3.0;
        let b = 2.0;
        let c = 1.0;
        let theta = PI / 4.0;
        let phi = 0.0;
        let xc = 1.0;
        let yc = 1.0;
        let zc = 4.0;
        let alpha_num = 10000;
        let beta_num = 10000;
        let vf_ana =
            diff_element_to_ellipsoid::tilted_offset(a, b, c, xc, yc, zc, theta, phi).unwrap();
        let vf_num = diff_element_to_ellipsoid_numerical(
            a, b, c, theta, phi, xc, yc, zc, alpha_num, beta_num,
        );
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
