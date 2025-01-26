mod utils;

// cargo test --test diff_element_to_sphere_test
#[cfg(test)]
mod tests {
    use crate::utils::diff_element_to_ellipsoid_numerical;
    use std::f64::consts::PI;
    use sterad_view_factor::diff_element_to_sphere;

    #[test]
    fn test_diff_element_to_sphere() {
        let r = 1.0;
        let theta = PI / 4.0;
        let phi = 0.0;
        let xc = 0.0;
        let yc = 0.0;
        let zc = 2.0;
        let alpha_num = 10000;
        let beta_num = 10000;
        let vf_ana = diff_element_to_sphere::tilted(zc, r, theta).unwrap();
        let vf_num = diff_element_to_ellipsoid_numerical(
            r, r, r, theta, phi, xc, yc, zc, alpha_num, beta_num,
        );
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
