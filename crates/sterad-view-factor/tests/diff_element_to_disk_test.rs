mod utils;

// cargo test --test diff_element_to_disk_test
#[cfg(test)]
mod tests {
    use crate::utils::diff_element_to_ellipse_numerical;
    use std::f64::consts::PI;
    use sterad_view_factor::diff_element_to_disk;

    #[test]
    fn test_diff_element_to_ellipse() {
        let r = 2.0;
        let theta = PI / 2.0;
        let phi = 0.0;
        let xc = 0.0;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let vf_ana = diff_element_to_disk::tilted_center(zc, r, theta).unwrap();
        let vf_ref = -r * zc / (PI * (r.powi(2) + zc.powi(2))) + (r / zc).atan() / PI;
        let vf_num = diff_element_to_ellipse_numerical(r, r, theta, phi, xc, yc, zc, a_num, b_num);
        assert!(
            (vf_ana - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana,
            vf_ref
        );
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
