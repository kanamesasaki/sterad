mod utils;

// cargo test --test diff_element_to_disk_test
#[cfg(test)]
mod tests {
    use crate::utils::diff_element_to_ellipse_numerical;
    use std::f64::consts::PI;
    use sterad_view_factor::diff_element_to_disk;

    #[test]
    fn test_diff_element_to_ellipse() {
        let a = 1.0;
        let b = 1.0;
        let theta = PI / 4.0;
        let phi = 0.0;
        let xc = 0.0;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let vf_ana = diff_element_to_disk::tilted_center(zc, a, theta).unwrap();
        let vf_num = diff_element_to_ellipse_numerical(a, b, theta, phi, xc, yc, zc, a_num, b_num);
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
