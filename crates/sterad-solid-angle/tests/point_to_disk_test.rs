mod utils;

// cargo test --test point_to_disk_test
#[cfg(test)]
mod tests {
    use crate::utils::point_to_ellipse_numerical;
    use sterad_solid_angle::point_to_disk;

    #[test]
    fn test_point_to_disk() {
        let a = 1.0;
        let b = 1.0;
        let xc = 0.0;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let sa_ana = point_to_disk::center(zc, a).unwrap();
        let sa_num = point_to_ellipse_numerical(a, b, xc, yc, zc, a_num, b_num);
        assert!(
            (sa_ana - sa_num).abs() < 1e-6,
            "sa_ana: {}, sa_num: {}",
            sa_ana,
            sa_num
        );
    }
}
