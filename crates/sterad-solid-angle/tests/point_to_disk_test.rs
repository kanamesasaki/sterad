mod utils;

// cargo test --test point_to_disk_test
#[cfg(test)]
mod tests {
    use crate::utils::point_to_ellipse_numerical;
    use sterad_solid_angle::point_to_disk;

    #[test]
    fn test_point_to_disk_center() {
        let r = 1.0;
        let xc = 0.0;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let sa_ana = point_to_disk::center(zc, r).unwrap();
        let sa_num = point_to_ellipse_numerical(r, r, xc, yc, zc, a_num, b_num);
        assert!(
            (sa_ana - sa_num).abs() < 1e-6,
            "sa_ana: {}, sa_num: {}",
            sa_ana,
            sa_num
        );
    }

    #[test]
    fn test_point_to_disk_offset_1() {
        let rm = 1.0;
        let ro = 0.5;
        let xc = ro;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let sa_ana = point_to_disk::offset(zc, rm, ro).unwrap();
        let sa_num = point_to_ellipse_numerical(rm, rm, xc, yc, zc, a_num, b_num);
        assert!(
            (sa_ana - sa_num).abs() < 1e-6,
            "sa_ana: {}, sa_num: {}",
            sa_ana,
            sa_num
        );
    }

    #[test]
    fn test_point_to_disk_offset_2() {
        let rm = 1.0;
        let ro = 1.0;
        let xc = ro;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let sa_ana = point_to_disk::offset(zc, rm, ro).unwrap();
        let sa_num = point_to_ellipse_numerical(rm, rm, xc, yc, zc, a_num, b_num);
        assert!(
            (sa_ana - sa_num).abs() < 1e-6,
            "sa_ana: {}, sa_num: {}",
            sa_ana,
            sa_num
        );
    }

    #[test]
    fn test_point_to_disk_offset_3() {
        let rm = 1.0;
        let ro = 2.0;
        let xc = ro;
        let yc = 0.0;
        let zc = 1.0;
        let a_num = 10000;
        let b_num = 10000;
        let sa_ana = point_to_disk::offset(zc, rm, ro).unwrap();
        let sa_num = point_to_ellipse_numerical(rm, rm, xc, yc, zc, a_num, b_num);
        assert!(
            (sa_ana - sa_num).abs() < 1e-6,
            "sa_ana: {}, sa_num: {}",
            sa_ana,
            sa_num
        );
    }
}
