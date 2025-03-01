use std::f64::consts::PI;
use sterad_view_factor::error::ViewFactorError;

#[allow(clippy::too_many_arguments)]
fn diff_element_to_sphericalcap_numerical(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    eta_steps: usize,
    lambda_steps: usize,
) -> Result<f64, ViewFactorError> {
    if !(0.0..=PI).contains(&omega) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "omega",
            message: "omega has to be 0 <= omega <= pi".to_string(),
        });
    }
    if rs < 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "rs",
            message: "rs has to be non-negative".to_string(),
        });
    }
    if d < rs {
        return Err(ViewFactorError::InvalidInput {
            param_name: "d",
            message: "d has to be greater than rs".to_string(),
        });
    }
    if !(0.0..=PI).contains(&phi) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "phi",
            message: "phi has to be 0 <= phi <= pi".to_string(),
        });
    }
    if psi <= 0.0 || psi > PI {
        return Err(ViewFactorError::InvalidInput {
            param_name: "psi",
            message: "psi has to be 0 < psi <= pi".to_string(),
        });
    }

    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let theta = cos_theta.acos();

    let deta = (PI / 2.0 - theta) / (eta_steps as f64);
    let dlambda = 2.0 * PI / (lambda_steps as f64);

    let mut vf = 0.0;
    for j in 0..eta_steps {
        let etaj = PI / 2.0 + theta + (j as f64 + 0.5) * deta;
        for k in 0..lambda_steps {
            let lambdak = -PI + (k as f64 + 0.5) * dlambda;

            let x = rs * etaj.sin() * lambdak.cos();
            let y = rs * etaj.sin() * lambdak.sin();
            let z = rs * etaj.cos() + d;

            let view =
                x * omega.sin() * gamma.cos() + y * omega.sin() * gamma.sin() + z * omega.cos();
            let cap = -(x + rs * psi.cos() * phi.sin()) * phi.sin()
                - (z - d + rs * psi.cos() * phi.cos()) * phi.cos();

            if view > 0.0 && cap > 0.0 {
                let s = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
                let cos_theta0 = (x * omega.sin() * gamma.cos()
                    + y * omega.sin() * gamma.sin()
                    + z * omega.cos())
                    / s;
                let cos_theta1 = (-x.powi(2) - y.powi(2) - z * (z - d)) / (s * rs);
                vf += cos_theta0 * cos_theta1 * rs.powi(2) * etaj.sin() * deta * dlambda
                    / (s.powi(2) * PI);
            }
        }
    }
    Ok(vf)
}

mod tests {
    use std::f64::consts::PI;
    use sterad_view_factor::diff_element_to_sphericalcap;

    #[test]
    fn test_diff_element_to_sphericalcap_case1() {
        let omega = 0.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 2.0 + f64::EPSILON;
        let gamma = 0.0;
        let psi = PI / 6.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = 0.0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.0, vf_ref);
        assert_eq!(vf_ana.1, 1);
        assert_eq!(vf_ana.0, vf_num);
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case2() {
        let omega = 2.0 / 3.0 * PI + f64::EPSILON;
        let d = 2.0;
        let rs = 1.0;
        let phi = 0.0;
        let gamma = 0.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = 0.0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.0, vf_ref);
        assert_eq!(vf_ana.1, 2);
        assert_eq!(vf_ana.0, vf_num);
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case3() {
        let omega = 0.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = 0.0;
        let gamma = 0.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = (rs / d).powi(2);
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 3);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-7,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case4() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = 0.0;
        let gamma = 0.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let r = f64::sqrt(3.0) / 2.0;
        let h = 1.5;
        let vf_ref = -r * h / (PI * (r.powi(2) + h.powi(2))) + (r / h).atan() / PI;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 4);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-7,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case5() {
        let omega = 0.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 2.0;
        let gamma = 0.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = (rs / d).powi(2) / 2.0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 5);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-7,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case6() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 4.0 + f64::EPSILON;
        let gamma = 0.0;
        let psi = PI / 4.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = 0.0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.0, vf_ref);
        assert_eq!(vf_ana.1, 6);
        assert_eq!(vf_ana.0, vf_num);
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case7() {
        let omega = PI * 3.0 / 8.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 3.0;
        let gamma = 0.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 3000, 3000,
        )
        .unwrap();
        assert_eq!(vf_ana.1, 7);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case10() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 3.0;
        let gamma = PI / 3.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.1, 10);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case11() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 3.0;
        let gamma = -PI / 3.0;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.1, 11);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case8() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 3.0;
        let gamma = PI;
        let psi = PI / 2.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let r = f64::sqrt(3.0) / 2.0;
        let h = 1.5;
        let vf_ref = -r * h / (PI * (r.powi(2) + h.powi(2))) + (r / h).atan() / PI;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 8);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-7,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case9() {
        let omega = PI * 3.0 / 8.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 3.0;
        let gamma = PI;
        let psi = PI / 6.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 3000, 3000,
        )
        .unwrap();
        assert_eq!(vf_ana.1, 9);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case16() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 10.0 + f64::EPSILON;
        let gamma = 0.0;
        let psi = PI / 10.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let vf_ref = 0.0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert_eq!(vf_ana.0, vf_ref);
        assert_eq!(vf_ana.1, 16);
        assert_eq!(vf_ana.0, vf_num);
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case18() {
        let omega = 0.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = 0.0;
        let gamma = 0.0;
        let psi = PI / 6.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let r: f64 = 0.5;
        let h: f64 = 2.0 - f64::sqrt(3.0) / 2.0;
        let vf_ref = r.powi(2) / (r.powi(2) + h.powi(2));
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 18);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-7,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }

    #[test]
    fn test_diff_element_to_sphericalcap_case19() {
        let omega = PI / 2.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = 0.0;
        let gamma = 0.0;
        let psi = PI / 6.0;
        let vf_ana =
            diff_element_to_sphericalcap::sphericalcap(omega, d, rs, phi, gamma, psi).unwrap();
        let r: f64 = 0.5;
        let h: f64 = 2.0 - f64::sqrt(3.0) / 2.0;
        let vf_ref = -r * h / (PI * (r.powi(2) + h.powi(2))) + (r / h).atan() / PI;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, 2000, 2000,
        )
        .unwrap();
        assert!(
            (vf_ana.0 - vf_ref).abs() < 1e-10,
            "vf_ana: {}, vf_ref: {}",
            vf_ana.0,
            vf_ref
        );
        assert_eq!(vf_ana.1, 19);
        assert!(
            (vf_ana.0 - vf_num).abs() < 1e-8,
            "vf_ana: {}, vf_num: {}",
            vf_ana.0,
            vf_num
        );
    }
}
