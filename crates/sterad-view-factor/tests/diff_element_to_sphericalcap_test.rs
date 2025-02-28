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
    fn test_diff_element_to_sphericalcap() {
        let omega = PI / 4.0;
        let d = 2.0;
        let rs = 1.0;
        let phi = PI / 4.0;
        let gamma = 0.0;
        let psi = PI / 4.0;
        let eta_num = 10000;
        let lambda_num = 10000;
        let vf_ana = diff_element_to_sphericalcap::sphericalcap(rs, d, omega, phi, gamma, psi)
            .unwrap()
            .0;
        let vf_num = super::diff_element_to_sphericalcap_numerical(
            omega, d, rs, phi, gamma, psi, eta_num, lambda_num,
        )
        .unwrap();
        assert!(
            (vf_ana - vf_num).abs() < 1e-6,
            "vf_ana: {}, vf_num: {}",
            vf_ana,
            vf_num
        );
    }
}
