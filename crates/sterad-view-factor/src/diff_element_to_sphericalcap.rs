use crate::diff_element_to_sphere;
use crate::error::ViewFactorError;
use crate::vecmath::Vector3f;
use std::f64::consts::{FRAC_PI_2, PI};

// Case Description

// Case 1 : Cap Orientation Mismatch
// Case 2 : Infinitesimal Surface Direction Mismatch
// Case 3 : Full Sphere
// Case 4 : Partial Sphere
// Case 5 : Full Spherical Cap
// Case 6 : Partial Cap, No Intersection 1
// Case 7 : Partial Cap, No Intersection 2
// Case 8 : Partial Cap, No Intersection 3
// Case 9 : Partial Cap, No Intersection 4
// Case 10: Partial Cap, One Intersection 1
// Case 11: Partial Cap, One Intersection 2
// Case 12: Partial Cap, Two Intersection 1
// Case 13: Partial Cap, Two Intersection 2
// Case 14: Partial Cap, Two Intersection 3
// Case 15: Partial Cap, Two Intersection 4
// Case 16: Small Cap, View Direction Mismatch
// Case 17: Small Cap, No Intersection
// Case 18: Small Cap, No Intersection phi=0
// Case 19: Small Cap, Two Intersections phi=0
// Case 20: Small Cap, Two Intersections 1
// Case 21: Small Cap, Two Intersections 2

fn l1(beta0: f64, beta1: f64, omega: f64, d: f64, rs: f64, gamma: f64) -> f64 {
    // Initialize n1 vector
    let n1 = Vector3f {
        x: omega.sin() * gamma.cos(),
        y: omega.sin() * gamma.sin(),
        z: omega.cos(),
    };

    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let r = rs * cos_theta;
    let h = d * cos_theta.powi(2);

    // Calculate denominator once
    let denom = 2.0 * PI * (r.powi(2) + h.powi(2));

    // Calculate l1_integration vector
    let l1_integration = Vector3f {
        x: h * r * (beta1.sin() - beta0.sin()) / denom,
        y: h * r * (-beta1.cos() + beta0.cos()) / denom,
        z: r.powi(2) * (-beta1 + beta0) / denom,
    };

    n1.dot(&l1_integration)
}

fn l2(
    alpha0: f64,
    alpha1: f64,
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<f64, ViewFactorError> {
    // Input validation
    if !(-PI..=PI).contains(&alpha0) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "alpha0",
            message: "alpha0 has to be between -pi and pi".to_string(),
        });
    }
    if !(-PI..=PI).contains(&alpha1) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "alpha1",
            message: "alpha1 has to be between -pi and pi".to_string(),
        });
    }

    // Initialize n1 vector
    let n1 = Vector3f {
        x: omega.sin() * gamma.cos(),
        y: omega.sin() * gamma.sin(),
        z: omega.cos(),
    };

    // Intermediate calculations
    let sin_theta = rs / d;
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    // Calculate A and related values
    let a = 2.0 * sin_theta * sin_psi * sin_phi
        / (1.0 + sin_theta.powi(2) - 2.0 * sin_theta * cos_psi * cos_phi);
    let l2 = (1.0 - sin_theta * cos_psi * cos_phi) * sin_theta * sin_psi;
    let l3 = sin_theta.powi(2) * sin_psi.powi(2) * sin_phi;
    let m2 = sin_theta * sin_psi * cos_phi - sin_theta.powi(2) * cos_psi * sin_psi;
    let n2 = sin_theta.powi(2) * cos_psi * sin_psi * sin_phi;
    let n3 = sin_theta.powi(2) * sin_psi.powi(2) * cos_phi;

    let denom = 4.0 * std::f64::consts::PI * sin_theta * sin_psi * sin_phi;
    // Calculate l2_integration components
    let l2_integration = Vector3f {
        x: a / denom
            * (-l2 / a * (alpha1 - alpha0)
                + (l2 / a - l3) * 2.0 / (1.0 - a.powi(2)).sqrt()
                    * (l2_arctan(a, alpha1) - l2_arctan(a, alpha0))),
        y: 1.0 / denom * m2 * ((1.0 - a * (alpha1.cos())).ln() - (1.0 - a * (alpha0.cos())).ln()),
        z: a / denom
            * (-n2 / a * (alpha1 - alpha0)
                + (n2 / a - n3) * 2.0 / (1.0 - a.powi(2)).sqrt()
                    * (l2_arctan(a, alpha1) - l2_arctan(a, alpha0))),
    };

    Ok(n1.dot(&l2_integration))
}

fn l2_arctan(a: f64, alpha: f64) -> f64 {
    if (alpha - PI).abs() < f64::EPSILON {
        FRAC_PI_2
    } else if (alpha + PI).abs() < f64::EPSILON {
        -FRAC_PI_2
    } else {
        let sqrt_term = ((1.0 + a) / (1.0 - a)).sqrt();
        let tan_term = (alpha / 2.0).tan();
        (sqrt_term * tan_term).atan()
    }
}

pub fn l4(x0: &Vector3f, x1: &Vector3f, omega: f64, gamma: f64) -> Result<f64, ViewFactorError> {
    // Input validation
    if (x0.z - x1.z).abs() > f64::EPSILON {
        return Err(ViewFactorError::InvalidInput {
            param_name: "x0.z, x1.z",
            message: "x0.z and x1.z have to be equal".to_string(),
        });
    }
    let h = x0.z;

    // Initialize n1 vector
    let n1 = Vector3f {
        x: omega.sin() * gamma.cos(),
        y: omega.sin() * gamma.sin(),
        z: omega.cos(),
    };

    // Calculate intermediate values
    let dx = Vector3f {
        x: x1.x - x0.x,
        y: x1.y - x0.y,
        z: x1.z - x0.z,
    };

    let a = dx.x.powi(2) + dx.y.powi(2);
    let b = x0.x * dx.x + x0.y * dx.y;
    let c = x0.x.powi(2) + x0.y.powi(2) + h.powi(2);
    let d = (a * c - b.powi(2)).sqrt();
    let e = ((a + b) / d).atan() - (b / d).atan();

    // Calculate l4_integration
    let l4_integration = Vector3f {
        x: h * dx.y / (2.0 * PI * d) * e,
        y: -h * dx.x / (2.0 * PI * d) * e,
        z: (x0.y * dx.x - x0.x * dx.y) / (2.0 * PI * d) * e,
    };

    Ok(n1.dot(&l4_integration))
}

fn f3(omega: f64, d: f64, rs: f64, gamma: f64) -> f64 {
    // return F3 view factor
    l1(PI, -PI, omega, d, rs, gamma)
}

fn f4(omega: f64, d: f64, rs: f64, gamma: f64) -> Result<f64, ViewFactorError> {
    // return F4 view factor
    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let r = rs * cos_theta;
    let h = d * cos_theta.powi(2);

    let cos_beta0 = -h / r * (PI - omega).tan();
    let beta0 = cos_beta0.acos();

    let x0 = Vector3f {
        x: r * (gamma - beta0).cos(),
        y: r * (gamma - beta0).sin(),
        z: h,
    };
    let x1 = Vector3f {
        x: r * (gamma + beta0).cos(),
        y: r * (gamma + beta0).sin(),
        z: h,
    };

    let vf = l1(gamma + beta0, gamma - beta0, omega, d, rs, gamma) + l4(&x0, &x1, omega, gamma)?;
    Ok(vf)
}

fn f5(omega: f64, d: f64, rs: f64, phi: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
    // return F5 view factor
    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();

    let cos_alpha0 = (sin_theta - cos_psi * cos_phi) / (sin_psi * sin_phi);
    let alpha0 = cos_alpha0.acos();
    let cos_beta0 = -cos_psi * sin_phi / cos_theta
        + (sin_theta - cos_psi * cos_phi) / cos_theta * cos_phi / sin_phi;
    let beta0 = cos_beta0.acos();

    let vf = l1(2.0 * PI - beta0, beta0, omega, d, rs, gamma)
        + l2(alpha0, -alpha0, omega, d, rs, phi, gamma, psi)?;
    Ok(vf)
}

fn f7(omega: f64, d: f64, rs: f64, phi: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
    // return F7 view factor
    let vf5 = f5(omega, d, rs, phi, gamma, psi)?;
    let vf4 = f4(omega, d, rs, gamma)?;
    let vf3 = f3(omega, d, rs, gamma);
    Ok(vf5 + vf4 - vf3)
}

fn f10(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    z2: f64,
) -> Result<f64, ViewFactorError> {
    // return F10 view factor
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let cos_gamma = gamma.cos();
    let sin_gamma = gamma.sin();
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();

    let h = d * (1.0 - sin_theta.powi(2));
    let r = d * sin_theta * cos_theta;

    let cos_alpha1 = (sin_theta - cos_psi * cos_phi) / (sin_psi * sin_phi);
    let alpha1 = cos_alpha1.acos();

    let x2 = (d - z2) * cos_phi / sin_phi - cos_psi / sin_phi * rs;
    let y2 = -cos_gamma / sin_gamma * x2 - cos_omega / (sin_omega * sin_gamma) * z2;
    let cos_alpha2 = (x2 + rs * cos_psi * sin_phi) / (rs * sin_psi * cos_phi);
    let sin_alpha2 = y2 / (rs * sin_psi);
    let alpha2 = if cos_phi.abs() < f64::EPSILON {
        sin_alpha2.asin()
    } else {
        sin_alpha2.atan2(cos_alpha2)
    };

    let cos_beta1 = -cos_psi * sin_phi / cos_theta
        + (sin_theta - cos_psi * cos_phi) / cos_theta * cos_phi / sin_phi;
    let beta1 = cos_beta1.acos();

    let cos_beta0 = -h / r * (PI / 2.0 - omega).tan();
    let beta0 = cos_beta0.acos();
    let x_start = Vector3f {
        x: x2 * h / z2,
        y: y2 * h / z2,
        z: z2 * h / z2,
    };
    let x_end = Vector3f {
        x: r * (gamma + beta0).cos(),
        y: r * (gamma + beta0).sin(),
        z: h,
    };

    let vf = l1(gamma + beta0, beta1, omega, d, rs, gamma)
        + l2(alpha1, alpha2, omega, d, rs, phi, gamma, psi)?
        + l4(&x_start, &x_end, omega, gamma)?;
    Ok(vf)
}

fn f12(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    // return F12 view factor
    let sin_theta = rs / d;
    let h = d * (1.0 - sin_theta.powi(2));
    let x_start = Vector3f {
        x: x2.x * h / x2.z,
        y: x2.y * h / x2.z,
        z: x2.z * h / x2.z,
    };
    let x_end = Vector3f {
        x: x1.x * h / x1.z,
        y: x1.y * h / x1.z,
        z: x1.z * h / x1.z,
    };

    let vf =
        l2(alpha1, alpha2, omega, d, rs, phi, gamma, psi)? + l4(&x_start, &x_end, omega, gamma)?;
    Ok(vf)
}

fn f13(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    let vf5 = f5(omega, d, rs, phi, gamma, psi)?;
    let vf12 = f12(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?;
    Ok(vf5 - vf12)
}

fn f14(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    let vf3 = f3(omega, d, rs, gamma);
    let vf4 = f4(omega, d, rs, gamma)?;
    let vf5 = f5(omega, d, rs, phi, gamma, psi)?;
    let vf12 = f12(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?;
    Ok(vf4 - (vf3 - vf5 - vf12))
}

fn f15(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    let vf4 = f4(omega, d, rs, gamma)?;
    let vf12 = f12(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?;
    Ok(vf4 - vf12)
}

fn f17(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<f64, ViewFactorError> {
    let vf = l2(PI, -PI, omega, d, rs, phi, gamma, psi)?;
    Ok(vf)
}

fn f18(omega: f64, d: f64, rs: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    let theta_ = (rs * sin_psi / (d - rs * cos_psi)).atan();
    let rs_transformed = rs * sin_psi / theta_.cos();
    let d_transformed = d - rs * cos_psi + rs_transformed * theta_.sin();
    let vf = f3(omega, d_transformed, rs_transformed, gamma);
    Ok(vf)
}

fn f19(omega: f64, d: f64, rs: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    let theta_ = (rs * sin_psi / (d - rs * cos_psi)).atan();
    let rs_transformed = rs * sin_psi / theta_.cos();
    let d_transformed = d - rs * cos_psi + rs_transformed * theta_.sin();
    let vf = f4(omega, d_transformed, rs_transformed, gamma)?;
    Ok(vf)
}

fn f20(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    let vf = l2(alpha1, alpha2, omega, d, rs, phi, gamma, psi)? + l4(&x1, &x2, omega, gamma)?;
    Ok(vf)
}

fn f21(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
    alpha1: f64,
    alpha2: f64,
) -> Result<f64, ViewFactorError> {
    let vf = l2(PI, -PI, omega, d, rs, phi, gamma, psi)?
        + l2(alpha1, alpha2, omega, d, rs, phi, gamma, psi)?
        + l4(&x1, &x2, omega, gamma)?;
    Ok(vf)
}
