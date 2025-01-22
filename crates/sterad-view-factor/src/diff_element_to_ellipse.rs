use crate::error::ViewFactorError;
use std::f64::consts::PI;

/// B-18: Planar element to elliptical plate in plane parallel to element. Normal to element passes through center of plate.
///
/// # Arguments
///
/// * `h` - A floating point number representing the distance between the planar element and the circular disk.
/// * `a` - A floating point number representing the one of the axes of the elliptical plate.
/// * `b` - A floating point number representing the one of the axes of the elliptical plate.
///
/// # Returns
///
/// The View Factor value.
///
/// # Errors
///
/// Returns an error if `h` or `a` or `b` is less than or equal to zero.
pub fn parallel_center(h: f64, a: f64, b: f64) -> Result<f64, ViewFactorError> {
    if h <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if a <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a",
            message: "r must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "b",
            message: "r must be greater than 0".to_string(),
        });
    }
    let a_over_h: f64 = a / h;
    let b_over_h: f64 = b / h;
    let f: f64 = a_over_h * b_over_h / ((1.0 + a_over_h.powi(2)) * (1.0 + b_over_h.powi(2))).sqrt();
    Ok(f)
}

pub fn tilted_center(h: f64, a: f64, b: f64, theta: f64, phi: f64) -> Result<f64, ViewFactorError> {
    if h <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "h",
            message: "h must be greater than 0".to_string(),
        });
    }
    if a <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a",
            message: "r must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "b",
            message: "r must be greater than 0".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }
    let a_over_h: f64 = a / h;
    let b_over_h: f64 = b / h;
    let f: f64 = a_over_h * b_over_h / ((1.0 + a_over_h.powi(2)) * (1.0 + b_over_h.powi(2))).sqrt();
    Ok(f)
}

fn tilted_center_full(h: f64, a: f64, b: f64, theta: f64) -> f64 {
    a * b * theta.cos() / ((h.powi(2) + a.powi(2)) * (h.powi(2) + b.powi(2))).sqrt()
}

fn tilted_center_partial(
    a: f64,
    b: f64,
    h: f64,
    theta: f64,
    phi: f64,
) -> Result<f64, ViewFactorError> {
    if !(0.0..=PI / 2.0).contains(&phi) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "phi",
            message: "phi must be between 0 and PI/2".to_string(),
        });
    }

    let (x0, y0, x1, y1) = if (phi - 0.0).abs() < 1e-10 {
        let x0: f64 = -h * theta.cos() / theta.sin();
        let y0: f64 = b / a * (a.powi(2) - x0.powi(2)).sqrt();
        (x0, y0, x0, -y0)
    } else if (phi - PI / 2.0).abs() < 1e-10 {
        let y0: f64 = -h * theta.cos() / theta.sin();
        let x0: f64 = -a / b * (b.powi(2) - y0.powi(2)).sqrt();
        (x0, y0, -x0, y0)
    } else {
        let capital_a: f64 = 1.0 / a.powi(2) + phi.cos().powi(2) / (b.powi(2) * phi.sin().powi(2));
        let capital_b: f64 =
            h * theta.cos() * phi.cos() / (b.powi(2) * theta.sin() * phi.sin().powi(2));
        let capital_c: f64 = h.powi(2) * theta.cos().powi(2)
            / (b.powi(2) * theta.sin().powi(2) * phi.sin().powi(2))
            - 1.0;
        let x0 = (-capital_b - (capital_b.powi(2) - capital_a * capital_c).sqrt()) / capital_a;
        let x1 = (-capital_b + (capital_b.powi(2) - capital_a * capital_c).sqrt()) / capital_a;
        let y0 = -x0 * phi.cos() / phi.sin() - h * theta.cos() / (phi.sin() * theta.sin());
        let y1 = -x1 * phi.cos() / phi.sin() - h * theta.cos() / (phi.sin() * theta.sin());
        (x0, y0, x1, y1)
    };

    let ang0 = (a * y0).atan2(b * x0);
    let ang1 = (a * y1).atan2(b * x1);

    let l = theta.sin() * phi.cos();
    let m = theta.sin() * phi.sin();
    let n = theta.cos();

    let vf: f64 = if y0 < 0.0 && y1 < 0.0 {
        let a2: f64 = (a.powi(2) - b.powi(2)) / (a.powi(2) + h.powi(2));
        let vf1: f64 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2)))
            * integral_cosx(PI, ang1, a2)
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(PI, ang1, a2)
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(PI, ang1, a2);
        let capital_k = (x0 - x1).powi(2) + (y0 - y1).powi(2);
        let capital_l = (x0 - x1) * x1 + (y0 - y1) * y1;
        let capital_m = x1.powi(2) + y1.powi(2) + h.powi(2);
        let capital_n = l * h * (y0 - y1) - m * h * (x0 - x1) + n * (y1 * x0 - x1 * y0);
        let vf2 = capital_n / (2.0 * PI * (-capital_l.powi(2) + capital_k * capital_m).sqrt())
            * ((capital_k + capital_l) / (-capital_l.powi(2) + capital_k * capital_m).sqrt())
                .atan()
            - (capital_l / (-capital_l.powi(2) + capital_k * capital_m).sqrt()).atan();

        let vf3 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_cosx(ang0, -PI, a2)
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(ang0, -PI, a2)
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(ang0, -PI, a2);

        vf1 + vf2 + vf3
    } else {
        let a2: f64 = (a.powi(2) - b.powi(2)) / (a.powi(2) + h.powi(2));
        let vf1: f64 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2)))
            * integral_cosx(ang0, ang1, a2)
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(ang0, ang1, a2)
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(ang0, ang1, a2);
        let capital_k = (x0 - x1).powi(2) + (y0 - y1).powi(2);
        let capital_l = (x0 - x1) * x1 + (y0 - y1) * y1;
        let capital_m = x1.powi(2) + y1.powi(2) + h.powi(2);
        let capital_n = l * h * (y0 - y1) - m * h * (x0 - x1) + n * (y1 * x0 - x1 * y0);
        let vf2 = capital_n / (2.0 * PI * (-capital_l.powi(2) + capital_k * capital_m).sqrt())
            * ((capital_k + capital_l) / (-capital_l.powi(2) + capital_k * capital_m).sqrt())
                .atan()
            - (capital_l / (-capital_l.powi(2) + capital_k * capital_m).sqrt()).atan();

        vf1 + vf2
    };

    Ok(vf)
}

/// Calculate the integral of 1 / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -pi <= x0 <= pi.
/// * `x1` - Upper limit of the integral, -pi <= x1 <= pi.
/// * `a2` - a^2, 0 <= a2 < 1.
///
/// # Returns
///
/// The integral value as f64.
///
fn integral_1(x0: f64, x1: f64, a2: f64) -> f64 {
    let val0: f64 = if x0 == -PI {
        -PI / (1.0 - a2).sqrt()
    } else if x0 == PI {
        PI / (1.0 - a2).sqrt()
    } else {
        1.0 / (1.0 - a2).sqrt() * ((1.0 - a2).sqrt() * x0.sin()).atan2(x0.cos())
    };

    let val1: f64 = if x1 == -PI {
        -PI / (1.0 - a2).sqrt()
    } else if x1 == PI {
        PI / (1.0 - a2).sqrt()
    } else {
        1.0 / (1.0 - a2).sqrt() * ((1.0 - a2).sqrt() * x1.sin()).atan2(x1.cos())
    };

    val1 - val0
}

/// Calculate the integral of cos(x) / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -pi <= x0 <= pi.
/// * `x1` - Upper limit of the integral, -pi <= x1 <= pi.
/// * `a2` - a^2, 0 <= a2 < 1.
///
/// # Returns
///
/// The integral value.
fn integral_cosx(x0: f64, x1: f64, a2: f64) -> f64 {
    let a: f64 = a2.sqrt();
    let val0: f64 = (a * x0.sin()).atanh() / a;
    let val1: f64 = (a * x1.sin()).atanh() / a;
    val1 - val0
}

/// Calculate the integral of sin(x) / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -pi <= x0 <= pi.
/// * `x1` - Upper limit of the integral, -pi <= x1 <= pi.
/// * `a2` - a^2, 0 <= a2 < 1.
///
/// # Returns
///
/// The integral value.
fn integral_sinx(x0: f64, x1: f64, a2: f64) -> f64 {
    let a: f64 = a2.sqrt();
    let val0: f64 = if x0 == -PI || x0 == PI {
        0.0
    } else {
        -1.0 / (a * (1.0 - a2).sqrt()) * ((a - (x0 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
            + ((a + (x0 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
    };

    let val1 = if x1 == -PI || x1 == PI {
        0.0
    } else {
        -1.0 / (a * (1.0 - a2).sqrt()) * ((a - (x1 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
            + ((a + (x1 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
    };

    val1 - val0
}
