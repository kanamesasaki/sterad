use crate::diff_element_to_disk;
use crate::error::ViewFactorError;
use nalgebra::{Matrix3, Vector3};
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
///
/// /// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_ellipse;
/// let h = 2.0;
/// let a = 1.0;
/// let b = 1.0;
/// let vf = diff_element_to_ellipse::parallel_center(h, a, b).unwrap();
/// assert!((vf - 0.2).abs() < 1e-8);
/// ```
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

/// B-18a: Tilted planar element to elliptical plate.
///
/// # Arguments
///
/// * `h` - Distance between the plate element and the ellipse.
/// * `a` - Semi-major axis of the ellipse.
/// * `b` - Semi-minor axis of the ellipse.
/// * `theta` - Polar angle of the plate element.
/// * `phi` - Azimuthal angle of the plate element.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the ellipse from the plate element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `a` or `b` is less than or equal to 0.
/// * `theta` is not between 0 and PI.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_ellipse;
/// let h = 2.0;
/// let a = 1.0;
/// let b = 1.0;
/// let theta = 0.0;
/// let phi = 0.3;
/// let vf = diff_element_to_ellipse::tilted_center(h, a, b, theta, phi).unwrap();
/// assert!((vf - 0.2).abs() < 1e-8);
/// ```
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

    // disk view factor
    if a == b {
        return diff_element_to_disk::tilted_center(h, a, theta);
    }
    // ellipse parallel view factor
    if theta == 0.0 {
        return parallel_center(h, a, b);
    }

    // set the semi-major axis and semi-minor axis
    let (a, b, mut phi) = if a > b {
        (a, b, phi)
    } else {
        (b, a, phi + PI / 2.0)
    };

    // transform phi to 0 <= phi <= PI/2
    phi %= 2.0 * PI;
    phi = if phi > PI / 2.0 && phi <= PI {
        PI - phi
    } else if phi > PI && phi <= PI * 3.0 / 2.0 {
        phi - PI
    } else if phi > PI * 3.0 / 2.0 && phi <= 2.0 * PI {
        2.0 * PI - phi
    } else {
        phi
    };

    let theta_tangent: f64 =
        (h / (b.powi(2) * phi.sin().powi(2) + a.powi(2) * phi.cos().powi(2)).sqrt()).atan();
    if theta <= theta_tangent + 1.0E-10 {
        Ok(tilted_center_full(h, a, b, theta))
    } else if theta < PI - theta_tangent - 1.0E-10 {
        tilted_center_partial(h, a, b, theta, phi)
    } else {
        Ok(0.0)
    }
}

fn tilted_center_full(h: f64, a: f64, b: f64, theta: f64) -> f64 {
    a * b * theta.cos() / ((h.powi(2) + a.powi(2)) * (h.powi(2) + b.powi(2))).sqrt()
}

fn tilted_center_partial(
    h: f64,
    a: f64,
    b: f64,
    theta: f64,
    phi: f64,
) -> Result<f64, ViewFactorError> {
    if !(0.0..=PI / 2.0).contains(&phi) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "phi",
            message: "phi must be between 0 and PI/2".to_string(),
        });
    }
    if a <= b {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a, b",
            message: "a must be greater than b".to_string(),
        });
    }

    let (x0, y0, x1, y1) = if (phi - 0.0).abs() < 1e-10 {
        let x0 = -h * theta.cos() / theta.sin();
        let y0 = b / a * (a.powi(2) - x0.powi(2)).sqrt();
        (x0, y0, x0, -y0)
    } else if (phi - PI / 2.0).abs() < 1e-10 {
        let y0 = -h * theta.cos() / theta.sin();
        let x0 = -a / b * (b.powi(2) - y0.powi(2)).sqrt();
        (x0, y0, -x0, y0)
    } else {
        let param_a = 1.0 / a.powi(2) + phi.cos().powi(2) / (b.powi(2) * phi.sin().powi(2));
        let param_b = h * theta.cos() * phi.cos() / (b.powi(2) * theta.sin() * phi.sin().powi(2));
        let param_c = h.powi(2) * theta.cos().powi(2)
            / (b.powi(2) * theta.sin().powi(2) * phi.sin().powi(2))
            - 1.0;
        let x0 = (-param_b - (param_b.powi(2) - param_a * param_c).sqrt()) / param_a;
        let x1 = (-param_b + (param_b.powi(2) - param_a * param_c).sqrt()) / param_a;
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
        // line integration sequence
        // vf1. ellipse edge: PI -> ang1
        // vf2. straight line: p1 -> p0
        // vf3. ellipse edge: ang0 -> -PI
        let a2 = (a.powi(2) - b.powi(2)) / (a.powi(2) + h.powi(2));
        let vf1 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_cosx(PI, ang1, a2)?
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(PI, ang1, a2)?
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(PI, ang1, a2)?;
        let param_k = (x0 - x1).powi(2) + (y0 - y1).powi(2);
        let param_l = (x0 - x1) * x1 + (y0 - y1) * y1;
        let param_m = x1.powi(2) + y1.powi(2) + h.powi(2);
        let param_n = l * h * (y0 - y1) - m * h * (x0 - x1) + n * (y1 * x0 - x1 * y0);
        let vf2 = param_n / (2.0 * PI * (-param_l.powi(2) + param_k * param_m).sqrt())
            * ((param_k + param_l) / (-param_l.powi(2) + param_k * param_m).sqrt()).atan()
            - (param_l / (-param_l.powi(2) + param_k * param_m).sqrt()).atan();

        let vf3 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_cosx(ang0, -PI, a2)?
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(ang0, -PI, a2)?
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(ang0, -PI, a2)?;

        vf1 + vf2 + vf3
    } else {
        // line integration sequence:
        // vf1. ellipse edge: ang0 -> ang1
        // vf2. straight line: p1 -> p0
        let a2 = (a.powi(2) - b.powi(2)) / (a.powi(2) + h.powi(2));
        let vf1 = l * h * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_cosx(ang0, ang1, a2)?
            + m * h * a / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_sinx(ang0, ang1, a2)?
            - n * a * b / (2.0 * PI * (a.powi(2) + h.powi(2))) * integral_1(ang0, ang1, a2)?;
        let param_k = (x0 - x1).powi(2) + (y0 - y1).powi(2);
        let param_l = (x0 - x1) * x1 + (y0 - y1) * y1;
        let param_m = x1.powi(2) + y1.powi(2) + h.powi(2);
        let param_n = l * h * (y0 - y1) - m * h * (x0 - x1) + n * (y1 * x0 - x1 * y0);
        let vf2 = param_n / (2.0 * PI * (-param_l.powi(2) + param_k * param_m).sqrt())
            * ((param_k + param_l) / (-param_l.powi(2) + param_k * param_m).sqrt()).atan()
            - (param_l / (-param_l.powi(2) + param_k * param_m).sqrt()).atan();

        vf1 + vf2
    };

    Ok(vf)
}

/// Calculate the integral of 1 / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -PI <= x0 <= PI.
/// * `x1` - Upper limit of the integral, -PI <= x1 <= PI.
/// * `a2` - a^2, 0 < a2 < 1.
///
/// # Returns
///
/// The integral value as f64.
///
fn integral_1(x0: f64, x1: f64, a2: f64) -> Result<f64, ViewFactorError> {
    if a2 <= 0.0 || a2 >= 1.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a2",
            message: "a2 must be in the range: 0 < a2 < 1".to_string(),
        });
    }

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

    Ok(val1 - val0)
}

/// Calculate the integral of cos(x) / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -PI <= x0 <= PI.
/// * `x1` - Upper limit of the integral, -PI <= x1 <= PI.
/// * `a2` - a^2, 0 < a2 < 1.
///
/// # Returns
///
/// The integral value.
fn integral_cosx(x0: f64, x1: f64, a2: f64) -> Result<f64, ViewFactorError> {
    if a2 <= 0.0 || a2 >= 1.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a2",
            message: "a2 must be in the range: 0 < a2 < 1".to_string(),
        });
    }

    let a: f64 = a2.sqrt();
    let val0: f64 = (a * x0.sin()).atanh() / a;
    let val1: f64 = (a * x1.sin()).atanh() / a;
    Ok(val1 - val0)
}

/// Calculate the integral of sin(x) / (1 - a^2 sin^2(x)) from x0 to x1
///
/// # Arguments
///
/// * `x0` - Lower limit of the integral, -PI <= x0 <= PI.
/// * `x1` - Upper limit of the integral, -PI <= x1 <= PI.
/// * `a2` - a^2, 0 < a2 < 1.
///
/// # Returns
///
/// The integral value.
fn integral_sinx(x0: f64, x1: f64, a2: f64) -> Result<f64, ViewFactorError> {
    if a2 <= 0.0 || a2 >= 1.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a2",
            message: "a2 must be in the range: 0 < a2 < 1".to_string(),
        });
    }

    let a: f64 = a2.sqrt();
    let val0: f64 = if x0 == -PI || x0 == PI {
        0.0
    } else {
        -1.0 / (a * (1.0 - a2).sqrt()) * ((a - (x0 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
            + ((a + (x0 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
    };

    let val1: f64 = if x1 == -PI || x1 == PI {
        0.0
    } else {
        -1.0 / (a * (1.0 - a2).sqrt()) * ((a - (x1 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
            + ((a + (x1 / 2.0).tan()) / (1.0 - a2).sqrt()).atan()
    };

    Ok(val1 - val0)
}

/// B-18b: Tilted planar element to elliptical plate with offset.
///
/// # Arguments
///
/// * `a` - Semi-major axis of the ellipse.
/// * `b` - Semi-minor axis of the ellipse.
/// * `xc` - X-coordinate of the ellipse center.
/// * `yc` - Y-coordinate of the ellipse center.
/// * `zc` - Z-coordinate of the ellipse center.
/// * `theta` - Polar angle of the plate element.
/// * `phi` - Azimuthal angle of the plate element.
///
/// # Returns
///
/// * `Result<f64, ViewFactorError>` - The view factor of the ellipse from the plate element, or an error if the input parameters are invalid.
///
/// # Errors
///
/// This function will return a `ViewFactorError` if:
///
/// * `a` or `b` is less than or equal to 0.
/// * `theta` is not between 0 and PI.
///
/// # Examples
///
/// ```
/// use sterad_view_factor::diff_element_to_ellipse;
/// let a = 1.0;
/// let b = 1.0;
/// let xc = 0.0;
/// let yc = 0.0;
/// let zc = 2.0;
/// let theta = 0.0;
/// let phi = 0.3;
/// let result = diff_element_to_ellipse::tilted_offset(a, b, xc, yc, zc, theta, phi);
/// assert!((result.unwrap() - 0.2).abs() < 1e-8);
/// ```
pub fn tilted_offset(
    a: f64,
    b: f64,
    xc: f64,
    yc: f64,
    zc: f64,
    theta: f64,
    phi: f64,
) -> Result<f64, ViewFactorError> {
    if a <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a",
            message: "a must be greater than 0".to_string(),
        });
    }
    if b <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "b",
            message: "b must be greater than 0".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let a2 = a * a;
    let b2 = b * b;
    let zc2 = zc * zc;
    let xc2 = xc * xc;
    let yc2 = yc * yc;

    let param_a = 1.0 - xc2 / a2 - yc2 / b2 - zc2 / a2 - zc2 / b2;
    let param_b = (xc2 + yc2 + zc2) * zc2 / (a2 * b2) - zc2 / a2 - zc2 / b2;
    let param_c = zc2 * zc2 / (a2 * b2);

    let p = param_b - param_a * param_a / 3.0;
    let q = param_c - param_a * param_b / 3.0 + 2.0 * param_a.powi(3) / 27.0;

    let discriminant = -4.0 * p.powi(3) - 27.0 * q.powi(2);
    if discriminant <= 0.0 {
        return Err(ViewFactorError::RuntimeError {
            message: "Discriminant has to be larger than zero.".to_string(),
        });
    }

    let acos_triple = (3.0 * q / (2.0 * p) * (-3.0 / p).sqrt()).acos() / 3.0;
    let lambda_1 = 2.0 * (-p / 3.0).sqrt() * (acos_triple).cos() - param_a / 3.0;
    let lambda_2 = 2.0 * (-p / 3.0).sqrt() * (acos_triple + 2.0 * PI / 3.0).cos() - param_a / 3.0;
    let lambda_3 = 2.0 * (-p / 3.0).sqrt() * (acos_triple + 4.0 * PI / 3.0).cos() - param_a / 3.0;

    let n2: Vector3<f64> = Vector3::new(
        -xc * zc / a2 * (lambda_2 - zc2 / b2),
        -yc * zc / b2 * (lambda_2 - zc2 / a2),
        (lambda_2 - zc2 / b2) * (lambda_2 - zc2 / a2),
    );
    let n3: Vector3<f64> = Vector3::new(
        -xc * zc / a2 * (lambda_3 - zc2 / b2),
        -yc * zc / b2 * (lambda_3 - zc2 / a2),
        (lambda_3 - zc2 / b2) * (lambda_3 - zc2 / a2),
    );

    let nz: Vector3<f64> = n2.normalize();
    let nx: Vector3<f64> = n3.normalize();
    let ny: Vector3<f64> = nz.cross(&nx);
    let m_trans: Matrix3<f64> = Matrix3::from_columns(&[nx, ny, nz]);

    let v: Vector3<f64> = Vector3::new(
        theta.sin() * phi.cos(),
        theta.sin() * phi.sin(),
        theta.cos(),
    );
    let v_new: Vector3<f64> = m_trans * v;
    let theta_new = v_new.z.acos();
    let phi_new = v_new.y.atan2(v_new.x);
    let a_new = (-lambda_2 / lambda_3).sqrt();
    let b_new = (-lambda_2 / lambda_1).sqrt();
    tilted_center(1.0, a_new, b_new, theta_new, phi_new)
}
