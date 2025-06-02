use crate::diff_element_to_ellipse;
use crate::error::ViewFactorError;

#[allow(clippy::too_many_arguments)]
pub fn tilted_offset(
    a: f64,
    b: f64,
    x0: f64,
    y0: f64,
    z0: f64,
    theta: f64,
    phi: f64,
) -> Result<f64, ViewFactorError> {
    if x0.powi(2) / a.powi(2) + y0.powi(2) / a.powi(2) + z0.powi(2) / b.powi(2) < 1.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "x0, y0, z0",
            message: "The point is inside the spheroid".to_string(),
        });
    }

    // Shift the plate element to XZ plane
    let phi: f64 = phi - y0.atan2(x0);
    let x0: f64 = (x0 * x0 + y0 * y0).sqrt();

    // Shift the plate element to z0 ≥ 0
    let theta: f64 = if z0 >= 0.0 { theta } else { -theta };
    let z0: f64 = z0.abs();

    let (cos_gamma0, sin_gamma0) = if z0 == 0.0 {
        // Special case when plate element is on the x-y plane
        (1.0, 0.0)
    } else if x0 == 0.0 {
        // Special case when plate element is on the y-z plane
        (0.0, 1.0)
    } else {
        // General case
        let tan_2gamma0 = 2.0 * x0 * z0 / (-a.powi(2) + x0.powi(2) + b.powi(2) - z0.powi(2));
        let tan_1gamma0 = -1.0 / tan_2gamma0 + (1.0 / tan_2gamma0.powi(2) + 1.0).sqrt();
        let cos_gamma02 = 1.0 / (1.0 + tan_1gamma0.powi(2));
        let sin_gamma02 = tan_1gamma0.powi(2) / (1.0 + tan_1gamma0.powi(2));
        let cos_gamma0 = cos_gamma02.sqrt();
        let sin_gamma0 = sin_gamma02.sqrt();

        (cos_gamma0, sin_gamma0)
    };

    // Common calculations
    let denom = -a.powi(2) * b.powi(2) + b.powi(2) * x0.powi(2) + a.powi(2) * z0.powi(2);
    let tan_gamma_delta = (a.powi(2) - x0.powi(2) + b.powi(2) - z0.powi(2))
        / (2.0 * (-a.powi(2) * b.powi(2) + b.powi(2) * x0.powi(2) + a.powi(2) * z0.powi(2)).sqrt())
        + ((a.powi(2) - x0.powi(2) + b.powi(2) - z0.powi(2)).powi(2)
            / (4.0 * (-a.powi(2) * b.powi(2) + b.powi(2) * x0.powi(2) + a.powi(2) * z0.powi(2)))
            + 1.0)
            .sqrt();

    let cs = 1.0 / b.powi(2)
        * (-b.powi(2) * cos_gamma0.powi(2) - a.powi(2) * sin_gamma0.powi(2)
            + b.powi(4) * x0.powi(2) * cos_gamma0.powi(2) / denom
            + a.powi(4) * z0.powi(2) * sin_gamma0.powi(2) / denom
            + 2.0 * a.powi(2) * b.powi(2) * cos_gamma0 * sin_gamma0 * x0 * z0 / denom);

    let theta_ellipse = (-cos_gamma0 * theta.cos() * phi.cos() - sin_gamma0 * theta.sin()).acos();
    let phi_ellipse = (-sin_gamma0 * theta.cos() * phi.cos() + cos_gamma0 * theta.sin())
        .atan2(theta.cos() * phi.sin());

    diff_element_to_ellipse::tilted_center(
        1.0,
        cs.sqrt(),
        tan_gamma_delta,
        theta_ellipse,
        phi_ellipse,
    )
}
