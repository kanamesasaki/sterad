use crate::diff_element_to_ellipse;
use crate::error::ViewFactorError;
use crate::vecmath::{Matrix3f, Vector3f};
use std::f64::consts::PI;

#[allow(clippy::too_many_arguments)]
pub fn tilted_offset(
    a: f64,
    b: f64,
    c: f64,
    xc: f64,
    yc: f64,
    zc: f64,
    theta: f64,
    phi: f64,
) -> Result<f64, ViewFactorError> {
    if a <= 0.0 || b <= 0.0 || c <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "a, b, c",
            message: "a, b, and c must be greater than 0".to_string(),
        });
    }
    if !(0.0..=PI).contains(&theta) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "theta",
            message: "theta must be between 0 and PI".to_string(),
        });
    }

    let a2 = a.powi(2);
    let b2 = b.powi(2);
    let c2 = c.powi(2);
    let xc2 = xc.powi(2);
    let yc2 = yc.powi(2);
    let zc2 = zc.powi(2);

    let param_a = (xc2 + yc2) / (a2 * b2) + (yc2 + zc2) / (b2 * b2) + (xc2 + zc2) / (a2 * b2)
        - 1.0 / a2
        - 1.0 / b2
        - 1.0 / b2;
    let param_b =
        (xc2 + yc2 + zc2 - a2 - b2 - b2) / (a2 * b2 * b2) * (xc2 / a2 + yc2 / b2 + zc2 / b2 - 1.0);
    let param_c = -(-a2 * b2 * b2 + a2 * b2 * zc2 + a2 * b2 * yc2 + b2 * b2 * xc2).powi(2)
        / (a2.powi(3) * b2.powi(3) * c2.powi(3));

    let p = param_b - param_a.powi(2) / 3.0;
    let q = param_c - param_a * param_b / 3.0 + 2.0 * param_a.powi(3) / 27.0;
    let discriminant = -4.0 * p.powi(3) - 27.0 * q.powi(2);
    if discriminant < 0.0 {
        return Err(ViewFactorError::RuntimeError {
            message: "Discriminant has to be zero or greater.".to_string(),
        });
    }

    let lambda_1 = 2.0
        * (-p / 3.0).sqrt()
        * (1.0 / 3.0 * (3.0 * q / (2.0 * p) * (-3.0 / p).sqrt()).acos()).cos()
        - param_a / 3.0;
    let lambda_2 = 2.0
        * (-p / 3.0).sqrt()
        * (1.0 / 3.0 * (3.0 * q / (2.0 * p) * (-3.0 / p).sqrt()).acos() + 2.0 * PI / 3.0).cos()
        - param_a / 3.0;
    let lambda_3 = 2.0
        * (-p / 3.0).sqrt()
        * (1.0 / 3.0 * (3.0 * q / (2.0 * p) * (-3.0 / p).sqrt()).acos() + 4.0 * PI / 3.0).cos()
        - param_a / 3.0;

    let vec_center = Vector3f::new(xc, yc, zc);
    let mut n1 = Vector3f::new(
        xc * yc2 * zc / (a2 * b2.powi(2) * c2)
            + xc * zc / (a2 * c2) * (lambda_1 + xc2 / (a2 * b2) + zc2 / (b2 * c2) - 1.0 / b2),
        xc2 * yc * zc / (a2.powi(2) * b2 * c2)
            + yc * zc / (b2 * c2) * (lambda_1 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
        (lambda_1 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2)
            * (lambda_1 + xc2 / (a2 * b2) + zc2 / (b2 * c2) - 1.0 / b2)
            - xc2 * yc2 / (a2.powi(2) * b2.powi(2)),
    );
    if n1.norm() < 1e-8 {
        n1 = Vector3f::new(
            -xc * yc * zc2 / (a2 * b2 * c2.powi(2))
                - xc * yc / (a2 * b2) * (lambda_1 + xc2 / (a2 * c2) + yc2 / (b2 * c2) - 1.0 / c2),
            xc2 * zc2 / (a2.powi(2) * c2.powi(2))
                - (lambda_1 + xc2 / (a2 * c2) + yc2 / (b2 * c2) - 1.0 / c2)
                    * (lambda_1 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
            -xc2 * yc * zc / (a2.powi(2) * b2 * c2)
                - yc * zc / (b2 * c2) * (lambda_1 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
        );
    }
    if n1.dot(&vec_center) < 0.0 {
        n1 = -n1;
    }

    let mut n3 = Vector3f::new(
        xc * yc2 * zc / (a2 * b2.powi(2) * c2)
            + xc * zc / (a2 * c2) * (lambda_3 + xc2 / (a2 * b2) + zc2 / (b2 * c2) - 1.0 / b2),
        xc2 * yc * zc / (a2.powi(2) * b2 * c2)
            + yc * zc / (b2 * c2) * (lambda_3 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
        (lambda_3 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2)
            * (lambda_3 + xc2 / (a2 * b2) + zc2 / (b2 * c2) - 1.0 / b2)
            - xc2 * yc2 / (a2.powi(2) * b2.powi(2)),
    );
    if n3.norm() < 1e-8 {
        n3 = Vector3f::new(
            -xc * yc * zc2 / (a2 * b2 * c2.powi(2))
                - xc * yc / (a2 * b2) * (lambda_3 + xc2 / (a2 * c2) + yc2 / (b2 * c2) - 1.0 / c2),
            xc2 * zc2 / (a2.powi(2) * c2.powi(2))
                - (lambda_3 + xc2 / (a2 * c2) + yc2 / (b2 * c2) - 1.0 / c2)
                    * (lambda_3 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
            -xc2 * yc * zc / (a2.powi(2) * b2 * c2)
                - yc * zc / (b2 * c2) * (lambda_3 + yc2 / (a2 * b2) + zc2 / (a2 * c2) - 1.0 / a2),
        );
    }

    let nz: Vector3f = n1.normalize();
    let nx: Vector3f = n3.normalize();
    let ny: Vector3f = nz.cross(&nx);
    let m_trans: Matrix3f = Matrix3f::from_rows(&nx, &ny, &nz);

    let v: Vector3f = Vector3f::new(
        theta.sin() * phi.cos(),
        theta.sin() * phi.sin(),
        theta.cos(),
    );
    let v_new: Vector3f = m_trans.mul_vector(&v);
    let theta_new = v_new.z.acos();
    let phi_new = v_new.y.atan2(v_new.x);
    let a_new = (-lambda_1 / lambda_3).sqrt();
    let b_new = (-lambda_1 / lambda_2).sqrt();
    let h_new = 1.0;

    diff_element_to_ellipse::tilted_center(h_new, a_new, b_new, theta_new, phi_new)
}
