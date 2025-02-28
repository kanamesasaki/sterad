// This module implements the view factor calculation from a differential element to a spherical cap.
// The implementation is based on the paper:
// "Analytical view factor solutions of a spherical cap from an infinitesimal surface"
// https://doi.org/10.1016/j.ijheatmasstransfer.2020.120477

// To build the ducumentation including the private items, run: cargo doc --document-private-items --open

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_else_if)]

use crate::error::ViewFactorError;
use crate::vecmath::Vector3f;
use std::f64::consts::{FRAC_PI_2, PI};

///////////////////////////////////////////////////////////////////////////////
// Case Description
///////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////
// Line Integral Calculation
///////////////////////////////////////////////////////////////////////////////

/// Calculates the L1 integration
///
/// <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
/// L1 is the line integration along the projection plane edge as described by
/// $$
/// \begin{equation}
/// \left[ \begin{array}{c} x \\\\ y \\\\ z \end{array} \right] =
/// \left[ \begin{array}{c} R \cos \beta \\\\ R \sin \beta \\\\ h \end{array} \right].
/// \end{equation}
/// $$
/// The analytical expression for L1 is given by
/// $$
/// \begin{equation}
/// L_1(\beta_0, \beta_1) = l_1 \frac{hR(\sin \beta_1 - \sin \beta_0)}{2\pi (R^2 + h^2)} + m_1 \frac{hR(- \cos \beta_1 + \cos \beta_0)}{2\pi (R^2 + h^2)} + n_1 \frac{R^2 (-\beta_1 + \beta_0)}{2\pi (R^2 + h^2)}.
/// \end{equation}
/// $$
///
/// # Arguments
///
/// * `beta0` - Starting angle of the arc segment in radians
/// * `beta1` - Ending angle of the arc segment in radians  
/// * `omega` - Polar angle of the differential element's normal vector in radians
/// * `d` - Distance from differential element to center of sphere
/// * `rs` - Radius of the sphere
/// * `gamma` - Azimuthal angle of the differential element's normal vector in radians
///
/// # Returns
///
/// The L1 integration value
///
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

/// Calculates the L2 integration
///
/// <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
/// L2 represents the line integration along the intersection curve between the spherical cap
/// and the projection plane. The intersection curve is described parametrically using angle α:
/// $$
/// \begin{equation}
///     \left[ \begin{array}{c} x \\\\ y \\\\ z \end{array} \right] =
///     \left[ \begin{array}{c}
///     -R_\mathrm{S} \cos \psi \sin \varphi + R_\mathrm{S} \cos \alpha \sin \psi \cos \varphi \\\\
///     R_\mathrm{S} \sin \alpha \sin \psi \\\\
///     d - R_\mathrm{S} \cos \psi \cos \varphi - R_\mathrm{S} \cos \alpha \sin \psi \sin \varphi
/// \end{array} \right]
/// \end{equation}
/// $$
/// The analytical expression for L2 is given by
/// $$
/// \begin{align}
///     L_2(\alpha_0, \alpha_1)
///     &= \frac{l_1 A}{4\pi \sin \theta \sin \psi \sin \varphi} \int_{\alpha_0}^{\alpha_1} \frac{l_2 \cos \alpha - l_3}{1-A \cos \alpha}d\alpha
///     + \frac{m_1 A}{4\pi \sin \theta \sin \psi \sin \varphi} \int_{\alpha_0}^{\alpha_1} \frac{m_2 \sin \alpha}{1-A \cos \alpha}d\alpha \\\\
///     &+ \frac{n_1 A}{4\pi \sin \theta \sin \psi \sin \varphi} \int_{\alpha_0}^{\alpha_1} \frac{n_2 \cos \alpha - n_3}{1-A \cos \alpha}d\alpha.
/// \end{align}
/// $$
///
/// # Arguments
///
/// * `alpha0` - Starting angle of the intersection curve in radians (-π ≤ alpha0 ≤ π)
/// * `alpha1` - Ending angle of the intersection curve in radians (-π ≤ alpha1 ≤ π)
/// * `omega` - Polar angle of the differential element's normal vector in radians
/// * `d` - Distance from differential element to center of sphere
/// * `rs` - Radius of the sphere
/// * `phi` - Polar angle of the spherical cap's axis in radians
/// * `gamma` - Azimuthal angle of the differential element's normal vector in radians
/// * `psi` - Half-angle of the spherical cap in radians
///
/// # Returns
///
/// * `Ok(f64)` - The L2 integration value
/// * `Err(ViewFactorError)` - If input parameters are invalid
///
/// # Errors
///
/// Returns `ViewFactorError::InvalidInput` if:
/// * `alpha0` is not between -π and π
/// * `alpha1` is not between -π and π
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

/// Calculates a part of the L2 integration
///
/// <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
/// This is a utility function, calculating an arctangent term in the L2 integration.
/// \begin{equation}
///    \arctan \left( \sqrt{\frac{1 + A}{1 - A}} \tan\frac{\alpha}{2} \right)
/// \end{equation}
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

/// Calculates the L4 integration term for the view factor computation.
///
/// <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
/// L4 represents the line integration along a straight line segment between two points on the projection plane.
/// The corresponding line segment is described by
/// $$
/// \begin{equation}
///     \left[ \begin{array}{c} x \\\\ y \\\\ z \end{array} \right]
///     = \left[ \begin{array}{c} x_0 + x_{\Delta}t \\\\ y_0 + y_{\Delta}t \\\\ h \end{array} \right],
/// \end{equation}
/// $$
/// where:
/// $$
/// \begin{equation}
///     x_{\Delta} = x_1 - x_0, \quad y_{\Delta} = y_1 - y_0, \quad t \in [0, 1].
/// \end{equation}
/// $$
///
/// The analytical expression of the L4 integration is given by
/// $$
/// \begin{align}
///     L_4(\boldsymbol{x_0}, \boldsymbol{x_1})
///     =& l_1 \frac{hy_{\Delta}}{2\pi \sqrt{AC - B^2}} \left\\{ \arctan \left( \frac{A+B}{\sqrt{AC - B^2}} \right) - \arctan \left( \frac{B}{\sqrt{AC - B^2}} \right) \right\\} \\\\
///     &+ m_1 \frac{-hx_{\Delta}}{2\pi \sqrt{AC - B^2}} \left\\{ \arctan \left( \frac{A+B}{\sqrt{AC - B^2}} \right) - \arctan \left( \frac{B}{\sqrt{AC - B^2}} \right) \right\\} \\\\
///     &+ n_1 \frac{y_0 x_{\Delta} - x_0 y_{\Delta}}{2\pi \sqrt{AC - B^2}} \left\\{ \arctan \left( \frac{A+B}{\sqrt{AC - B^2}} \right)  - \arctan \left( \frac{B}{\sqrt{AC - B^2}} \right) \right\\},
/// \end{align}
/// $$
/// where:
/// $$
/// \begin{equation}
///     A = x_{\Delta}^2 + y_{\Delta}^2, \quad
///     B = x_0 x_{\Delta} + y_0 y_{\Delta}, \quad
///     C = x_0^2 + y_0^2 + h^2.
/// \end{equation}
/// $$
///
/// # Arguments
///
/// * `x0` - Start point of the line segment
/// * `x1` - End point of the line segment
/// * `omega` - Polar angle of the differential element's normal vector in radians
/// * `gamma` - Azimuthal angle of the differential element's normal vector in radians
///
/// # Returns
///
/// * `Ok(f64)` - The L4 integration value
/// * `Err(ViewFactorError)` - If the input points have different z-coordinates
///
/// # Errors
///
/// Returns `ViewFactorError::InvalidInput` if:
/// * `x0.z` and `x1.z` are not equal (points must lie in same plane)
fn l4(x0: &Vector3f, x1: &Vector3f, omega: f64, gamma: f64) -> Result<f64, ViewFactorError> {
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

///////////////////////////////////////////////////////////////////////////////
// Spherical View Factor Calculation
///////////////////////////////////////////////////////////////////////////////

/// Calculates the view factor from a differential element to a spherical cap
///
/// This function serves as the main entry point for spherical cap view factor calculations,
/// handling both the standard case (0 < psi <= π/2) and the extended case (π/2 < psi < π).
///
/// # Arguments
///
/// * `omega` - Polar angle of the differential element's normal vector in radians
/// * `d` - Distance from differential element to center of sphere
/// * `rs` - Radius of the sphere
/// * `phi` - Polar angle of spherical cap's axis in radians
/// * `gamma` - Azimuthal angle of differential element's normal vector in radians
/// * `psi` - Half-angle of the spherical cap in radians
///
/// # Returns
///
/// * `Ok((f64, i32))` - Tuple with (view factor, case number)
/// * `Err(ViewFactorError)` - If input parameters are invalid
///
pub fn sphericalcap(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    // Input validation
    if !(0.0..=PI).contains(&omega) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "omega",
            message: "omega must be between 0 and π".to_string(),
        });
    }

    if d <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "d",
            message: "d must be positive".to_string(),
        });
    }

    if rs <= 0.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "rs",
            message: "rs must be positive".to_string(),
        });
    }

    if !(0.0..=PI).contains(&phi) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "phi",
            message: "phi must be between 0 and π".to_string(),
        });
    }

    if !(-PI..=PI).contains(&gamma) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "gamma",
            message: "gamma must be between -π and π".to_string(),
        });
    }

    if !(0.0..PI).contains(&psi) {
        return Err(ViewFactorError::InvalidInput {
            param_name: "psi",
            message: "psi must be between 0 and π (exclusive)".to_string(),
        });
    }

    if rs >= d {
        return Err(ViewFactorError::InvalidInput {
            param_name: "rs, d",
            message: "rs must be less than d (differential element must be outside the sphere)"
                .to_string(),
        });
    }

    // For standard half-angle cases (0 < psi <= π/2)
    if psi <= FRAC_PI_2 {
        sphericalcap_half(omega, d, rs, phi, gamma, psi)
    } else {
        // For extended half-angle cases (π/2 < psi < π)
        // Calculate the complementary angle (π - psi)
        let psi_complementary = PI - psi;

        // Calculate view factor for the complementary spherical cap
        let (vf_complementary, case_complementary) =
            sphericalcap_half(omega, d, rs, phi + PI, gamma, psi_complementary)?;

        // Calculate view factor for a full sphere
        let vf_sphere = sphere(omega, d, rs, gamma)?;

        // The view factor is: full sphere minus complementary spherical cap
        Ok((vf_sphere - vf_complementary, case_complementary))
    }
}

/// Calculates view factor from a differential element to a full sphere
///
/// # Arguments
///
/// * `omega` - Polar angle of the differential element's normal vector in radians
/// * `d` - Distance from differential element to center of sphere
/// * `rs` - Radius of the sphere
/// * `gamma` - Azimuthal angle of differential element's normal vector in radians
///
/// # Returns
///
/// * `Ok(f64)` - The view factor value
/// * `Err(ViewFactorError)` - If input parameters are invalid
fn sphere(omega: f64, d: f64, rs: f64, gamma: f64) -> Result<f64, ViewFactorError> {
    let sin_theta = rs / d;

    if sin_theta >= 1.0 {
        return Err(ViewFactorError::InvalidInput {
            param_name: "rs, d",
            message: "rs must be less than d (differential element must be outside the sphere)"
                .to_string(),
        });
    }

    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let theta = cos_theta.acos();

    if omega <= FRAC_PI_2 - theta || (omega + theta - FRAC_PI_2).abs() < f64::EPSILON {
        // Full sphere is visible
        Ok(sin_theta.powi(2))
    } else if omega >= FRAC_PI_2 + theta || (omega - theta - FRAC_PI_2).abs() < f64::EPSILON {
        // Full sphere is not visible
        Ok(0.0)
    } else {
        // Partially visible sphere
        let vf = f4(omega, d, rs, gamma)?;
        Ok(vf)
    }
}

/// Returns the spherical cap view factor for 0 < psi <= π/2
///
/// # Arguments
///
/// * `omega` - Direction of the infinitesimal surface in radians
/// * `d` - Distance from the infinitesimal surface to the sphere center
/// * `rs` - Radius of the sphere
/// * `phi` - Direction of the spherical cap in radians
/// * `gamma` - Direction of the infinitesimal surface in radians
/// * `psi` - Size of the spherical cap in radians
///
/// # Returns
///
/// * `Ok((f64, i32))` - Tuple with (view factor, case number)
/// * `Err(ViewFactorError)` - If input parameters are invalid
///
fn sphericalcap_half(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let theta = cos_theta.acos();

    if phi - psi >= FRAC_PI_2 - theta || (phi - psi + theta - FRAC_PI_2).abs() < f64::EPSILON {
        // Case 1: cap orientation mismatch
        Ok((0.0, 1))
    } else if omega >= theta + FRAC_PI_2 || (omega - theta - FRAC_PI_2).abs() < f64::EPSILON {
        // Case 2: surface direction mismatch
        Ok((0.0, 2))
    } else if psi - phi + theta >= FRAC_PI_2 || (psi - phi + theta - FRAC_PI_2).abs() < f64::EPSILON
    {
        // Case 3 and Case 4
        if omega <= FRAC_PI_2 - theta || (omega + theta - FRAC_PI_2).abs() < f64::EPSILON {
            // Case 3: full sphere
            let vf = f3(omega, d, rs, gamma);
            Ok((vf, 3))
        } else {
            // Case 4: partial sphere
            let vf = f4(omega, d, rs, gamma)?;
            Ok((vf, 4))
        }
    } else if phi + psi + theta <= FRAC_PI_2 {
        // Case 16 - Case 21: small cap
        sphericalcap_small(omega, d, rs, phi, gamma, psi)
    } else {
        // Case 5 - Case 15
        if omega <= FRAC_PI_2 - theta || (omega + theta - FRAC_PI_2).abs() < f64::EPSILON {
            // Case 5: full cap
            let vf = f5(omega, d, rs, phi, gamma, psi)?;
            Ok((vf, 5))
        } else {
            // Case 6 - Case 15
            sphericalcap_part(omega, d, rs, phi, gamma, psi)
        }
    }
}

/// Handles cases 6-15 for partial spherical cap view factor calculation
///
/// This function determines the appropriate case and calculates the view factor
/// for a partial spherical cap based on the geometry.
///
/// # Arguments
///
/// * `omega` - Polar angle of differential element's normal vector in radians
/// * `d` - Distance from differential element to center of sphere
/// * `rs` - Radius of the sphere
/// * `phi` - Polar angle of spherical cap's axis in radians
/// * `gamma` - Azimuthal angle of differential element's normal vector in radians
/// * `psi` - Half-angle of the spherical cap in radians
///
/// # Returns
///
/// * `Ok((f64, i32))` - Tuple with (view factor, case number)
/// * `Err(ViewFactorError)` - If input parameters are invalid
fn sphericalcap_part(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let cos_gamma = gamma.cos();
    let sin_gamma = gamma.sin();
    let cos_psi = psi.cos();
    // let sin_psi = psi.sin();

    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();

    let h = d * cos_theta.powi(2);
    // let r = rs * cos_theta;

    if gamma.abs() < f64::EPSILON {
        // gamma == 0
        if (cos_phi * sin_omega - sin_phi * cos_omega).abs() < f64::EPSILON {
            intersection0(omega, d, rs, phi, gamma, psi)
        } else {
            let z = (-rs * cos_psi + d * cos_phi) * sin_omega
                / (cos_phi * sin_omega - cos_omega * sin_phi);
            let x = (rs * cos_psi - d * cos_phi) * cos_omega
                / (cos_phi * sin_omega - cos_omega * sin_phi);
            let y2 = rs.powi(2) - x.powi(2) - (z - d).powi(2);

            if y2.abs() < f64::EPSILON {
                intersection0(omega, d, rs, phi, gamma, psi)
            } else if z <= h && y2 > 0.0 {
                let w1 = Vector3f { x, y: y2.sqrt(), z };
                let w2 = Vector3f {
                    x,
                    y: -y2.sqrt(),
                    z,
                };
                intersection2(omega, d, rs, phi, gamma, psi, w1, w2)
            } else {
                intersection0(omega, d, rs, phi, gamma, psi)
            }
        }
    } else if (gamma - PI).abs() < f64::EPSILON || (gamma + PI).abs() < f64::EPSILON {
        if (cos_phi * sin_omega + sin_phi * cos_omega).abs() < f64::EPSILON {
            intersection0(omega, d, rs, phi, gamma, psi)
        } else {
            let z = (-rs * cos_psi + d * cos_phi) * sin_omega
                / (cos_phi * sin_omega + cos_omega * sin_phi);
            let x = (-rs * cos_psi + d * cos_phi) * cos_omega
                / (cos_phi * sin_omega + cos_omega * sin_phi);
            let y2 = rs.powi(2) - x.powi(2) - (z - d).powi(2);

            if y2.abs() < f64::EPSILON {
                intersection0(omega, d, rs, phi, gamma, psi)
            } else if z <= h && y2 > 0.0 {
                let w1 = Vector3f { x, y: y2.sqrt(), z };
                let w2 = Vector3f {
                    x,
                    y: -y2.sqrt(),
                    z,
                };
                intersection2(omega, d, rs, phi, gamma, psi, w1, w2)
            } else {
                intersection0(omega, d, rs, phi, gamma, psi)
            }
        }
    } else {
        // Calculate coefficients for the quadratic equation
        let a = (1.0 + cos_gamma.powi(2) / sin_gamma.powi(2)) * cos_phi.powi(2) / sin_phi.powi(2)
            - 2.0 * cos_omega * cos_gamma * cos_phi / (sin_omega * sin_gamma.powi(2) * sin_phi)
            + cos_omega.powi(2) / (sin_omega.powi(2) * sin_gamma.powi(2))
            + 1.0;

        let b = -(1.0 + cos_gamma.powi(2) / sin_gamma.powi(2))
            * (d * cos_phi.powi(2) / sin_phi.powi(2) - rs * cos_phi * cos_psi / sin_phi.powi(2))
            + cos_omega * cos_gamma / (sin_omega * sin_gamma.powi(2))
                * (d * cos_phi / sin_phi - rs * cos_psi / sin_phi)
            - d;

        let c = (1.0 + cos_gamma.powi(2) / sin_gamma.powi(2))
            * (d * cos_phi / sin_phi - rs * cos_psi / sin_phi).powi(2)
            + d.powi(2)
            - rs.powi(2);

        let discriminant = b.powi(2) - a * c;

        if discriminant <= 0.0 || discriminant.abs() < f64::EPSILON {
            intersection0(omega, d, rs, phi, gamma, psi)
        } else {
            let zp = (-b + discriminant.sqrt()) / a;
            let zm = (-b - discriminant.sqrt()) / a;

            if zp <= h && zm <= h {
                let xp = (d - zp) * cos_phi / sin_phi - cos_psi / sin_phi * rs;
                let yp = -cos_gamma / sin_gamma * xp - cos_omega / (sin_omega * sin_gamma) * zp;
                let xm = (d - zm) * cos_phi / sin_phi - cos_psi / sin_phi * rs;
                let ym = -cos_gamma / sin_gamma * xm - cos_omega / (sin_omega * sin_gamma) * zm;

                let w1 = Vector3f {
                    x: xp,
                    y: yp,
                    z: zp,
                };
                let w2 = Vector3f {
                    x: xm,
                    y: ym,
                    z: zm,
                };

                intersection2(omega, d, rs, phi, gamma, psi, w1, w2)
            } else if zp > h && zm > h {
                intersection0(omega, d, rs, phi, gamma, psi)
            } else {
                intersection1(omega, d, rs, phi, gamma, psi, zm)
            }
        }
    }
}

/// Calculates F3 view factor: Full Sphere
fn f3(omega: f64, d: f64, rs: f64, gamma: f64) -> f64 {
    l1(PI, -PI, omega, d, rs, gamma)
}

/// Calculates F4 view factor: Partial Sphere
fn f4(omega: f64, d: f64, rs: f64, gamma: f64) -> Result<f64, ViewFactorError> {
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

/// Calculates F5 view factor: Full Spherical Cap
fn f5(omega: f64, d: f64, rs: f64, phi: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
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

/// Calculates F7 view factor: Partial Cap, No Intersection 2
fn f7(omega: f64, d: f64, rs: f64, phi: f64, gamma: f64, psi: f64) -> Result<f64, ViewFactorError> {
    // return F7 view factor
    let vf5 = f5(omega, d, rs, phi, gamma, psi)?;
    let vf4 = f4(omega, d, rs, gamma)?;
    let vf3 = f3(omega, d, rs, gamma);
    Ok(vf5 + vf4 - vf3)
}

/// Calculates F10 view factor: Partial Cap, One Intersection 1
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

fn intersection0(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_psi = psi.cos();

    let sin_theta = rs / d;
    let cos_theta = (1.0 - sin_theta.powi(2)).sqrt();
    let h = d * cos_theta.powi(2);
    let r = rs * cos_theta;

    let cos_beta0 = -cos_psi * sin_phi / cos_theta
        + (sin_theta - cos_psi * cos_phi) / cos_theta * cos_phi / sin_phi;
    let beta0 = cos_beta0.acos();

    let (vf, case) = if beta0 <= PI / 2.0 {
        if gamma < beta0 && gamma > -beta0 {
            if h * (omega - FRAC_PI_2).tan() >= r * cos_beta0 {
                (0.0, 6)
            } else {
                (f7(omega, d, rs, phi, gamma, psi)?, 7)
            }
        } else if PI - beta0 < gamma || gamma < -PI + beta0 {
            if h * (FRAC_PI_2 - omega).tan() >= r * cos_beta0 {
                (f5(omega, d, rs, phi, gamma, psi)?, 9)
            } else {
                (f4(omega, d, rs, gamma)?, 8)
            }
        } else {
            if h * (omega - FRAC_PI_2).tan() >= r * (gamma.abs() - beta0).cos() {
                (f4(omega, d, rs, gamma)?, 8)
            } else if h * (omega - FRAC_PI_2).tan() <= r * (gamma.abs() + beta0).cos() {
                (f7(omega, d, rs, phi, gamma, psi)?, 7)
            } else {
                (f64::NAN, 99)
            }
        }
    } else {
        if gamma < PI - beta0 && gamma > -PI + beta0 {
            if h * (omega - FRAC_PI_2).tan() <= r * cos_beta0 {
                (f7(omega, d, rs, phi, gamma, psi)?, 7)
            } else {
                (0.0, 6)
            }
        } else if beta0 < gamma || gamma < -beta0 {
            if h * (FRAC_PI_2 - omega).tan() <= r * cos_beta0 {
                (f4(omega, d, rs, gamma)?, 8)
            } else {
                (f5(omega, d, rs, phi, gamma, psi)?, 9)
            }
        } else {
            if h * (omega - FRAC_PI_2).tan() >= r * (beta0 - gamma.abs()).cos() {
                (0.0, 6)
            } else if h * (omega - FRAC_PI_2).tan() <= r * (-beta0 - gamma.abs()).cos() {
                (f5(omega, d, rs, phi, gamma, psi)?, 9)
            } else {
                (f64::NAN, 99)
            }
        }
    };

    Ok((vf, case))
}

fn intersection1(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    z1: f64,
) -> Result<(f64, i32), ViewFactorError> {
    let (vf, case) = if gamma < PI && gamma > 0.0 {
        // Case 10: partial cap, one intersection 1
        (f10(omega, d, rs, phi, gamma, psi, z1)?, 10)
    } else {
        // Case 11: partial cap, one intersection 2
        (f10(omega, d, rs, phi, -gamma, psi, z1)?, 11)
    };

    Ok((vf, case))
}

fn intersection2(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
) -> Result<(f64, i32), ViewFactorError> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    let alpha1 = if cos_phi.abs() < f64::EPSILON {
        (x1.y / (rs * sin_psi)).asin()
    } else {
        (x1.y / (rs * sin_psi)).atan2((x1.x + rs * cos_psi * sin_phi) / (rs * sin_psi * cos_phi))
    };

    let alpha2 = if cos_phi.abs() < f64::EPSILON {
        (x2.y / (rs * sin_psi)).asin()
    } else {
        (x2.y / (rs * sin_psi)).atan2((x2.x + rs * cos_psi * sin_phi) / (rs * sin_psi * cos_phi))
    };

    let alpha_mid = 0.5 * (alpha2 + alpha1);
    let vec = Vector3f {
        x: -rs * cos_psi * sin_phi + rs * alpha_mid.cos() * sin_psi * cos_phi,
        y: rs * alpha_mid.sin() * sin_psi,
        z: d - rs * cos_psi * cos_phi - rs * alpha_mid.cos() * sin_psi * sin_phi,
    };

    let n1 = Vector3f {
        x: omega.sin() * gamma.cos(),
        y: omega.sin() * gamma.sin(),
        z: omega.cos(),
    };

    let (vf, case) = if d * cos_phi >= rs * cos_psi {
        if vec.dot(&n1) > 0.0 {
            if alpha1 > alpha2 {
                (
                    f12(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?,
                    12,
                )
            } else {
                (
                    f12(omega, d, rs, phi, gamma, psi, x2, x1, alpha2, alpha1)?,
                    12,
                )
            }
        } else {
            if alpha1 > alpha2 {
                (
                    f13(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?,
                    13,
                )
            } else {
                (
                    f13(omega, d, rs, phi, gamma, psi, x2, x1, alpha2, alpha1)?,
                    13,
                )
            }
        }
    } else {
        if vec.dot(&n1) > 0.0 {
            if alpha1 < alpha2 {
                (
                    f15(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?,
                    15,
                )
            } else {
                (
                    f15(omega, d, rs, phi, gamma, psi, x2, x1, alpha2, alpha1)?,
                    15,
                )
            }
        } else {
            if alpha1 < alpha2 {
                (
                    f14(omega, d, rs, phi, gamma, psi, x1, x2, alpha1, alpha2)?,
                    14,
                )
            } else {
                (
                    f14(omega, d, rs, phi, gamma, psi, x2, x1, alpha2, alpha1)?,
                    14,
                )
            }
        }
    };

    Ok((vf, case))
}

fn sphericalcap_small(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let cos_gamma = gamma.cos();
    let sin_gamma = gamma.sin();
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    if phi.abs() < f64::EPSILON {
        let theta_ = (rs * sin_psi / (d - rs * cos_psi)).atan();
        if omega >= FRAC_PI_2 + theta_ {
            Ok((0.0, 16)) // Case 16: small cap, view direction mismatch
        } else if 0.0 <= omega && omega <= FRAC_PI_2 - theta_ {
            let vf = f18(omega, d, rs, gamma, psi)?;
            Ok((vf, 18)) // Case 18: small cap, no intersection phi=0
        } else {
            let vf = f19(omega, d, rs, gamma, psi)?;
            Ok((vf, 19)) // Case 19: small cap, two intersections phi=0
        }
    } else if omega.abs() < f64::EPSILON {
        let vf = f17(omega, d, rs, phi, gamma, psi)?;
        Ok((vf, 17)) // Case 17: small cap, no intersection
    } else if gamma.abs() < f64::EPSILON {
        if (cos_phi * sin_omega - sin_phi * cos_omega).abs() < f64::EPSILON {
            small0(omega, d, rs, phi, gamma, psi)
        } else {
            let z = (-rs * cos_psi + d * cos_phi) * sin_omega
                / (cos_phi * sin_omega - cos_omega * sin_phi);
            let x = (rs * cos_psi - d * cos_phi) * cos_omega
                / (cos_phi * sin_omega - cos_omega * sin_phi);
            let y2 = rs.powi(2) - x.powi(2) - (z - d).powi(2);

            if y2 > 0.0 {
                let x1 = Vector3f { x, y: y2.sqrt(), z };
                let x2 = Vector3f {
                    x,
                    y: -y2.sqrt(),
                    z,
                };
                small2(omega, d, rs, phi, gamma, psi, x1, x2)
            } else {
                small0(omega, d, rs, phi, gamma, psi)
            }
        }
    } else if (gamma - PI).abs() < f64::EPSILON || (gamma + PI).abs() < f64::EPSILON {
        if (cos_phi * sin_omega + sin_phi * cos_omega).abs() < f64::EPSILON {
            small0(omega, d, rs, phi, gamma, psi)
        } else {
            let z = (-rs * cos_psi + d * cos_phi) * sin_omega
                / (cos_phi * sin_omega + cos_omega * sin_phi);
            let x = (-rs * cos_psi + d * cos_phi) * cos_omega
                / (cos_phi * sin_omega + cos_omega * sin_phi);
            let y2 = rs.powi(2) - x.powi(2) - (z - d).powi(2);

            if y2 > 0.0 {
                let x1 = Vector3f { x, y: y2.sqrt(), z };
                let x2 = Vector3f {
                    x,
                    y: -y2.sqrt(),
                    z,
                };
                small2(omega, d, rs, phi, gamma, psi, x1, x2)
            } else {
                small0(omega, d, rs, phi, gamma, psi)
            }
        }
    } else {
        let a = (1.0 + cos_gamma.powi(2) / sin_gamma.powi(2)) * cos_phi.powi(2) / sin_phi.powi(2)
            - 2.0 * cos_omega * cos_gamma * cos_phi / (sin_omega * sin_gamma.powi(2) * sin_phi)
            + cos_omega.powi(2) / (sin_omega.powi(2) * sin_gamma.powi(2))
            + 1.0;

        let b = -(1.0 + cos_gamma.powi(2) / sin_gamma.powi(2))
            * (d * cos_phi.powi(2) / sin_phi.powi(2) - rs * cos_phi * cos_psi / sin_phi.powi(2))
            + cos_omega * cos_gamma / (sin_omega * sin_gamma.powi(2))
                * (d * cos_phi / sin_phi - rs * cos_psi / sin_phi)
            - d;

        let c = (1.0 + cos_gamma.powi(2) / sin_gamma.powi(2))
            * (d * cos_phi / sin_phi - rs * cos_psi / sin_phi).powi(2)
            + d.powi(2)
            - rs.powi(2);

        let d_disc = b.powi(2) - a * c;

        if d_disc <= 0.0 {
            small0(omega, d, rs, phi, gamma, psi)
        } else {
            let z1 = (-b + d_disc.sqrt()) / a;
            let z2 = (-b - d_disc.sqrt()) / a;
            let x1 = cos_phi / sin_phi * (d - z1) - cos_psi / sin_phi * rs;
            let x2 = cos_phi / sin_phi * (d - z2) - cos_psi / sin_phi * rs;
            let y1 = -cos_gamma / sin_gamma * x1 - cos_omega / (sin_omega * sin_gamma) * z1;
            let y2 = -cos_gamma / sin_gamma * x2 - cos_omega / (sin_omega * sin_gamma) * z2;

            let w1 = Vector3f {
                x: x2,
                y: y2,
                z: z2,
            };
            let w2 = Vector3f {
                x: x1,
                y: y1,
                z: z1,
            };
            small2(omega, d, rs, phi, gamma, psi, w1, w2)
        }
    }
}

fn small0(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
) -> Result<(f64, i32), ViewFactorError> {
    // Initialize vec vector with x and z components
    let vec = Vector3f {
        x: -rs * psi.cos() * phi.sin(),
        y: 0.0,
        z: d - rs * psi.cos() * phi.cos(),
    };

    // Initialize n1 vector
    let n1 = Vector3f {
        x: omega.sin() * gamma.cos(),
        y: omega.sin() * gamma.sin(),
        z: omega.cos(),
    };

    // Check dot product and determine case
    let (vf, case) = if vec.dot(&n1) <= 0.0 {
        (0.0, 16) // Case 16: small cap, view direction mismatch
    } else {
        (f17(omega, d, rs, phi, gamma, psi)?, 17) // Case 17: small cap, no intersection
    };

    Ok((vf, case))
}

fn small2(
    omega: f64,
    d: f64,
    rs: f64,
    phi: f64,
    gamma: f64,
    psi: f64,
    x1: Vector3f,
    x2: Vector3f,
) -> Result<(f64, i32), ViewFactorError> {
    // Calculate trigonometric values
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let cos_gamma = gamma.cos();
    let sin_gamma = gamma.sin();
    let cos_psi = psi.cos();
    let sin_psi = psi.sin();

    // Calculate theta values
    let sin_theta = rs / d;
    let h = d * (1.0 - sin_theta.powi(2));

    // Initialize n1 vector
    let n1 = Vector3f {
        x: sin_omega * cos_gamma,
        y: sin_omega * sin_gamma,
        z: cos_omega,
    };

    // Calculate alpha values
    let alpha1_temp =
        (x1.y / (rs * sin_psi)).atan2((x1.x + rs * cos_psi * sin_phi) / (rs * sin_psi * cos_phi));
    let alpha2_temp =
        (x2.y / (rs * sin_psi)).atan2((x2.x + rs * cos_psi * sin_phi) / (rs * sin_psi * cos_phi));

    // Calculate mid-point alpha and vector
    let alpha_mid = 0.5 * (alpha2_temp + alpha1_temp);
    let vec = Vector3f {
        x: -rs * cos_psi * sin_phi + rs * alpha_mid.cos() * sin_psi * cos_phi,
        y: rs * alpha_mid.sin() * sin_psi,
        z: d - rs * cos_psi * cos_phi - rs * alpha_mid.cos() * sin_psi * sin_phi,
    };

    let (vf, case) = if n1.dot(&vec) > 0.0 {
        if alpha1_temp > alpha2_temp {
            let x_start = Vector3f {
                x: x2.x * h / x2.z,
                y: x2.y * h / x2.z,
                z: h,
            };
            let x_end = Vector3f {
                x: x1.x * h / x1.z,
                y: x1.y * h / x1.z,
                z: h,
            };
            (
                f20(
                    omega,
                    d,
                    rs,
                    phi,
                    gamma,
                    psi,
                    x_start,
                    x_end,
                    alpha1_temp,
                    alpha2_temp,
                )?,
                20,
            )
        } else {
            let x_start = Vector3f {
                x: x1.x * h / x1.z,
                y: x1.y * h / x1.z,
                z: h,
            };
            let x_end = Vector3f {
                x: x2.x * h / x2.z,
                y: x2.y * h / x2.z,
                z: h,
            };
            (
                f20(
                    omega,
                    d,
                    rs,
                    phi,
                    gamma,
                    psi,
                    x_start,
                    x_end,
                    alpha2_temp,
                    alpha1_temp,
                )?,
                20,
            )
        }
    } else {
        if alpha1_temp > alpha2_temp {
            let x_start = Vector3f {
                x: x1.x * h / x1.z,
                y: x1.y * h / x1.z,
                z: h,
            };
            let x_end = Vector3f {
                x: x2.x * h / x2.z,
                y: x2.y * h / x2.z,
                z: h,
            };
            (
                f21(
                    omega,
                    d,
                    rs,
                    phi,
                    gamma,
                    psi,
                    x_start,
                    x_end,
                    alpha2_temp,
                    alpha1_temp,
                )?,
                21,
            )
        } else {
            let x_start = Vector3f {
                x: x2.x * h / x2.z,
                y: x2.y * h / x2.z,
                z: h,
            };
            let x_end = Vector3f {
                x: x1.x * h / x1.z,
                y: x1.y * h / x1.z,
                z: h,
            };
            (
                f21(
                    omega,
                    d,
                    rs,
                    phi,
                    gamma,
                    psi,
                    x_start,
                    x_end,
                    alpha1_temp,
                    alpha2_temp,
                )?,
                21,
            )
        }
    };

    Ok((vf, case))
}
