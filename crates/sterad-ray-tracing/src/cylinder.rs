use crate::error::NumericalError;
use crate::float::{float::PI, Float};
use crate::ray::{QuadricIntersection, Ray};
use crate::transform::Transform;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Cylinder {
    pub radius: Float,
    pub height: Float,
    pub start_angle: Float,
    pub end_angle: Float,
    pub trans: Transform,
}

impl Default for Cylinder {
    fn default() -> Self {
        Cylinder {
            radius: 1.0,
            height: 1.0,
            start_angle: 0.0,
            end_angle: 2.0 * PI,
            trans: Transform::default(),
        }
    }
}

impl Cyliner {
    pub fn new(
        radius: Float,
        height: Float,
        start_angle: Float,
        end_angle: Float,
        trans: Transform,
    ) -> Self {
        Cylinder {
            radius,
            height,
            start_angle,
            end_angle,
            trans,
        }
    }

    pub fn sample(&self, rand0: Float, rand1: Float) -> Result<Point3f, NumericalError> {
        let phi: Float = 2.0 * PI * rand0;
        let z: Float = self.height * rand1;
        let x: Float = self.radius * phi.cos();
        let y: Float = self.radius * phi.sin();
        let p_local: Point3f = Point3f { x, y, z };
        let p: Point3f = self.trans.object_from_render_point(p_local)?;
        Ok(p)
    }

    fn diffuse(
        &self,
        alpha: Float,
        rand0: Float,
        rand1: Float,
    ) -> Result<Vector3f, NumericalError> {
        // x: (-sin(alpha), cos(alpha), 0)
        // y: (0, 0, 1)
        // z: (cos(alpha), sin(alpha), 0)
        // v = (x, y, z)

        let phi: Float = (self.end_angle - self.start_angle) * rand0 + self.start_angle;
        let theta: Float = Float::acos(1.0 - 2.0 * rand1) / 2.0;
        let v_local: Vector3f = Vector3f {
            x: -alpha.sin() * theta.sin() * phi.cos() + alpha.cos() * theta.cos(),
            y: alpha.cos() * theta.sin() * phi.cos() + alpha.sin() * theta.cos(),
            z: theta.sin() * phi.sin(),
        };

        let v: Vector3f = self.trans.render_from_object_vector(v_local)?;
        Ok(v)
    }

    fn basic_intersect(
        &self,
        ray: &Ray,
        t_max: Float,
    ) -> Result<Option<QuadricIntersection>, NumericalError> {
        let oi: Point3f = self.trans.object_from_render_point(ray.o)?;
        let di: Vector3f = self.trans.object_from_render_vector(ray.d)?;

        let a: Float = di.x.powi(2) + di.y.powi(2);
        let b: Float = 2.0 * (di.x * oi.x + di.y * oi.y);
        
    }
}
