use crate::error::NumericalError;
use crate::float::{float::PI, Float};
use crate::ray::Ray;
use crate::transform::Transform;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Default, Clone, PartialEq)]
pub struct PlateElement {
    pub trans: Transform,
}

impl PlateElement {
    fn random_vector(&self, rand1: Float, rand2: Float) -> Result<Vector3f, NumericalError> {
        let phi: Float = 2.0 * PI * rand1;
        let theta: Float = Float::acos(1.0 - 2.0 * rand2) / 2.0;
        let n_local: Vector3f = Vector3f {
            x: theta.sin() * phi.cos(),
            y: theta.sin() * phi.sin(),
            z: theta.cos(),
        };
        let n: Vector3f = self.trans.render_from_object_vector(n_local)?;
        Ok(n)
    }

    pub fn spawn_ray(&self, rand1: Float, rand2: Float) -> Result<Ray, NumericalError> {
        let n: Vector3f = self.random_vector(rand1, rand2)?;
        let p: Point3f = Point3f {
            x: (self.trans.m[0][3]),
            y: (self.trans.m[1][3]),
            z: (self.trans.m[2][3]),
        };
        let ray: Ray = Ray {
            o: p,
            d: n,
            time: 0.0,
        };
        Ok(ray)
    }
}
