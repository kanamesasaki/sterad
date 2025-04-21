use crate::error::NumericalError;
use crate::float::{float::PI, Float};
use crate::ray::{QuadricIntersection, Ray};
use crate::transform::Transform;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Rectangle {
    pub x_max: Float,
    pub y_max: Float,
    pub trans: Transform,
}

impl Default for Rectangle {
    fn default() -> Self {
        Rectangle {
            x_max: 1.0,
            y_max: 1.0,
            trans: Transform::default(),
        }
    }
}

impl Rectangle {
    fn random_point(&self, rand1: Float, rand2: Float) -> Result<Point3f, NumericalError> {
        let p_local: Point3f = Point3f {
            x: self.x_max * rand1,
            y: self.y_max * rand2,
            z: 0.0,
        };
        let p: Point3f = self.trans.object_from_render_point(p_local)?;
        Ok(p)
    }

    fn random_vector(&self, rand1: Float, rand2: Float) -> Result<Vector3f, NumericalError> {
        let phi: Float = 2.0 * PI * rand1;
        let theta: Float = Float::acos(1.0 - 2.0 * rand2) / 2.0;
        let n_local: Vector3f = Vector3f {
            x: theta.sin() * phi.cos(),
            y: theta.sin() * phi.sin(),
            z: theta.cos(),
        };
        let n: Vector3f = self.trans.object_from_render_vector(n_local)?;
        Ok(n)
    }

    pub fn spawn_ray(&self, rand1: Float, rand2: Float) -> Result<Ray, NumericalError> {
        let p: Point3f = self.random_point(rand1, rand2)?;
        let n: Vector3f = self.random_vector(rand1, rand2)?;
        let ray: Ray = Ray {
            o: p,
            d: n,
            time: 0.0,
        };
        Ok(ray)
    }

    fn basic_intersect(
        &self,
        ray: &Ray,
        t_max: Float,
    ) -> Result<Option<QuadricIntersection>, NumericalError> {
        let oi: Point3f = self.trans.object_from_render_point(ray.o)?;
        let di: Vector3f = self.trans.object_from_render_vector(ray.d)?;

        // Reject disk intersections for rays parallel to the disk's plane
        if di.z == 0.0 {
            return Ok(None);
        }

        let t_shape_hit: Float = -oi.z / di.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= t_max {
            return Ok(None);
        }

        let p_hit: Point3f = oi + (t_shape_hit * di);
        if (p_hit.x > self.x_max) || (p_hit.x < 0.0) {
            return Ok(None);
        }
        if (p_hit.y > self.y_max) || (p_hit.y < 0.0) {
            return Ok(None);
        }
        let intersection: QuadricIntersection = QuadricIntersection {
            t_hit: t_shape_hit,
            p_obj: p_hit,
            phi: 0.0,
        };
        Ok(Some(intersection))
    }

    pub fn intersect_p(&self, ray: &Ray, t_max: Float) -> Result<bool, NumericalError> {
        let intersection: Option<QuadricIntersection> = self.basic_intersect(ray, t_max)?;
        if intersection.is_some() {
            return Ok(true);
        }
        Ok(false)
    }
}
