use crate::error::NumericalError;
use crate::float::{float::PI, Float};
use crate::ray::{QuadricIntersection, Ray};
use crate::transform::Transform;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Disk {
    pub max_radius: Float,
    pub min_radius: Float,
    pub start_angle: Float,
    pub end_angle: Float,
    pub trans: Transform,
}

impl Default for Disk {
    fn default() -> Self {
        Disk {
            max_radius: 1.0,
            min_radius: 0.0,
            start_angle: 0.0,
            end_angle: 360.0,
            trans: Transform::default(),
        }
    }
}

impl Disk {
    fn random_point(&self, rand1: Float, rand2: Float) -> Result<Point3f, NumericalError> {
        let r: Float =
            rand1 * (self.max_radius.powi(2) - self.min_radius.powi(2)) + self.min_radius.powi(2);
        let theta: Float = rand2 * (self.end_angle - self.start_angle) + self.start_angle;
        let p_local: Point3f = Point3f {
            x: r * theta.cos(),
            y: r * theta.sin(),
            z: 0.0,
        };
        let p: Point3f = self.trans.render_from_object_point(p_local)?;
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
        let n: Vector3f = self.trans.render_from_object_vector(n_local)?;
        Ok(n)
    }

    pub fn spawn_ray(
        &self,
        rand1: Float,
        rand2: Float,
        rand3: Float,
        rand4: Float,
    ) -> Result<Ray, NumericalError> {
        let p: Point3f = self.random_point(rand1, rand2)?;
        let n: Vector3f = self.random_vector(rand3, rand4)?;
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
        let dist2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if (dist2 > self.max_radius.powi(2)) || (dist2 < self.min_radius.powi(2)) {
            return Ok(None);
        }
        let phi_atan: Float = Float::atan2(p_hit.y, p_hit.x);
        let phi: Float = if phi_atan < 0.0 {
            phi_atan + 2.0 * PI
        } else {
            phi_atan
        };

        if (phi < self.start_angle) || (phi > self.end_angle) {
            return Ok(None);
        }

        let intersection: QuadricIntersection = QuadricIntersection {
            t_hit: t_shape_hit,
            p_obj: p_hit,
            phi,
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
