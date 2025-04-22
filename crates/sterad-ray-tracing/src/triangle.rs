use crate::error::NumericalError;
use crate::float::Float;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Triangle {
    pub p0: Point3f,
    pub p1: Point3f,
    pub p2: Point3f,
}

impl Default for Triangle {
    fn default() -> Self {
        Triangle {
            p0: Point3f::new(0.0, 0.0, 0.0).unwrap(),
            p1: Point3f::new(1.0, 0.0, 0.0).unwrap(),
            p2: Point3f::new(0.0, 1.0, 0.0).unwrap(),
        }
    }
}

impl Triangle {
    pub fn new(p0: Point3f, p1: Point3f, p2: Point3f) -> Result<Self, NumericalError> {
        let u: Vector3f = p1 - p0;
        let v: Vector3f = p2 - p0;
        let n: Vector3f = u.cross(&v);
        if n.dot(&n) == 0.0 {
            return Err(NumericalError::TriangleDegenerate);
        }
        Ok(Triangle { p0, p1, p2 })
    }

    // Generate a random point inside the triangle
    fn sample(&self, rand0: Float, rand1: Float) -> Point3f {
        let b0: Float;
        let b1: Float;
        if rand0 < rand1 {
            b0 = rand0;
            b1 = rand1;
        } else {
            b0 = rand1;
            b1 = rand0;
        }
        let b2: Float = 1.0 - b0 - b1;
        let p: Point3f = self.p0 * b0 + self.p1 * b1 + self.p2 * b2;
        p
    }
}
