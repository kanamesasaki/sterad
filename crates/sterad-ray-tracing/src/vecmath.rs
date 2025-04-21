use crate::error::NumericalError;
use crate::float::Float;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Vector3f {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Vector3f {
    pub fn new(x: Float, y: Float, z: Float) -> Result<Self, NumericalError> {
        if x.is_nan() || y.is_nan() || z.is_nan() {
            return Err(NumericalError::Vector3fErrorNaN);
        }
        if x.is_infinite() || y.is_infinite() || z.is_infinite() {
            return Err(NumericalError::Vector3fErrorInf);
        }
        Ok(Vector3f { x, y, z })
    }

    pub fn dot(&self, other: &Vector3f) -> Float {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(&self, other: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn norm(&self) -> Float {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Vector3f {
        let norm: Float = self.norm();
        Vector3f {
            x: self.x / norm,
            y: self.y / norm,
            z: self.z / norm,
        }
    }

    pub fn has_nan(&self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }
}

// Vector3f + Vector3f = Vector3f
impl Add<Vector3f> for Vector3f {
    type Output = Vector3f;

    fn add(self, rhs: Vector3f) -> Self::Output {
        Vector3f {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

// Vector3f - Vector3f = Vector3f
impl Sub<Vector3f> for Vector3f {
    type Output = Vector3f;

    fn sub(self, rhs: Vector3f) -> Self::Output {
        Vector3f {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

// Vector3f * Float = Vector3f
impl Mul<Float> for Vector3f {
    type Output = Vector3f;

    fn mul(self, rhs: Float) -> Self::Output {
        Vector3f {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

// Float * Vector3f = Vector3f
impl Mul<Vector3f> for Float {
    type Output = Vector3f;

    fn mul(self, rhs: Vector3f) -> Self::Output {
        Vector3f {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Point3f {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Point3f {
    pub fn new(x: Float, y: Float, z: Float) -> Result<Self, NumericalError> {
        if x.is_nan() || y.is_nan() || z.is_nan() {
            return Err(NumericalError::Point3fErrorNaN);
        }
        if x.is_infinite() || y.is_infinite() || z.is_infinite() {
            return Err(NumericalError::Point3fErrorInf);
        }
        Ok(Point3f { x, y, z })
    }
}

// Point3f + Vector3f = Point3f
impl Add<Vector3f> for Point3f {
    type Output = Point3f;

    fn add(self, other: Vector3f) -> Self::Output {
        Point3f {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// Point3f - Vector3f = Point3f
impl Sub<Vector3f> for Point3f {
    type Output = Point3f;

    fn sub(self, other: Vector3f) -> Self::Output {
        Point3f {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// Point3f - Point3f = Vector3f
impl Sub<Point3f> for Point3f {
    type Output = Vector3f;

    fn sub(self, other: Point3f) -> Self::Output {
        Vector3f {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

// Point3f + Point3f = Point3f
impl Add<Point3f> for Point3f {
    type Output = Point3f;

    fn add(self, other: Point3f) -> Self::Output {
        Point3f {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// Point3f * Float = Point3f
impl Mul<Float> for Point3f {
    type Output = Point3f;

    fn mul(self, rhs: Float) -> Self::Output {
        Point3f {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

// Float * Point3f = Point3f
impl Mul<Point3f> for Float {
    type Output = Point3f;

    fn mul(self, rhs: Point3f) -> Self::Output {
        Point3f {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}
