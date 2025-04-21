use crate::float::Float;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Ray {
    pub o: Point3f,
    pub d: Vector3f,
    pub time: f64,
}

impl Ray {
    pub fn new(o: Point3f, d: Vector3f, time: f64) -> Self {
        Ray { o, d, time }
    }

    pub fn operator(&self, t: Float) -> Point3f {
        self.o + (self.d * t)
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct RayDifferential {
    pub ray: Ray,
    pub has_differentials: bool,
    pub rx_origin: Point3f,
    pub ry_origin: Point3f,
    pub rx_direction: Vector3f,
    pub ry_direction: Vector3f,
}

impl RayDifferential {
    pub fn new(o: Point3f, d: Vector3f, time: f64) -> Self {
        RayDifferential {
            ray: Ray::new(o, d, time),
            has_differentials: false,
            rx_origin: Point3f::default(),
            ry_origin: Point3f::default(),
            rx_direction: Vector3f::default(),
            ry_direction: Vector3f::default(),
        }
    }

    pub fn scale_differentials(&mut self, s: Float) {
        self.rx_origin = self.ray.o + ((self.rx_origin - self.ray.o) * s);
        self.ry_origin = self.ray.o + ((self.ry_origin - self.ray.o) * s);
        self.rx_direction = self.ray.d + ((self.rx_direction - self.ray.d) * s);
        self.ry_direction = self.ray.d + ((self.ry_direction - self.ray.d) * s);
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct QuadricIntersection {
    pub t_hit: Float,
    pub p_obj: Point3f,
    pub phi: Float,
}
