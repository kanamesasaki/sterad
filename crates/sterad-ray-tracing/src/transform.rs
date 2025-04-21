use crate::error::NumericalError;
use crate::float::Float;
use crate::vecmath::{Point3f, Vector3f};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Transform {
    pub m: [[Float; 4]; 4],
}

impl Default for Transform {
    fn default() -> Self {
        Transform {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }
}

impl Transform {
    pub fn translate(&mut self, x: Float, y: Float, z: Float) {
        self.m[0][3] = x;
        self.m[1][3] = y;
        self.m[2][3] = z;
    }

    // rotate base vectors around local x axis
    pub fn rot_z(&mut self, angle: Float) {
        let cos_angle: Float = angle.cos();
        let sin_angle: Float = angle.sin();
        self.m[0][0] = cos_angle * self.m[0][0] + sin_angle * self.m[0][1];
        self.m[1][0] = cos_angle * self.m[1][0] + sin_angle * self.m[1][1];
        self.m[2][0] = cos_angle * self.m[2][0] + sin_angle * self.m[2][1];
        self.m[0][1] = -sin_angle * self.m[0][0] + cos_angle * self.m[0][1];
        self.m[1][1] = -sin_angle * self.m[1][0] + cos_angle * self.m[1][1];
        self.m[2][1] = -sin_angle * self.m[2][0] + cos_angle * self.m[2][1];
        // self.m[0][2] = self.m[0][2];
        // self.m[1][2] = self.m[1][2];
        // self.m[2][2] = self.m[2][2];
    }

    pub fn render_from_object_point(&self, p: Point3f) -> Result<Point3f, NumericalError> {
        let x: Float = self.m[0][0] * p.x + self.m[0][1] * p.y + self.m[0][2] * p.z + self.m[0][3];
        let y: Float = self.m[1][0] * p.x + self.m[1][1] * p.y + self.m[1][2] * p.z + self.m[1][3];
        let z: Float = self.m[2][0] * p.x + self.m[2][1] * p.y + self.m[2][2] * p.z + self.m[2][3];
        Point3f::new(x, y, z)
    }
    pub fn render_from_object_vector(&self, v: Vector3f) -> Result<Vector3f, NumericalError> {
        let x: Float = self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z;
        let y: Float = self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z;
        let z: Float = self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z;
        Vector3f::new(x, y, z)
    }

    pub fn object_from_render_vector(&self, v: Vector3f) -> Result<Vector3f, NumericalError> {
        let x: Float = self.m[0][0] * v.x + self.m[1][0] * v.y + self.m[2][0] * v.z;
        let y: Float = self.m[0][1] * v.x + self.m[1][1] * v.y + self.m[2][1] * v.z;
        let z: Float = self.m[0][2] * v.x + self.m[1][2] * v.y + self.m[2][2] * v.z;
        Vector3f::new(x, y, z)
    }

    pub fn object_from_render_point(&self, p: Point3f) -> Result<Point3f, NumericalError> {
        let x_rel: Float = p.x - self.m[0][3];
        let y_rel: Float = p.y - self.m[1][3];
        let z_rel: Float = p.z - self.m[2][3];
        let x: Float = self.m[0][0] * x_rel + self.m[1][0] * y_rel + self.m[2][0] * z_rel;
        let y: Float = self.m[0][1] * x_rel + self.m[1][1] * y_rel + self.m[2][1] * z_rel;
        let z: Float = self.m[0][2] * x_rel + self.m[1][2] * y_rel + self.m[2][2] * z_rel;
        Point3f::new(x, y, z)
    }
}
