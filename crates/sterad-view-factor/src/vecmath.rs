#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3f {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3f {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vector3f { x, y, z }
    }

    // dot product
    pub fn dot(self, other: &Vector3f) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    // cross product
    pub fn cross(self, other: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    // norm
    pub fn norm(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    // normalize
    pub fn normalize(self) -> Vector3f {
        let norm = self.norm();
        Vector3f {
            x: self.x / norm,
            y: self.y / norm,
            z: self.z / norm,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix3f {
    pub x_axis: Vector3f,
    pub y_axis: Vector3f,
    pub z_axis: Vector3f,
}

impl Matrix3f {
    pub fn from_columns(col_x: &Vector3f, col_y: &Vector3f, col_z: &Vector3f) -> Self {
        Matrix3f {
            x_axis: *col_x,
            y_axis: *col_y,
            z_axis: *col_z,
        }
    }

    pub fn from_rows(row_x: &Vector3f, row_y: &Vector3f, row_z: &Vector3f) -> Self {
        Matrix3f {
            x_axis: Vector3f::new(row_x.x, row_y.x, row_z.x),
            y_axis: Vector3f::new(row_x.y, row_y.y, row_z.y),
            z_axis: Vector3f::new(row_x.z, row_y.z, row_z.z),
        }
    }

    pub fn mul_vector(self, v: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.x_axis.x * v.x + self.y_axis.x * v.y + self.z_axis.x * v.z,
            y: self.x_axis.y * v.x + self.y_axis.y * v.y + self.z_axis.y * v.z,
            z: self.x_axis.z * v.x + self.y_axis.z * v.y + self.z_axis.z * v.z,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Vector3f;

    #[test]
    fn test_dot_product() {
        let v1 = Vector3f {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        let v2 = Vector3f {
            x: 4.0,
            y: 5.0,
            z: 6.0,
        };
        let result = v1.dot(&v2);
        assert_eq!(result, 32.0);
    }

    #[test]
    fn test_cross_product() {
        let v1 = Vector3f {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        };
        let v2 = Vector3f {
            x: 4.0,
            y: 5.0,
            z: 6.0,
        };
        let result = v1.cross(&v2);
        assert_eq!(
            result,
            Vector3f {
                x: -3.0,
                y: 6.0,
                z: -3.0
            }
        );
    }

    #[test]
    fn test_norm() {
        let v = Vector3f {
            x: 1.0,
            y: 2.0,
            z: 2.0,
        };
        let result = v.norm();
        assert_eq!(result, 3.0);
    }
}
