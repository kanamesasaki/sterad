use crate::error::NumericalError;
use crate::float::Float;
use crate::vecmath::{Point3f, Vector3f};

/// Represents a 4x4 homogeneous transformation matrix for 3D coordinate transformations.
/// <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
///
/// $$
/// \\begin{gather}
/// \boldsymbol{m} = \left[ \begin{array}{cccc}
/// m[0][0] & m[0][1] & m[0][2] & m[0][3] \\\\
/// m[1][0] & m[1][1] & m[1][2] & m[1][3] \\\\
/// m[2][0] & m[2][1] & m[2][2] & m[2][3] \\\\
/// m[3][0] & m[3][1] & m[3][2] & m[3][3] \\\\
/// \\end{array} \right]
/// \\end{gather}
/// $$
///
/// The following column vectors represent the base vectors of object x, y, and z axes.
/// And the last column vector represents the origin of the object.
/// $$
/// \begin{gather}
/// \boldsymbol{e}_x = \left[ \begin{array}{c} m[0][0] \\\\ m[1][0] \\\\ m[2][0] \end{array} \right], \quad
/// \boldsymbol{e}_y = \left[ \begin{array}{c} m[0][1] \\\\ m[1][1] \\\\ m[2][1] \end{array} \right], \quad
/// \boldsymbol{e}_z = \left[ \begin{array}{c} m[0][2] \\\\ m[1][2] \\\\ m[2][2] \end{array} \right], \quad
/// \boldsymbol{o} = \left[ \begin{array}{c} m[0][3] \\\\ m[1][3] \\\\ m[2][3] \end{array} \right].
/// \end{gather}
/// $$
///
/// The last row of the matrix is fixed to the following values.
/// $$
/// \begin{gather}
/// m[3][0] = 0, \quad m[3][1] = 0, \quad m[3][2] = 0, \quad m[3][3] = 1.
/// \end{gather}
/// $$
///
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
    // translate the object origin to (x, y, z)
    pub fn translate(&mut self, x: Float, y: Float, z: Float) {
        self.m[0][3] = x;
        self.m[1][3] = y;
        self.m[2][3] = z;
    }

    /// Rotate base vectors around local z axis
    ///
    /// $$
    /// \begin{gather}
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right] =
    /// \left[ \begin{array}{ccc}
    /// \cos\theta & \sin\theta & 0 \\\\
    /// -\sin\theta & \cos\theta & 0 \\\\
    /// 0 & 0 & 1
    /// \end{array} \right]
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right]
    /// \end{gather}
    /// $$
    ///
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

    /// Rotate base vectors around local y axis
    ///
    /// $$
    /// \begin{gather}
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right] =
    /// \left[ \begin{array}{ccc}
    /// \cos\theta & 0 & -\sin\theta \\\\
    /// 0 & 1 & 0 \\\\
    /// \sin\theta & 0 & \cos\theta
    /// \end{array} \right]
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right]
    /// \end{gather}
    /// $$
    ///
    pub fn rot_y(&mut self, angle: Float) {
        let cos_angle: Float = angle.cos();
        let sin_angle: Float = angle.sin();
        self.m[0][0] = cos_angle * self.m[0][0] - sin_angle * self.m[0][2];
        self.m[1][0] = cos_angle * self.m[1][0] - sin_angle * self.m[1][2];
        self.m[2][0] = cos_angle * self.m[2][0] - sin_angle * self.m[2][2];
        // self.m[0][1] = self.m[0][1];
        // self.m[1][1] = self.m[1][1];
        // self.m[2][1] = self.m[2][1];
        self.m[0][2] = sin_angle * self.m[0][0] + cos_angle * self.m[0][2];
        self.m[1][2] = sin_angle * self.m[1][0] + cos_angle * self.m[1][2];
        self.m[2][2] = sin_angle * self.m[2][0] + cos_angle * self.m[2][2];
    }

    /// Rotate base vectors around local x axis
    ///
    /// $$
    /// \begin{gather}
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right] =
    /// \left[ \begin{array}{ccc}
    /// 1 & 0 & 0 \\\\
    /// 0 & \cos\theta & \sin\theta \\\\
    /// 0 & -\sin\theta & \cos\theta
    /// \end{array} \right]
    /// \left[ \begin{array}{c} \boldsymbol{e}^\mathrm{T}_x \\\\ \boldsymbol{e}^\mathrm{T}_y \\\\ \boldsymbol{e}^\mathrm{T}_z \end{array} \right]
    /// \end{gather}
    /// $$
    ///
    pub fn rot_x(&mut self, angle: Float) {
        let cos_angle: Float = angle.cos();
        let sin_angle: Float = angle.sin();
        // self.m[0][0] = self.m[0][0];
        // self.m[1][0] = self.m[1][0];
        // self.m[2][0] = self.m[2][0];
        self.m[0][1] = cos_angle * self.m[0][1] + sin_angle * self.m[0][2];
        self.m[1][1] = cos_angle * self.m[1][1] + sin_angle * self.m[1][2];
        self.m[2][1] = cos_angle * self.m[2][1] + sin_angle * self.m[2][2];
        self.m[0][2] = -sin_angle * self.m[0][1] + cos_angle * self.m[0][2];
        self.m[1][2] = -sin_angle * self.m[1][1] + cos_angle * self.m[1][2];
        self.m[2][2] = -sin_angle * self.m[2][1] + cos_angle * self.m[2][2];
    }

    // move a point from object coordinate to render coordinate
    pub fn render_from_object_point(&self, p: Point3f) -> Result<Point3f, NumericalError> {
        let x: Float = self.m[0][0] * p.x + self.m[0][1] * p.y + self.m[0][2] * p.z + self.m[0][3];
        let y: Float = self.m[1][0] * p.x + self.m[1][1] * p.y + self.m[1][2] * p.z + self.m[1][3];
        let z: Float = self.m[2][0] * p.x + self.m[2][1] * p.y + self.m[2][2] * p.z + self.m[2][3];
        Point3f::new(x, y, z)
    }

    // move a vector from object coordinate to render coordinate
    pub fn render_from_object_vector(&self, v: Vector3f) -> Result<Vector3f, NumericalError> {
        let x: Float = self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z;
        let y: Float = self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z;
        let z: Float = self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z;
        Vector3f::new(x, y, z)
    }

    // move a point from render coordinate to object coordinate
    pub fn object_from_render_vector(&self, v: Vector3f) -> Result<Vector3f, NumericalError> {
        let x: Float = self.m[0][0] * v.x + self.m[1][0] * v.y + self.m[2][0] * v.z;
        let y: Float = self.m[0][1] * v.x + self.m[1][1] * v.y + self.m[2][1] * v.z;
        let z: Float = self.m[0][2] * v.x + self.m[1][2] * v.y + self.m[2][2] * v.z;
        Vector3f::new(x, y, z)
    }

    // move a point from render coordinate to object coordinate
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
