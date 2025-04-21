use std::fmt;

#[derive(Debug)]
pub enum NumericalError {
    Vector3fErrorNaN,
    Vector3fErrorInf,
    Point3fErrorNaN,
    Point3fErrorInf,
    TriangleDegenerate,
}

impl fmt::Display for NumericalError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NumericalError::Vector3fErrorNaN => {
                write!(f, "NaN value encountered at Vector3f initialization")
            }
            NumericalError::Vector3fErrorInf => {
                write!(f, "Infinity value encountered at Vector3f initialization")
            }
            NumericalError::Point3fErrorNaN => {
                write!(f, "NaN value encountered at Point3f initialization")
            }
            NumericalError::Point3fErrorInf => {
                write!(f, "Infinity value encountered at Point3f initialization")
            }
            NumericalError::TriangleDegenerate => {
                write!(f, "Triangle is degenerate")
            }
        }
    }
}
