use std::fmt;

#[derive(Debug)]
pub enum SolidAngleError {
    InvalidInput {
        param_name: &'static str,
        message: String,
    },
    RuntimeError {
        message: String,
    },
}

impl fmt::Display for SolidAngleError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SolidAngleError::InvalidInput {
                param_name,
                message,
            } => {
                write!(
                    f,
                    "Invalid input for parameter '{}': {}",
                    param_name, message
                )
            }
            SolidAngleError::RuntimeError { message } => {
                write!(f, "Runtime error: {}", message)
            }
        }
    }
}
