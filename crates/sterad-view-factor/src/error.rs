use std::fmt;

#[derive(Debug)]
pub enum ViewFactorError {
    InvalidInput {
        param_name: &'static str,
        message: String,
    },
    RuntimeError {
        message: String,
    },
}

impl fmt::Display for ViewFactorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ViewFactorError::InvalidInput {
                param_name,
                message,
            } => {
                write!(
                    f,
                    "Invalid input for parameter '{}': {}",
                    param_name, message
                )
            }
            ViewFactorError::RuntimeError { message } => {
                write!(f, "Runtime error: {}", message)
            }
        }
    }
}
