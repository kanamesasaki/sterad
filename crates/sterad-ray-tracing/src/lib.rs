pub mod cylinder;
pub mod disk;
pub mod error;
pub mod float;
pub mod interval;
pub mod plate_element;
pub mod ray;
pub mod rectangle;
pub mod transform;
pub mod triangle;
pub mod vecmath;
pub mod view_factor;

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
