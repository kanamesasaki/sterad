use crate::float::{float, Float};
use std::ops::{Add, Sub};

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Interval {
    pub low: Float,
    pub high: Float,
}

impl Interval {
    fn new(low: Float, high: Float) -> Self {
        Self {
            low: low.min(high),
            high: low.max(high),
        }
    }

    fn lower_bound(&self) -> Float {
        self.low
    }

    fn upper_bound(&self) -> Float {
        self.high
    }

    fn mid_point(&self) -> Float {
        (self.low + self.high) / 2.0
    }

    fn width(&self) -> Float {
        self.high - self.low
    }
}

impl Add<&Interval> for &Interval {
    type Output = Interval;

    fn add(self, other: &Interval) -> Self::Output {
        Interval {
            low: float::add_round_down(self.low, other.low),
            high: float::add_round_up(self.high, other.high),
        }
    }
}

impl Add<Float> for &Interval {
    type Output = Interval;

    fn add(self, other: Float) -> Self::Output {
        Interval {
            low: self.low + other,
            high: self.high + other,
        }
    }
}

impl Add<&Interval> for Float {
    type Output = Interval;

    fn add(self, other: &Interval) -> Self::Output {
        Interval {
            low: self + other.low,
            high: self + other.high,
        }
    }
}
