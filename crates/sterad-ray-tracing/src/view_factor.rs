use crate::disk::Disk;
use crate::error::NumericalError;
use crate::float::Float;
use crate::plate_element::PlateElement;
use crate::ray::Ray;
use rand::Rng;

pub struct ViewFactorMonteCarlo {
    pub num_samples: usize,
    pub max_distance: Float,
}

impl Default for ViewFactorMonteCarlo {
    fn default() -> Self {
        ViewFactorMonteCarlo {
            num_samples: 100000,
            max_distance: Float::INFINITY,
        }
    }
}

impl ViewFactorMonteCarlo {
    /// モンテカルロ法でPlateElementからDiskへのView Factorを計算
    pub fn plate_element_to_disk(
        &self,
        plate: &PlateElement,
        disk: &Disk,
    ) -> Result<Float, NumericalError> {
        let mut rng = rand::rng();
        let mut hits = 0;

        for _ in 0..self.num_samples {
            // プレートからランダムな方向にレイを生成
            let rand1: Float = rng.random();
            let rand2: Float = rng.random();

            let ray: Ray = plate.spawn_ray(rand1, rand2)?;

            // レイがディスクと交差するか判定
            if disk.intersect_p(&ray, self.max_distance)? {
                hits += 1;
            }
        }

        // View Factorを計算 (ヒット数 / サンプル数)
        let view_factor = hits as Float / self.num_samples as Float;

        Ok(view_factor)
    }
}

// 既存のコードの末尾に追加

#[cfg(test)]
mod tests {
    use super::*;
    use crate::float::float::PI;
    use crate::transform::Transform;

    #[test]
    fn test_plate_element_to_disk() -> Result<(), NumericalError> {
        // 原点にあるPlateElement
        let plate = PlateElement {
            trans: Transform::default(),
        };

        // X方向に距離1だけ離れた位置にあるディスク（Y-Z平面に垂直）
        let mut disk_trans = Transform::default();
        disk_trans.translate(1.0, 0.0, 1.0);

        let disk = Disk {
            max_radius: 2.0,
            min_radius: 0.0,
            start_angle: 0.0,
            end_angle: 2.0 * PI,
            trans: disk_trans,
        };

        let vf_calc = ViewFactorMonteCarlo {
            num_samples: 100000,
            max_distance: 10.0,
        };

        let vf: Float = vf_calc.plate_element_to_disk(&plate, &disk)?;
        let expected: Float = 0.7236067977499789;
        let tolerance: Float = 0.0001;
        assert!(
            (vf - expected).abs() < tolerance,
            "vf calculated: {}, vf expected: {}",
            vf,
            expected
        );
        Ok(())
    }
}
