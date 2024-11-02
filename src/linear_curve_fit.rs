/*
A no_std and no heap memory library for linear curve coefficents calculation.
Copyright (C) 2024  joker2770

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#![allow(unused_imports)]

pub mod linear_curve {

    use nalgebra::{ComplexField, SMatrix};
    const POINT_CNT: usize = 8;

    #[derive(Debug, Default)]
    pub struct LinearCoefficents2D {
        k: f32,
        b: f32,
    }

    impl LinearCoefficents2D {
        pub fn value(&self, x: f32) -> f32 {
            self.k * x + self.b
        }

        pub fn coefficents(&self) -> (f32, f32) {
            (self.k, self.b)
        }

        pub fn get_matrix_data_from_8_points(
            x: &[f32; POINT_CNT],
            y: &[f32; POINT_CNT],
        ) -> ([f32; POINT_CNT * 2], [f32; POINT_CNT]) {
            let mut tmp_x = [0.0f32; POINT_CNT * 2];
            let mut tmp_y = [0.0f32; POINT_CNT];
            for i in 0..POINT_CNT {
                tmp_x[(i * 2 + 0) as usize] = 1.0f32;
                tmp_x[(i * 2 + 1) as usize] = x[i] as f32;
                tmp_y[i] = y[i] as f32;
            }
            (tmp_x, tmp_y)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficients_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_CNT * 2],
            b: &[f32; POINT_CNT],
            eps: f32,
        ) {
            if a.len() == 2 * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_CNT, 1>;
                type MatrixXx2f32 = SMatrix<f32, POINT_CNT, 2>;
                let mb = MatrixXx1f32::from_row_slice(b);
                let ma = MatrixXx2f32::from_row_slice(a);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, eps);
                if let Ok(r) = x {
                    self.b = r[0];
                    self.k = r[1];
                }
            }
        }
    }
}
