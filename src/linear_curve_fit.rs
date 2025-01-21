/*
A `no_std` and no `alloc` library for linear curve coefficents calculation.
Copyright (C) 2024-2025 joker2770

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
    use anyhow::Result;
    use nalgebra::{ComplexField, SMatrix};
    const POINT_NUM_2D: usize = 8;
    const POINT_NUM_3D: usize = 8;

    const MATRIX_COLUMNS_2D: usize = 2;
    const MATRIX_COLUMNS_3D: usize = 3;

    #[derive(Debug, Default)]
    pub struct LinearCoefficents2D {
        k: f32,
        b: f32,
    }

    impl LinearCoefficents2D {
        /// $$f(x) = kx + b$$
        pub fn value(&self, x: f32) -> f32 {
            self.k * x + self.b
        }

        /// Return function ceofficents `k` and `b`.
        pub fn coefficents(&self) -> (f32, f32) {
            (self.k, self.b)
        }

        pub fn get_matrix_data_from_8_points(
            x: &[f32; POINT_NUM_2D],
            y: &[f32; POINT_NUM_2D],
        ) -> ([f32; POINT_NUM_2D * MATRIX_COLUMNS_2D], [f32; POINT_NUM_2D]) {
            let mut tmp_x = [0.0f32; POINT_NUM_2D * MATRIX_COLUMNS_2D];
            let mut tmp_y = [0.0f32; POINT_NUM_2D];
            for i in 0..POINT_NUM_2D {
                tmp_x[(i * MATRIX_COLUMNS_2D + 0) as usize] = 1.0f32;
                tmp_x[(i * MATRIX_COLUMNS_2D + 1) as usize] = x[i] as f32;
                tmp_y[i] = y[i] as f32;
            }
            (tmp_x, tmp_y)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficients_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_NUM_2D * MATRIX_COLUMNS_2D],
            b: &[f32; POINT_NUM_2D],
            eps: f32,
        ) -> Result<Self> {
            if a.len() == 2 * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_NUM_2D, 1>;
                type MatrixXx2f32 = SMatrix<f32, POINT_NUM_2D, MATRIX_COLUMNS_2D>;
                let mb = MatrixXx1f32::from_row_slice(b);
                let ma = MatrixXx2f32::from_row_slice(a);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, eps);
                match x {
                    Ok(r) => {
                        self.b = r[0];
                        self.k = r[1];
                    }
                    Err(_) => {
                        self.b = 0.0f32;
                        self.k = 0.0f32;
                        return Err(anyhow::anyhow!("SVD solve failed"));
                    }
                }
            } else {
                self.b = 0.0f32;
                self.k = 0.0f32;
                return Err(anyhow::anyhow!("Matrix size not match"));
            }

            Ok(Self {
                b: self.b,
                k: self.k,
            })
        }
    }

    #[derive(Debug, Default)]
    pub struct LinearCoefficents3D {
        a: f32,
        b: f32,
        c: f32,
    }

    impl LinearCoefficents3D {
        /// $$f(x, y) = ax + by + c$$
        pub fn value(&self, x: f32, y: f32) -> f32 {
            self.a * x + self.b * y + self.c
        }

        /// Return function ceofficents `a`, `b` and `c`.
        pub fn coefficents(&self) -> (f32, f32, f32) {
            (self.a, self.b, self.c)
        }

        pub fn get_matrix_data_from_8_points(
            x: &[f32; POINT_NUM_3D],
            y: &[f32; POINT_NUM_3D],
            z: &[f32; POINT_NUM_3D],
        ) -> ([f32; POINT_NUM_3D * MATRIX_COLUMNS_3D], [f32; POINT_NUM_3D]) {
            let mut tmp_x_y = [0.0f32; POINT_NUM_3D * MATRIX_COLUMNS_3D];
            let mut tmp_z = [0.0f32; POINT_NUM_3D];
            for i in 0..POINT_NUM_3D {
                tmp_x_y[(i * MATRIX_COLUMNS_3D + 0) as usize] = 1.0f32;
                tmp_x_y[(i * MATRIX_COLUMNS_3D + 1) as usize] = x[i] as f32;
                tmp_x_y[(i * MATRIX_COLUMNS_3D + 2) as usize] = y[i] as f32;
                tmp_z[i] = z[i] as f32;
            }
            (tmp_x_y, tmp_z)
        }

        /// Solve ax = b.
        /// Any singular value smaller than `eps` is assumed to be zero.
        pub fn get_coefficients_from_8_matrix_data(
            &mut self,
            a: &[f32; POINT_NUM_3D * MATRIX_COLUMNS_3D],
            b: &[f32; POINT_NUM_3D],
            eps: f32,
        ) -> Result<Self> {
            if a.len() == MATRIX_COLUMNS_3D * b.len() {
                type MatrixXx1f32 = SMatrix<f32, POINT_NUM_3D, 1>;
                type MatrixXx3f32 = SMatrix<f32, POINT_NUM_3D, MATRIX_COLUMNS_3D>;
                let mb = MatrixXx1f32::from_row_slice(b);
                let ma = MatrixXx3f32::from_row_slice(a);

                let decomp = ma.svd(true, true);

                let x = decomp.solve(&mb, eps);
                match x {
                    Ok(r) => {
                        self.c = r[0];
                        self.a = r[1];
                        self.b = r[2];
                    }
                    Err(_) => {
                        self.c = 0.0f32;
                        self.a = 0.0f32;
                        self.b = 0.0f32;
                        return Err(anyhow::anyhow!("SVD solve failed"));
                    }
                }
            } else {
                self.c = 0.0f32;
                self.a = 0.0f32;
                self.b = 0.0f32;
                return Err(anyhow::anyhow!("Matrix size not match"));
            }

            Ok(Self {
                a: self.a,
                b: self.b,
                c: self.c,
            })
        }
    }
}
