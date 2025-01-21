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

#![no_std]

pub mod linear_curve_fit;

pub use linear_curve_fit::linear_curve;

#[cfg(test)]
mod tests {
    use linear_curve::LinearCoefficents2D;

    use super::*;

    #[test]
    fn it_works() {
        let mut linear_coefficents = LinearCoefficents2D::default();
        let x_test_arr = [
            -2.8f32, -1.6f32, -0.5f32, 5.0f32, 5.4f32, 6.7f32, 10.3f32, 13.8f32,
        ];
        let y_test_arr = [
            33.1f32, 21.1f32, 9.9f32, -45.2f32, -49.1f32, -61.9f32, -98.1f32, -132.99f32,
        ];
        let (xm_data, ym_data) =
            LinearCoefficents2D::get_matrix_data_from_8_points(&x_test_arr, &y_test_arr);
        let _ = linear_coefficents.get_coefficients_from_8_matrix_data(&xm_data, &ym_data, 0.0001f32);

        let (k, b) = linear_coefficents.coefficents();
        assert!(k < -9.5f32 && k > -10.5f32);
        assert!(b > 4.5f32 && b < 5.5f32);
        assert!((linear_coefficents.value(5.0) + 45.2f32).abs() < 1.0f32);
        assert!((linear_coefficents.value(10.3) + 98.1f32).abs() < 1.0f32);
    }
}
