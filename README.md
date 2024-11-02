# linear_curve_fit

[![Rust](https://github.com/Joker2770/linear_curve_fit/actions/workflows/rust.yml/badge.svg)](https://github.com/Joker2770/linear_curve_fit/actions/workflows/rust.yml)

A `no_std` and no heap memory library for linear curve coefficents calculation.

## example

```rust
    use linear_curve_fit::linear_curve::LinearCoefficents2D;

    let mut linear_coefficents = LinearCoefficents2D::default();
    let x_test_arr = [-2.8f32, -1.6f32, -0.5f32, 5.0f32, 5.4f32, 6.7f32, 10.3f32, 13.8f32];
    let y_test_arr = [33.1f32, 21.1f32, 9.9f32, -45.2f32, -49.1f32, -61.9f32, -98.1f32,-132.99f32];
    let (xm_data, ym_data) = LinearCoefficents2D::get_matrix_data_from_8_points(&x_test_arr, &y_test_arr);
    linear_coefficents.get_coefficients_from_8_matrix_data(&xm_data, &ym_data, 0.0001f32);

    assert!(linear_coefficents.coefficents().0 < -9.5f32);
    assert!(linear_coefficents.coefficents().1 > 4.5f32);
    assert!((linear_coefficents.value(5.0)  + 45.2f32).abs() < 1.0f32);
    assert!((linear_coefficents.value(10.3)  + 98.1f32).abs() < 1.0f32);

```
