use thiserror::Error;

#[derive(Debug, Error)]
pub enum CustomError {
    #[error("SVD solve failed")]
    SvdFailed,
    #[error("Matrix size not match")]
    MatrixSizeNotMatch,
}
