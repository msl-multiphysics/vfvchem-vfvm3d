use crate::base::prelude::*;
use faer::Col;
use faer::sparse::{SparseColMat, Triplet};

pub trait LinearBase {
    // to implement in specific solvers
    fn update_sparse(&mut self, a_triplet: &Vec<Triplet<usize, usize, f64>>, mat_size: usize) -> Result<(), VFVChemError>;
    fn build_matrix(&self, a_triplet: &Vec<Triplet<usize, usize, f64>>) -> Result<SparseColMat<usize, f64>, VFVChemError>;
    fn solve_matrix(&self, a_mat: &SparseColMat<usize, f64>, b_vec: &Col<f64>, x_init: &Col<f64>, mat_size: usize) -> Result<Col<f64>, VFVChemError>;
}
