use crate::base::prelude::*;
use crate::linear::prelude::*;
use faer::linalg::solvers::Solve;
use faer::sparse::{Argsort, Pair, SparseColMat, SymbolicSparseColMat, Triplet};
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
use faer::{Col, Par};
use std::num::NonZeroUsize;

pub struct LinearLu {
    pub symbolic: SymbolicSparseColMat<usize, usize, usize>,
    pub argsort: Argsort<usize>,
    pub symbolic_lu: SymbolicLu<usize>,
}

impl LinearLu {
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // placeholder default
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(0, 0, &[]).map_err(|e| _symbolic_err(e))?;
        let symbolic_lu = SymbolicLu::try_new(symbolic.as_ref()).map_err(|e| _symbolic_lu_err(e))?;
        let lin = LinearLu { symbolic, argsort, symbolic_lu };

        // set number of threads
        let num_thread_nz = NonZeroUsize::new(num_thread).unwrap();
        faer::set_global_parallelism(Par::Rayon(num_thread_nz));

        // return
        Ok(Box::new(lin))
    }
}

impl LinearBase for LinearLu {
    fn update_sparse(&mut self, a_triplet: &Vec<Triplet<usize, usize, f64>>, mat_size: usize) -> Result<(), VFVChemError> {
        // collect (row, col) pairs
        let a_pair: Vec<Pair<usize, usize>> = a_triplet.iter().map(|t| Pair::new(t.row, t.col)).collect();

        // symbolic sparse matrix and LU factorization
        // reused for every later matrix with this sparsity pattern
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(mat_size, mat_size, &a_pair).map_err(|e| _symbolic_err(e))?;
        let symbolic_lu = SymbolicLu::try_new(symbolic.as_ref()).map_err(|e| _symbolic_lu_err(e))?;
        self.symbolic = symbolic;
        self.argsort = argsort;
        self.symbolic_lu = symbolic_lu;

        Ok(())
    }

    fn build_matrix(&self, a_triplet: &Vec<Triplet<usize, usize, f64>>) -> Result<SparseColMat<usize, f64>, VFVChemError> {
        // collect values
        let a_vals: Vec<f64> = a_triplet.iter().map(|t| t.val).collect();

        // assemble matrix using symbolic factorization
        let a_mat = SparseColMat::new_from_argsort(self.symbolic.clone(), &self.argsort, &a_vals).map_err(|e| _matrix_err(e))?;

        Ok(a_mat)
    }

    fn solve_matrix(&self, a_mat: &SparseColMat<usize, f64>, b_vec: &Col<f64>, _x_init: &Col<f64>, _mat_size: usize) -> Result<Col<f64>, VFVChemError> {
        // numeric factorization using the cached symbolic LU
        // cloning SymbolicLu only copies its reference count
        let lu = Lu::try_new_with_symbolic(self.symbolic_lu.clone(), a_mat.as_ref()).map_err(|e| _lu_err(e))?;

        // matrix solve
        let x_vec = lu.solve(b_vec);

        Ok(x_vec)
    }
}

fn _symbolic_err(source: faer::sparse::CreationError) -> VFVChemError {
    VFVChemError::LinearLuSymbolicError { source }
}

fn _matrix_err(source: faer::sparse::FaerError) -> VFVChemError {
    VFVChemError::LinearLuMatrixError { source }
}

fn _symbolic_lu_err(source: faer::sparse::FaerError) -> VFVChemError {
    VFVChemError::LinearLuSymbolicLuError { source }
}

fn _lu_err(source: faer::sparse::linalg::LuError) -> VFVChemError {
    VFVChemError::LinearLuLuError { source }
}
