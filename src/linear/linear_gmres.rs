use crate::base::prelude::*;
use crate::linear::prelude::*;
use faer::sparse::{Argsort, Pair, SparseColMat, SymbolicSparseColMat, Triplet};
use faer::{Col, Par};
use faer_gmres::{JacobiPreconLinOp, restarted_gmres};
use std::num::NonZeroUsize;

pub struct LinearGmres {
    pub symbolic: SymbolicSparseColMat<usize, usize, usize>,
    pub argsort: Argsort<usize>,
    tol: f64,
    restart: usize,
    max_outer: usize,
}

impl LinearGmres {
    // default: Jacobi-preconditioned restarted GMRES
    // restart = Krylov subspace size per outer iteration (memory ~ O(n * restart))
    // max_outer = number of restart cycles
    pub fn new(num_thread: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        Self::new_with_options(num_thread, 1e-6, 50, 200)
    }

    pub fn new_with_options(num_thread: usize, tol: f64, restart: usize, max_outer: usize) -> Result<Box<dyn LinearBase>, VFVChemError> {
        // error handling
        if num_thread == 0 {
            return Err(VFVChemError::LinearGmresNonPositiveThreads { num_thread });
        }
        if tol.is_nan() || tol <= 0.0 {
            return Err(VFVChemError::LinearGmresNonPositiveTolerance { tol });
        }
        if restart == 0 {
            return Err(VFVChemError::LinearGmresNonPositiveRestart { restart });
        }
        if max_outer == 0 {
            return Err(VFVChemError::LinearGmresNonPositiveMaxOuter { max_outer });
        }

        // placeholder default
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(0, 0, &[]).map_err(|e| _symbolic_err(e))?;
        let lin = LinearGmres { symbolic, argsort, tol, restart, max_outer };

        // set number of threads
        let num_thread_nz = NonZeroUsize::new(num_thread).unwrap();
        faer::set_global_parallelism(Par::Rayon(num_thread_nz));

        // return
        Ok(Box::new(lin))
    }
}

impl LinearBase for LinearGmres {
    fn update_sparse(&mut self, a_triplet: &Vec<Triplet<usize, usize, f64>>, mat_size: usize) -> Result<(), VFVChemError> {
        // collect (row, col) pairs
        let a_pair: Vec<Pair<usize, usize>> = a_triplet.iter().map(|t| Pair::new(t.row, t.col)).collect();

        // symbolic sparse matrix
        // reused for every later matrix with this sparsity pattern
        let (symbolic, argsort) = SymbolicSparseColMat::try_new_from_indices(mat_size, mat_size, &a_pair).map_err(|e| _symbolic_err(e))?;
        self.symbolic = symbolic;
        self.argsort = argsort;

        Ok(())
    }

    fn build_matrix(&self, a_triplet: &Vec<Triplet<usize, usize, f64>>) -> Result<SparseColMat<usize, f64>, VFVChemError> {
        // collect values
        let a_vals: Vec<f64> = a_triplet.iter().map(|t| t.val).collect();

        // assemble matrix using symbolic factorization
        let a_mat = SparseColMat::new_from_argsort(self.symbolic.clone(), &self.argsort, &a_vals).map_err(|e| _matrix_err(e))?;

        Ok(a_mat)
    }

    fn solve_matrix(&self, a_mat: &SparseColMat<usize, f64>, b_vec: &Col<f64>, x_init: &Col<f64>, _mat_size: usize) -> Result<Col<f64>, VFVChemError> {
        // Jacobi-preconditioned restarted GMRES, starting from x_init
        let mut x_vec = x_init.clone();
        let jacobi = JacobiPreconLinOp::new(a_mat.as_ref());

        restarted_gmres(a_mat.as_ref(), b_vec.as_mat(), x_vec.as_mat_mut(), self.restart, self.max_outer, self.tol, Some(&jacobi)).map_err(|e| VFVChemError::LinearGmresDidNotConverge { message: format!("{e:?}") })?;

        Ok(x_vec)
    }
}

fn _symbolic_err(source: faer::sparse::CreationError) -> VFVChemError {
    VFVChemError::LinearGmresSymbolicError { source }
}

fn _matrix_err(source: faer::sparse::FaerError) -> VFVChemError {
    VFVChemError::LinearGmresMatrixError { source }
}
