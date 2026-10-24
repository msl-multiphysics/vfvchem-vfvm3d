use crate::base::prelude::*;
use crate::linear::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::time::{Duration, Instant};

pub trait SteadyBase {
    // to implement in specific solvers
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError>;  // set start indices, create derived data, and compute matrix sizes
    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError>;  // initializes operators
    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError>;  // assemble global matrix
    fn initial_dirichlet(&self, vars: &mut Variables);  // load Dirichlet boundary values into unknowns
    fn update_flow(&self, vars: &mut Variables);  // update the face mass flux of the flow data; empty for solvers without flow

    // solver function
    fn solve(
        &mut self, vars: &mut Variables, lin: &mut Box<dyn LinearBase>,
        max_iter: usize, tol: f64, damp: f64
    ) -> Result<(), VFVChemError> {
        // error handling
        if max_iter == 0 { return Err(VFVChemError::SolverNonPositiveMaxIter { max_iter }); }
        if tol.is_nan() || tol <= 0.0 { return Err(VFVChemError::SolverNonPositiveTolerance { tol }); }
        if damp.is_nan() || damp <= 0.0 || damp > 1.0 { return Err(VFVChemError::SolverInvalidDamping { damp }); }

        // start solver
        let time_start = Instant::now();
        println!("Starting steady-state solver.");

        // initialize time measurement
        let mut time_dual = Duration::ZERO;
        let mut time_asm = Duration::ZERO;
        let mut time_lin = Duration::ZERO;
        let mut time_wrt = Duration::ZERO;

        // initialize dual control volumes
        let time_0 = Instant::now();
        Self::_update_tet(vars)?;
        let time_1 = Instant::now();
        time_dual += time_1 - time_0;

        // initialize operators
        let mat_size = self.initial_start(vars)?;
        self.initial_operator(vars)?;
        let fields = field_norms(vars, &[]);  // unknown fields for the scaled residual

        // initialize matrices
        let mut a_triplet: Vec<Triplet<usize, usize, f64>> = Vec::new();
        let mut b_vec: Col<f64> = Col::zeros(mat_size);
        let mut x_udmp_vec: Col<f64>;  // undamped result
        let mut x_iter_vec: Col<f64> = Col::zeros(mat_size);  // result between iterations

        // load Dirichlet boundary values into the initial unknowns
        Self::_update_func(vars, 0.0);
        self.initial_dirichlet(vars);

        // load initial unknown value to solution vector
        Self::_unkn_to_solu(vars, &mut x_iter_vec);

        // initial assembly of A and b
        // used for sparsity pattern and first iteration
        Self::_update_func(vars, 0.0);
        Self::_update_grad(vars);
        self.update_flow(vars);
        self.assemble_matrix(vars, &mut a_triplet, &mut b_vec)?;
        lin.update_sparse(&a_triplet, mat_size)?;
        let mut a_mat = lin.build_matrix(&a_triplet)?;
        let num_triplet = a_triplet.len();
        let time_2 = Instant::now();
        time_asm += time_2 - time_1;

        // iterate to convergence
        let mut is_conv = false;
        let mut res = 0.0;
        for iter in 0..max_iter {
            // solve for x_udmp in A_k x_udmp = b_k
            // solve for x_{k+1} in x_damp = (1 - damp) * x_k + damp * x_udmp
            let time_i0 = Instant::now();
            x_udmp_vec = lin.solve_matrix(&a_mat, &b_vec, &x_iter_vec, mat_size)?;
            let time_i1 = Instant::now();
            let x_damp_vec = (1.0 - damp) * &x_iter_vec + damp * &x_udmp_vec;
            time_lin += time_i1 - time_i0;

            // assemble A_{k+1} and b_{k+1} using x_{k+1}
            a_triplet = Vec::with_capacity(num_triplet);  // reset A
            b_vec = Col::zeros(mat_size);  // reset b
            Self::_solu_to_unkn(vars, &x_damp_vec);
            Self::_update_func(vars, 0.0);
            Self::_update_grad(vars);
            self.update_flow(vars);
            self.assemble_matrix(vars, &mut a_triplet, &mut b_vec)?;
            a_mat = lin.build_matrix(&a_triplet)?;
            let time_i2 = Instant::now();
            time_asm += time_i2 - time_i1;

            // compute tolerance: max over fields of rms(x_new - x_old) / scale (see scaling.rs)
            res = scaled_residual(vars, &fields, &x_iter_vec, &x_damp_vec).0;
            println!("Iteration: {iter}; Residual: {res:.4e}");

            // error if diverged; NaN fails every comparison, so it would otherwise run to max_iter
            if !res.is_finite() {
                return Err(VFVChemError::SteadyNonFiniteResidual { iter });
            }

            // update iteration vector
            x_iter_vec = x_damp_vec;

            // tolerance check
            if res < tol {
                is_conv = true;
                break;
            }
        }

        // error if not converged
        if !is_conv {
            return Err(VFVChemError::SteadyDidNotConverge { max_iter, res });
        }

        // write variables
        let time_3  = Instant::now();
        Self::_write(vars, 0)?;
        let time_end  = Instant::now();
        time_wrt += time_end - time_3;

        // output time measurement
        let time_total = time_end - time_start;
        println!("Solution completed!");
        println!("Total: {:.6} s", time_total.as_secs_f64());
        println!("  Dual Mesh       : {:.6} s", time_dual.as_secs_f64());
        println!("  Matrix Assembly : {:.6} s", time_asm.as_secs_f64());
        println!("  Matrix Solve    : {:.6} s", time_lin.as_secs_f64());
        println!("  File Output     : {:.6} s", time_wrt.as_secs_f64());

        Ok(())
    }

    fn _unkn_to_solu(vars: &Variables, x_vec: &mut Col<f64>) {
        // iterate through unknown scalars
        for scldom in &vars.scldom {
            match scldom.scldom_type {
                SclDomType::Unknown { start } => {
                    // iterate through points
                    for (tid, &value) in scldom.value.iter().enumerate() {
                        let xid = start + tid;
                        x_vec[xid] = value;
                    }
                }
                _ => {continue;}  // skip if not unknown
            }
        }

        // iterate through unknown vectors
        for vecdom in &vars.vecdom {
            match vecdom.vecdom_type {
                VecDomType::Unknown { start } => {
                    // iterate through points
                    let num_tot = vecdom.value.len();
                    for (tid, &value) in vecdom.value.iter().enumerate() {
                        let xid_x = start + tid;
                        let xid_y = start + num_tot + tid;
                        let xid_z = start + 2 * num_tot + tid;
                        let (value_x, value_y, value_z) = value;
                        x_vec[xid_x] = value_x;
                        x_vec[xid_y] = value_y;
                        x_vec[xid_z] = value_z;
                    }
                }
                _ => {continue;}  // skip if not unknown
            }
        }

    }

    fn _solu_to_unkn(vars: &mut Variables, x_vec: &Col<f64>) {
        // iterate through unknown scalars
        for scldom in &mut vars.scldom {
            match scldom.scldom_type {
                SclDomType::Unknown { start } => {
                    // iterate through points
                    let num_tot = scldom.value.len();
                    for tid in 0..num_tot {
                        let xid = start + tid;
                        let value = x_vec[xid];
                        scldom.value[tid] = value;
                    }
                }
                _ => {continue;}  // skip if not unknown
            }
        }

        // iterate through unknown vectors
        for vecdom in &mut vars.vecdom {
            match vecdom.vecdom_type {
                VecDomType::Unknown { start } => {
                    // iterate through points
                    let num_tot = vecdom.value.len();
                    for tid in 0..num_tot {
                        let xid_x = start + tid;
                        let xid_y = start + num_tot + tid;
                        let xid_z = start + 2 * num_tot + tid;
                        let value_x = x_vec[xid_x];
                        let value_y = x_vec[xid_y];
                        let value_z = x_vec[xid_z];
                        vecdom.value[tid] = (value_x, value_y, value_z);
                    }
                }
                _ => {continue;}  // skip if not unknown
            }
        }

    }

    fn _write(vars: &Variables, ts: usize) -> Result<(), VFVChemError> {
        for scldom in vars.scldom.iter() {
            scldom.write(&vars.pts, &vars.nds, &vars.dom[scldom.dom_id], &vars.tet[scldom.dom_id], ts)?;
        }
        for sclbnd in vars.sclbnd.iter() {
            sclbnd.write(&vars.nds, &vars.bnd[sclbnd.bnd_id], ts)?;
        }
        for sclitf in vars.sclitf.iter() {
            sclitf.write(&vars.nds, &vars.itf[sclitf.itf_id], ts)?;
        }
        for vecdom in vars.vecdom.iter() {
            vecdom.write(&vars.pts, &vars.nds, &vars.dom[vecdom.dom_id], &vars.tet[vecdom.dom_id], ts)?;
        }
        for vecbnd in vars.vecbnd.iter() {
            vecbnd.write(&vars.nds, &vars.bnd[vecbnd.bnd_id], ts)?;
        }
        for vecitf in vars.vecitf.iter() {
            vecitf.write(&vars.nds, &vars.itf[vecitf.itf_id], ts)?;
        }
        Ok(())
    }

    fn _update_tet(vars: &mut Variables) -> Result<(), VFVChemError> {
        for tet in vars.tet.iter_mut() {
            tet.update_tet(&vars.nds, &vars.pts, &vars.dom[tet.dom_id])?;
        }
        Ok(())
    }

    fn _update_func(vars: &mut Variables, t: f64) {
        for scldom_id in 0..vars.scldom.len() {
            ScalarDomain::update_func(vars, scldom_id, t);
        }
        for sclbnd_id in 0..vars.sclbnd.len() {
            ScalarBoundary::update_func(vars, sclbnd_id, t);
        }
        for sclitf_id in 0..vars.sclitf.len() {
            ScalarInterface::update_func(vars, sclitf_id, t);
        }
        for vecdom_id in 0..vars.vecdom.len() {
            VectorDomain::update_func(vars, vecdom_id, t);
        }
        for vecbnd_id in 0..vars.vecbnd.len() {
            VectorBoundary::update_func(vars, vecbnd_id, t);
        }
        for vecitf_id in 0..vars.vecitf.len() {
            VectorInterface::update_func(vars, vecitf_id, t);
        }
    }

    fn _update_grad(vars: &mut Variables) {
        for scldom in vars.scldom.iter_mut() {
            scldom.update_grad(&vars.tet[scldom.dom_id]);
        }
        for vecdom in vars.vecdom.iter_mut() {
            vecdom.update_grad(&vars.tet[vecdom.dom_id]);
        }
    }

}
