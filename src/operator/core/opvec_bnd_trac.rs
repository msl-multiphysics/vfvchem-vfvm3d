use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndTraction {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p); unknown
    pub val_id: usize,   // outflow pressure (p_set)
}

impl OpVecBndTraction {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize, val_id: usize) -> Result<OpVecBndTraction, VFVChemError> {
        // adds the pressure gradient term to the control volumes of nodes on a traction outflow boundary
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^
        // sets the traction on the boundary sub-faces
        // (-p * I + mu * grad(vel)) . n = -p_set * n
        // the viscous boundary operators omit the boundary sub-faces, so the viscous flux there is zero
        // replaces OpSclBndDirichlet on pressure and OpVecBndGradient; pressure stays unknown on the boundary
        // and the continuity equation of the boundary nodes is kept

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;
        check_bnd_sclbnd(vars, bnd_id, val_id)?;

        // create struct
        let mut oper = OpVecBndTraction::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.val_id = val_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndTraction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -grad(p) dV
        // int_surf -p n dA
        // sum_ngh -p_f * n * dA + sum_bnd -p_b * n * dA
        // sum_ngh -p_f * n * dA + sum_bnd -(w_v * p_v + w_j * p_j + w_k * p_k) * n * dA
        //
        // expand the least-squares gradients; see Tetrahedron
        // sum_m -(c_x, c_y, c_z) * p_m
        //
        // on sub-faces along this boundary, replace the interpolated pressure with the outflow pressure
        // sum_bnd (w_v * p_v + w_j * p_j + w_k * p_k) * n * dA - (w_v * p_set_v + w_j * p_set_j + w_k * p_set_k) * n * dA

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let val_bnd = &vars.sclbnd[self.val_id];

        // get pressure gradient coefficients
        let grad_vec = &tet.grad_coeff;

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_b = &tet.bnd_vid[vid_v];
            let area_b = &tet.bnd_area[vid_v];
            let norm_b = &tet.bnd_norm[vid_v];
            let wgt_b = &tet.bnd_wgt[vid_v];

            // get vertex properties
            let val_v = val_bnd.value[idx_v];

            // load gradient into matrix
            for &(vid_m, c_x, c_y, c_z) in &grad_vec[vid_v] {
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_m, -c_x)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_m, -c_y)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 2, vid_v, self.pres_id, vid_m, -c_z)?;
            }

            // iterate through boundary sub-faces
            for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                // get indices
                let gnd_j = dom.nds_gnd[vid_j - dom.num_pts()];
                let gnd_k = dom.nds_gnd[vid_k - dom.num_pts()];

                // skip if a neighbor is on another boundary
                let Some(&idx_j) = bnd.nds_idx.get(&gnd_j) else { continue };
                let Some(&idx_k) = bnd.nds_idx.get(&gnd_k) else { continue };

                // get geometry
                let area = area_b[blc];
                let (nx, ny, nz) = norm_b[blc];
                let (w_v, w_j, w_k) = wgt_b[blc];

                // get neighbor properties
                let val_j = val_bnd.value[idx_j];
                let val_k = val_bnd.value[idx_k];

                // load into matrix and vector
                let val_b = w_v * val_v + w_j * val_j + w_k * val_k;
                for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                    self.add_vecscl_a(vars, a_triplet, self.vel_id, comp, vid_v, self.pres_id, vid_v, w_v * area * n_c)?;
                    self.add_vecscl_a(vars, a_triplet, self.vel_id, comp, vid_v, self.pres_id, vid_j, w_j * area * n_c)?;
                    self.add_vecscl_a(vars, a_triplet, self.vel_id, comp, vid_v, self.pres_id, vid_k, w_k * area * n_c)?;
                    self.add_vec_b(vars, b_vec, self.vel_id, comp, vid_v, val_b * area * n_c)?;
                }
            }
        }

        Ok(())
    }
}
