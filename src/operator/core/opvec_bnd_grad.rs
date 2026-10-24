use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndGradient {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p); unknown
}

impl OpVecBndGradient {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize) -> Result<OpVecBndGradient, VFVChemError> {
        // adds the pressure gradient term to the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;

        // create struct
        let mut oper = OpVecBndGradient::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndGradient {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -grad(p) dV
        // int_surf -p n dA
        // sum_ngh -p_f * n * dA + sum_bnd -p_b * n * dA
        // sum_ngh -p_f * n * dA + sum_bnd -(w_v * p_v + w_j * p_j + w_k * p_k) * n * dA
        //
        // apply face skewness correction
        // p_f = (p_v + p_n) / 2 + (grad(p_v) + grad(p_n)) / 2 . skew
        //
        // expand the least-squares gradients; see Tetrahedron
        // sum_m -(c_x, c_y, c_z) * p_m

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];

        // get pressure gradient coefficients
        let grad_vec = &tet.grad_coeff;

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // load into matrix
            for &(vid_m, c_x, c_y, c_z) in &grad_vec[vid_v] {
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_m, -c_x)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_m, -c_y)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 2, vid_v, self.pres_id, vid_m, -c_z)?;
            }

        }

        Ok(())
    }
}
