use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfGradient {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub pres_a_id: usize,  // pressure (p^a); unknown

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub pres_b_id: usize,  // pressure (p^b); unknown
}

impl OpVecItfGradient {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, pres_a_id: usize, vel_b_id: usize, pres_b_id: usize) -> Result<OpVecItfGradient, VFVChemError> {
        // adds the pressure gradient term to the combined control volumes of interface nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b
        // the sub-faces on the interface cancel since the pressure is continuous

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, pres_a_id, pres_b_id)?;

        // create struct
        let mut oper = OpVecItfGradient::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.pres_a_id = pres_a_id;
        oper.vel_b_id = vel_b_id;
        oper.pres_b_id = pres_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfGradient {
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
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for pres_id in [self.pres_a_id, self.pres_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[pres_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];

            // get pressure gradient coefficients
            let grad_vec = &tet.grad_coeff;

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // load into matrix
                for &(vid_m, c_x, c_y, c_z) in &grad_vec[vid_v] {
                    self.add_vecscl_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, pres_id, vid_m, -c_x)?;
                    self.add_vecscl_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, pres_id, vid_m, -c_y)?;
                    self.add_vecscl_a(vars, a_triplet, self.vel_b_id, 2, vid_b_v, pres_id, vid_m, -c_z)?;
                }

            }
        }

        Ok(())
    }
}
