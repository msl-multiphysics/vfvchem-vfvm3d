use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfDiffusion {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub visc_a_id: usize,  // viscosity (mu^a)

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub visc_b_id: usize,  // viscosity (mu^b)
}

impl OpVecItfDiffusion {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, visc_a_id: usize, vel_b_id: usize, visc_b_id: usize) -> Result<OpVecItfDiffusion, VFVChemError> {
        // adds the diffusion term to the combined control volumes of interface nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // the interface flux leaves one half and enters the other, so it cancels

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;

        // create struct
        let mut oper = OpVecItfDiffusion::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.visc_a_id = visc_a_id;
        oper.vel_b_id = vel_b_id;
        oper.visc_b_id = visc_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( mu * grad(vel) ) dV
        // int_surf mu * grad(vel) . n dA
        // sum_ngh mu * grad(vel) . n dA
        //
        // apply normal derivative
        // sum_ngh mu_f * (vel_n - vel_v) / dn * dA

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id) in [(self.vel_a_id, self.visc_a_id), (self.vel_b_id, self.visc_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let visc_dom = &vars.scldom[visc_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let list_n = &tet.ngh_vid[vid_v];
                let dist_n = &tet.ngh_dist[vid_v];
                let face_n = &tet.face_area[vid_v];
                let skew_n = &tet.face_skew[vid_v];

                // get vertex properties
                let visc_v = visc_dom.value[vid_v];
                let (visc_x_v, visc_y_v, visc_z_v) = visc_dom.grad[vid_v];

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (sx, sy, sz) = skew_n[loc];

                    // get neighbor properties
                    let visc_n = visc_dom.value[vid_n];
                    let (visc_x_n, visc_y_n, visc_z_n) = visc_dom.grad[vid_n];

                    // get face viscosity
                    let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                    let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n) + sz * (visc_z_v + visc_z_n));  // skew correction
                    let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                    // load into matrix
                    let coeff = visc_f * face / dist;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_v, -coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_v, -coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 2, vid_b_v, vel_id, 2, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 2, vid_b_v, vel_id, 2, vid_v, -coeff)?;
                }

            }
        }

        Ok(())
    }
}
