use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomDiffusion {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecDomDiffusion, VFVChemError> {
        // adds the diffusion term to the Navier-Stokes equation
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = OpVecDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( mu * grad(vel) ) dV
        // int_surf mu * grad(vel) . n dA
        // sum_ngh mu * grad(vel) . n dA
        //
        // apply normal derivative
        // sum_ngh mu_f * (vel_n - vel_v) / dn * dA

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_n,  coeff)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -coeff)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_n,  coeff)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -coeff)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 2, vid_v, self.vel_id, 2, vid_n,  coeff)?;
                self.add_vecvec_a(vars, a_triplet, self.vel_id, 2, vid_v, self.vel_id, 2, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
