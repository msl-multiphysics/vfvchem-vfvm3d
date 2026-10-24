use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct OpVecDomDilatatory {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecDomDilatatory {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecDomDilatatory, VFVChemError> {
        // adds the dilatatory diffusion term to the Navier-Stokes equation
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^                              ^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = OpVecDomDilatatory::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( lam * div(vel) * I ) dV
        // int_surf lam * div(vel) * I . n dA
        // int_surf lam * div(vel) * n dA
        // sum_ngh lam * div(vel) * n dA
        // sum_ngh lam_f * div(vel)_f * n * dA
        //
        // apply face skewness correction
        // lam_f = (lam_v + lam_n) / 2 + (grad(lam_v) + grad(lam_n)) / 2 . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tet.ngh_vid[vid_v];
            let face_n = &tet.face_area[vid_v];
            let norm_n = &tet.face_norm[vid_v];
            let skew_n = &tet.face_skew[vid_v];

            // get vertex properties
            let visc_v = visc_dom.value[vid_v];
            let (visc_x_v, visc_y_v, visc_z_v) = visc_dom.grad[vid_v];
            let (u_x_v, _, _, _, v_y_v, _, _, _, w_z_v) = vel_dom.grad[vid_v];
            let div_v = u_x_v + v_y_v + w_z_v;

            // iterate through faces
            let mut sumx = 0.0;
            let mut sumy = 0.0;
            let mut sumz = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny, nz) = norm_n[loc];
                let (sx, sy, sz) = skew_n[loc];

                // get neighbor properties
                let visc_n = visc_dom.value[vid_n];
                let (visc_x_n, visc_y_n, visc_z_n) = visc_dom.grad[vid_n];
                let (u_x_n, _, _, _, v_y_n, _, _, _, w_z_n) = vel_dom.grad[vid_n];
                let div_n = u_x_n + v_y_n + w_z_n;

                // get face second viscosity
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n) + sz * (visc_z_v + visc_z_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity
                let lam_f = NEG_FRAC_TWO_THREE * visc_f;

                // get face divergence
                let div_f = 0.5 * (div_v + div_n);

                // add to sum
                let coeff = lam_f * div_f * face;
                sumx += coeff * nx;
                sumy += coeff * ny;
                sumz += coeff * nz;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 2, vid_v, -sumz)?;

        }

        Ok(())
    }
}
