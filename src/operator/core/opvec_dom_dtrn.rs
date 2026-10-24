use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomTranspose {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecDomTranspose {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecDomTranspose, VFVChemError> {
        // adds the transpose diffusion term to the Navier-Stokes equation
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^              ^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, visc_id)?;

        // create struct
        let mut oper = OpVecDomTranspose::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( mu * grad(vel)^T ) dV
        // int_surf mu * grad(vel)^T . n dA
        // sum_ngh mu * grad(vel)^T . n dA
        // sum_ngh mu_f * grad(vel)^T_f . n * dA
        //
        // apply face skewness correction
        // mu_f = (mu_v + mu_n) / 2 + (grad(mu_v) + grad(mu_n)) / 2 . skew

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
            let (u_x_v, u_y_v, u_z_v, v_x_v, v_y_v, v_z_v, w_x_v, w_y_v, w_z_v) = vel_dom.grad[vid_v];

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
                let (u_x_n, u_y_n, u_z_n, v_x_n, v_y_n, v_z_n, w_x_n, w_y_n, w_z_n) = vel_dom.grad[vid_n];

                // get face viscosity
                let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n) + sz * (visc_z_v + visc_z_n));  // skew correction
                let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                // get face velocity gradient; average of both sides
                let u_x_f = 0.5 * (u_x_v + u_x_n);
                let u_y_f = 0.5 * (u_y_v + u_y_n);
                let u_z_f = 0.5 * (u_z_v + u_z_n);
                let v_x_f = 0.5 * (v_x_v + v_x_n);
                let v_y_f = 0.5 * (v_y_v + v_y_n);
                let v_z_f = 0.5 * (v_z_v + v_z_n);
                let w_x_f = 0.5 * (w_x_v + w_x_n);
                let w_y_f = 0.5 * (w_y_v + w_y_n);
                let w_z_f = 0.5 * (w_z_v + w_z_n);

                // add to sum
                sumx += visc_f * (u_x_f * nx + v_x_f * ny + w_x_f * nz) * face;
                sumy += visc_f * (u_y_f * nx + v_y_f * ny + w_y_f * nz) * face;
                sumz += visc_f * (u_z_f * nx + v_z_f * ny + w_z_f * nz) * face;
            }

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -sumx)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -sumy)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 2, vid_v, -sumz)?;

        }

        Ok(())
    }
}
