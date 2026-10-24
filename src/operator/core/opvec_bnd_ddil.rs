use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

const NEG_FRAC_TWO_THREE: f64 = -2.0/3.0;

#[derive(Default)]
pub struct OpVecBndDilatatory {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub visc_id: usize,  // viscosity (mu)
}

impl OpVecBndDilatatory {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, visc_id: usize) -> Result<OpVecBndDilatatory, VFVChemError> {
        // adds the dilatatory diffusion term to the control volumes of boundary nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^                              ^^^^^^^^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, visc_id)?;

        // create struct
        let mut oper = OpVecBndDilatatory::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.visc_id = visc_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecBndDilatatory {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( lam * div(vel) * I ) dV
        // int_surf lam * div(vel) * I . n dA
        // int_surf lam * div(vel) * n dA
        // sum_ngh lam * div(vel) * n dA + sum_bnd lam * div(vel) * n dA
        // sum_ngh lam_f * div(vel)_f * n * dA + sum_bnd lam_b * div(vel)_b * n * dA
        //
        // apply face skewness correction
        // lam_f = (lam_v + lam_n) / 2 + (grad(lam_v) + grad(lam_n)) / 2 . skew
        //
        // apply sub-face interpolation
        // lam_b = w_v * lam_v + w_j * lam_j + w_k * lam_k

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let visc_dom = &vars.scldom[self.visc_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tet.ngh_vid[vid_v];
            let face_n = &tet.face_area[vid_v];
            let norm_n = &tet.face_norm[vid_v];
            let skew_n = &tet.face_skew[vid_v];
            let list_b = &tet.bnd_vid[vid_v];
            let area_b = &tet.bnd_area[vid_v];
            let norm_b = &tet.bnd_norm[vid_v];
            let wgt_b = &tet.bnd_wgt[vid_v];

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

            // iterate through boundary sub-faces
            for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                // get geometry
                let area = area_b[blc];
                let (nx, ny, nz) = norm_b[blc];
                let (wgt_v, wgt_j, wgt_k) = wgt_b[blc];

                // get neighbor properties
                let visc_j = visc_dom.value[vid_j];
                let visc_k = visc_dom.value[vid_k];
                let (u_x_j, _, _, _, v_y_j, _, _, _, w_z_j) = vel_dom.grad[vid_j];
                let (u_x_k, _, _, _, v_y_k, _, _, _, w_z_k) = vel_dom.grad[vid_k];
                let div_j = u_x_j + v_y_j + w_z_j;
                let div_k = u_x_k + v_y_k + w_z_k;

                // get sub-face second viscosity
                let visc_b = wgt_v * visc_v + wgt_j * visc_j + wgt_k * visc_k;
                let lam_b = NEG_FRAC_TWO_THREE * visc_b;

                // get sub-face divergence
                let div_b = wgt_v * div_v + wgt_j * div_j + wgt_k * div_k;

                // add to sum
                let coeff = lam_b * div_b * area;
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
