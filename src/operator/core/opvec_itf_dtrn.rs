use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfTranspose {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub visc_a_id: usize,  // viscosity (mu^a)

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub visc_b_id: usize,  // viscosity (mu^b)
}

impl OpVecItfTranspose {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, visc_a_id: usize, vel_b_id: usize, visc_b_id: usize) -> Result<OpVecItfTranspose, VFVChemError> {
        // adds the transpose diffusion term to the combined control volumes of interface nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                     ^^^^^^^^^^              ^^^^^^^^^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b
        // each half is integrated over its closed surface as in OpVecBndTranspose, then the interface flux is removed
        // the one-sided vertex gradient cancels on a closed surface, which keeps the deferred iteration stable

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;

        // create struct
        let mut oper = OpVecItfTranspose::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.visc_a_id = visc_a_id;
        oper.vel_b_id = vel_b_id;
        oper.visc_b_id = visc_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfTranspose {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol div( mu * grad(vel)^T ) dV
        // int_surf mu * grad(vel)^T . n dA
        // sum_ngh mu * grad(vel)^T . n dA + sum_bnd mu * grad(vel)^T . n dA - sum_itf mu * grad(vel)^T . n dA
        // sum_ngh mu_f * grad(vel)^T_f . n * dA + sum_bnd mu_b * grad(vel)^T_b . n * dA - sum_itf mu_b * grad(vel)^T_i . n * dA
        //
        // apply face skewness correction
        // mu_f = (mu_v + mu_n) / 2 + (grad(mu_v) + grad(mu_n)) / 2 . skew
        //
        // apply sub-face interpolation
        // mu_b = w_v * mu_v + w_j * mu_j + w_k * mu_k
        //
        // use linear interpolation on the interface triangle with surface gradient grad_s
        // grad(vel)^T_i . n = grad_s(vel . n) + n * (div(vel)_b - div_s(vel))

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id) in [(self.vel_a_id, self.visc_a_id), (self.vel_b_id, self.visc_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let vel_dom = &vars.vecdom[vel_id];
            let visc_dom = &vars.scldom[visc_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

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
                let vel_v = vel_dom.value[vid_v];
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

                // iterate through boundary sub-faces
                for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                    // get geometry
                    let area = area_b[blc];
                    let (nx, ny, nz) = norm_b[blc];
                    let (wgt_v, wgt_j, wgt_k) = wgt_b[blc];

                    // get neighbor properties
                    let visc_j = visc_dom.value[vid_j];
                    let visc_k = visc_dom.value[vid_k];
                    let (u_x_j, u_y_j, u_z_j, v_x_j, v_y_j, v_z_j, w_x_j, w_y_j, w_z_j) = vel_dom.grad[vid_j];
                    let (u_x_k, u_y_k, u_z_k, v_x_k, v_y_k, v_z_k, w_x_k, w_y_k, w_z_k) = vel_dom.grad[vid_k];

                    // get sub-face viscosity
                    let visc_b = wgt_v * visc_v + wgt_j * visc_j + wgt_k * visc_k;

                    // get sub-face velocity gradient
                    let u_x_b = wgt_v * u_x_v + wgt_j * u_x_j + wgt_k * u_x_k;
                    let u_y_b = wgt_v * u_y_v + wgt_j * u_y_j + wgt_k * u_y_k;
                    let u_z_b = wgt_v * u_z_v + wgt_j * u_z_j + wgt_k * u_z_k;
                    let v_x_b = wgt_v * v_x_v + wgt_j * v_x_j + wgt_k * v_x_k;
                    let v_y_b = wgt_v * v_y_v + wgt_j * v_y_j + wgt_k * v_y_k;
                    let v_z_b = wgt_v * v_z_v + wgt_j * v_z_j + wgt_k * v_z_k;
                    let w_x_b = wgt_v * w_x_v + wgt_j * w_x_j + wgt_k * w_x_k;
                    let w_y_b = wgt_v * w_y_v + wgt_j * w_y_j + wgt_k * w_y_k;
                    let w_z_b = wgt_v * w_z_v + wgt_j * w_z_j + wgt_k * w_z_k;

                    // add to sum
                    sumx += visc_b * (u_x_b * nx + v_x_b * ny + w_x_b * nz) * area;
                    sumy += visc_b * (u_y_b * nx + v_y_b * ny + w_y_b * nz) * area;
                    sumz += visc_b * (u_z_b * nx + v_z_b * ny + w_z_b * nz) * area;

                    // remove the interface flux
                    let gnd_j = dom.nds_gnd[vid_j - dom.num_pts()];
                    let gnd_k = dom.nds_gnd[vid_k - dom.num_pts()];
                    if !(itf.nds_nid.contains_key(&gnd_j) && itf.nds_nid.contains_key(&gnd_k)) { continue; }

                    // get shape function gradients on the interface triangle
                    // grad(phi_v) = cross x (pos_k - pos_j) / |cross|^2 with cross = (pos_j - pos_v) x (pos_k - pos_v)
                    let pos_v = vars.nds.nds_pos[gnd_v];
                    let pos_j = vars.nds.nds_pos[gnd_j];
                    let pos_k = vars.nds.nds_pos[gnd_k];
                    let cross = cross3(sub3(pos_j, pos_v), sub3(pos_k, pos_v));
                    let cross_sq = dot3(cross, cross);
                    let phi_v = scale3(1.0 / cross_sq, cross3(cross, sub3(pos_k, pos_j)));
                    let phi_j = scale3(1.0 / cross_sq, cross3(cross, sub3(pos_v, pos_k)));
                    let phi_k = scale3(1.0 / cross_sq, cross3(cross, sub3(pos_j, pos_v)));

                    // get surface derivatives and sub-face divergence
                    let vel_j = vel_dom.value[vid_j];
                    let vel_k = vel_dom.value[vid_k];
                    let n = (nx, ny, nz);
                    let dnorm = add3(add3(scale3(dot3(vel_v, n), phi_v), scale3(dot3(vel_j, n), phi_j)), scale3(dot3(vel_k, n), phi_k));
                    let dtang = dot3(vel_v, phi_v) + dot3(vel_j, phi_j) + dot3(vel_k, phi_k);
                    let div_b = u_x_b + v_y_b + w_z_b;

                    // subtract from sum
                    sumx -= visc_b * (dnorm.0 + nx * (div_b - dtang)) * area;
                    sumy -= visc_b * (dnorm.1 + ny * (div_b - dtang)) * area;
                    sumz -= visc_b * (dnorm.2 + nz * (div_b - dtang)) * area;
                }

                // load into vector
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -sumx)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -sumy)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 2, vid_b_v, -sumz)?;

            }
        }

        Ok(())
    }
}
