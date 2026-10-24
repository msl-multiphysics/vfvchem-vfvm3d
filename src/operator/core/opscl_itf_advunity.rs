use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfAdvectionUnity {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub unk_b_id: usize,  // unknown scalar (c_i^b)

    // scalars - flowing side
    pub unk_id: usize,           // unknown scalar (c_i); unknown
    pub diff_id: Option<usize>,  // diffusion coefficient (G_ii); for the Peclet number

    // vectors - flowing side
    pub vel_id: usize,  // velocity (vel)

    // flow - flowing side
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl OpSclItfAdvectionUnity {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize, diff_id: Option<usize>, vel_id: usize, flw_id: usize) -> Result<OpSclItfAdvectionUnity, VFVChemError> {
        // adds the advection term in terms of the mass flux to the combined control volumes of interface nodes
        // rho * dc_i/dt + rho * vel . grad(c_i) = -div( -G_ij * grad(c_j) ) + R_i
        //                 ^^^^^^^^^^^^^^^^^^^^^
        // covers the half of the control volume on the side of the flow data; loads into the rows of domain b
        // use once for each side with flow; a side without flow has no advection
        // must be used together with OpSclItfDiffusion; the power law scheme corrects its central difference
        // a missing diffusion coefficient sets G_f to zero and falls back to pure upwinding
        // the boundary sub-faces are upwinded from the node and contribute nothing

        // input checking
        // the unknown on the flowing side is either unk_a_id or unk_b_id
        let dom_id = vars.sclflw[flw_id].dom_id;
        let unk_id = if vars.scldom[unk_a_id].dom_id == dom_id { unk_a_id } else { unk_b_id };
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_dom_scldom(vars, dom_id, unk_id)?;
        if let Some(diff_id) = diff_id {
            check_dom_scldom(vars, dom_id, diff_id)?;
        }
        check_dom_vecdom(vars, dom_id, vel_id)?;

        // create struct
        let mut oper = OpSclItfAdvectionUnity::default();
        oper.itf_id = itf_id;
        oper.unk_b_id = unk_b_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.vel_id = vel_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfAdvectionUnity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -rho * vel . grad(c) dV
        // int_surf -rho * (c - c_v) * (vel . n) dA
        // sum_ngh -rho * (c - c_v) * (vel . n) dA + sum_bnd -rho * (c - c_v) * (vel . n) dA
        //
        // introduce the mass flux
        // F = rho * vel * dA
        // sum_ngh -F_f * (c_f - c_v) + sum_bnd -F_b * (c_v - c_v)
        // sum_ngh -F_f * (c_f - c_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v)
        //
        // conservative form minus c_v times the discrete continuity residual of the same F_f
        // uniform c is preserved even if continuity is not satisfied; c is conserved once it is

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];
        let tet = &vars.tet[dom_id];
        let diff_dom = self.diff_id.map(|diff_id| &vars.scldom[diff_id]);
        let flow_dom = &vars.sclflw[self.flw_id];

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
            let flux_n = &flow_dom.flux[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy, sz) = skew_n[loc];

                // get face diffusion coefficient; zero if absent
                let diff_f = match diff_dom {
                    Some(diff_dom) => {
                        let diff_v = diff_dom.value[vid_v];
                        let (diff_x_v, diff_y_v, diff_z_v) = diff_dom.grad[vid_v];
                        let diff_n = diff_dom.value[vid_n];
                        let (diff_x_n, diff_y_n, diff_z_n) = diff_dom.grad[vid_n];
                        let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                        let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n) + sz * (diff_z_v + diff_z_n));  // skew correction
                        if diff_f > 0.0 { diff_f } else { diff_m }  // guard against negative diffusion
                    }
                    None => 0.0,
                };

                // get power law term
                let cond = diff_f * face / dist;
                let flux = flux_n[loc];
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_id, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
