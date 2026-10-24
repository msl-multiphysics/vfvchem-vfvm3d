use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfAdvection {
    // interface
    pub itf_id: usize,

    // scalars - rows
    pub unk_b_id: usize,  // unknown scalar (c_i^b)

    // scalars - flowing side
    pub unk_id: usize,   // unknown scalar (c_i); unknown
    pub wgt_id: usize,   // weight coefficient (m_i)
    pub diff_id: usize,  // diffusion coefficient (G_ii); for the Peclet number
    pub den_id: usize,   // density (rho); converts the mass flux

    // vectors - flowing side
    pub vel_id: usize,  // velocity (vel)

    // flow - flowing side
    pub flw_id: usize,  // flow data; for the Rhie-Chow mass flux
}

impl OpSclItfAdvection {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize, wgt_id: usize, diff_id: usize, vel_id: usize, den_id: usize, flw_id: usize) -> Result<OpSclItfAdvection, VFVChemError> {
        // adds the advection term to the combined control volumes of interface nodes
        // m_i * dc_i/dt + m_i * vel . grad(c_i) = -div( -G_ii * grad(c_i) ) + R_i
        //                 ^^^^^^^^^^^^^^^^^^^^^
        // covers the half of the control volume on the side of the flow data; loads into the rows of domain b
        // use once for each side with flow; a side without flow has no advection
        // must be used together with OpSclItfDiffusion; the power law scheme corrects its central difference
        // the boundary sub-faces are upwinded from the node and contribute nothing

        // input checking
        // the unknown on the flowing side is either unk_a_id or unk_b_id
        let dom_id = vars.sclflw[flw_id].dom_id;
        let unk_id = if vars.scldom[unk_a_id].dom_id == dom_id { unk_a_id } else { unk_b_id };
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;
        check_dom_vecdom(vars, dom_id, vel_id)?;

        // create struct
        let mut oper = OpSclItfAdvection::default();
        oper.itf_id = itf_id;
        oper.unk_b_id = unk_b_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;
        oper.diff_id = diff_id;
        oper.den_id = den_id;
        oper.vel_id = vel_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -m * vel . grad(c) dV
        // int_surf -m * (c - c_v) * (vel . n) dA
        // sum_ngh -m * (c - c_v) * (vel . n) dA + sum_bnd -m * (c - c_v) * (vel . n) dA
        //
        // introduce weight flux from the mass flux
        // F = m * vel * dA = (m / rho) * rho * vel * dA
        // sum_ngh -F_f * (c_f - c_v) + sum_bnd -F_b * (c_v - c_v)
        // sum_ngh -F_f * (c_f - c_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v)
        //
        // non-conservative form; uniform c is preserved even if continuity is not satisfied

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_id = vars.scldom[self.unk_id].dom_id;
        let dom = &vars.dom[dom_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];
        let tet = &vars.tet[dom_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        let diff_dom = &vars.scldom[self.diff_id];
        let den_dom = &vars.scldom[self.den_id];
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

            // get vertex properties
            let wgt_v = wgt_dom.value[vid_v];
            let (wgt_x_v, wgt_y_v, wgt_z_v) = wgt_dom.grad[vid_v];
            let diff_v = diff_dom.value[vid_v];
            let (diff_x_v, diff_y_v, diff_z_v) = diff_dom.grad[vid_v];
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v, den_z_v) = den_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let dist = dist_n[loc];
                let (sx, sy, sz) = skew_n[loc];

                // get neighbor properties
                let wgt_n = wgt_dom.value[vid_n];
                let (wgt_x_n, wgt_y_n, wgt_z_n) = wgt_dom.grad[vid_n];
                let diff_n = diff_dom.value[vid_n];
                let (diff_x_n, diff_y_n, diff_z_n) = diff_dom.grad[vid_n];
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n, den_z_n) = den_dom.grad[vid_n];

                // get face weight, density, and diffusion coefficient
                let wgt_f = 0.5 * (wgt_v + wgt_n) + 0.5 * (sx * (wgt_x_v + wgt_x_n) + sy * (wgt_y_v + wgt_y_n) + sz * (wgt_z_v + wgt_z_n));  // skew-corrected arithmetic mean
                let den_f = 0.5 * (den_v + den_n) + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n) + sz * (den_z_v + den_z_n));  // skew-corrected arithmetic mean
                let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n) + sz * (diff_z_v + diff_z_n));  // skew correction
                let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion

                // get power law term
                let cond = diff_f * face / dist;
                let flux = if den_f != 0.0 { wgt_f / den_f * flux_n[loc] } else { 0.0 };
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
