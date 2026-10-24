use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfDiffusion {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,   // unknown scalar (c_i^a)
    pub diff_a_id: usize,  // diffusion coefficient (G_ij^a)
    pub drv_a_id: usize,   // driving scalar (c_j^a)

    // scalars - domain b
    pub unk_b_id: usize,   // unknown scalar (c_i^b)
    pub diff_b_id: usize,  // diffusion coefficient (G_ij^b)
    pub drv_b_id: usize,   // driving scalar (c_j^b)
}

impl OpSclItfDiffusion {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, diff_a_id: usize, drv_a_id: usize, unk_b_id: usize, diff_b_id: usize, drv_b_id: usize) -> Result<OpSclItfDiffusion, VFVChemError> {
        // adds the diffusion term to the combined control volumes of interface nodes
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        //                                          ^^^^^^^^^^^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // the interface flux leaves one half and enters the other, so it cancels

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, diff_a_id, diff_b_id)?;
        check_itf_scldom(vars, itf_id, drv_a_id, drv_b_id)?;

        // create struct
        let mut oper = OpSclItfDiffusion::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.diff_a_id = diff_a_id;
        oper.drv_a_id = drv_a_id;
        oper.unk_b_id = unk_b_id;
        oper.diff_b_id = diff_b_id;
        oper.drv_b_id = drv_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div( -G * grad(c) ) dV
        // int_surf G * grad(c) . n dA
        // sum_ngh G * grad(c) . n dA
        //
        // apply normal derivative
        // sum_ngh G_f * (c_n - c_v) / dn * dA

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for (diff_id, drv_id) in [(self.diff_a_id, self.drv_a_id), (self.diff_b_id, self.drv_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[drv_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let diff_dom = &vars.scldom[diff_id];

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
                let diff_v = diff_dom.value[vid_v];
                let (diff_x_v, diff_y_v, diff_z_v) = diff_dom.grad[vid_v];

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (sx, sy, sz) = skew_n[loc];

                    // get neighbor properties
                    let diff_n = diff_dom.value[vid_n];
                    let (diff_x_n, diff_y_n, diff_z_n) = diff_dom.grad[vid_n];

                    // get face diffusion coefficient
                    let diff_m = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean
                    let diff_f = diff_m + 0.5 * (sx * (diff_x_v + diff_x_n) + sy * (diff_y_v + diff_y_n) + sz * (diff_z_v + diff_z_n));  // skew correction
                    let diff_f = if diff_f > 0.0 { diff_f } else { diff_m };  // guard against negative diffusion

                    // load into matrix
                    let coeff = diff_f * face / dist;
                    self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, drv_id, vid_n,  coeff)?;
                    self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, drv_id, vid_v, -coeff)?;
                }

            }
        }

        Ok(())
    }
}
