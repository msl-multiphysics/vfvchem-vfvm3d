use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub diff_id: usize,  // diffusion coefficient (G_ij)
    pub drv_id: usize,   // driving scalar (c_j)
}

impl OpSclDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, diff_id: usize, drv_id: usize) -> Result<OpSclDomDiffusion, VFVChemError> {
        // adds the diffusion term to scalar transport equations
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        //                                          ^^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;
        check_dom_scldom(vars, dom_id, drv_id)?;

        // create struct
        let mut oper = OpSclDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.drv_id = drv_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div( -G * grad(c) ) dV
        // int_surf G * grad(c) . n dA
        // sum_ngh G * grad(c) . n dA
        //
        // apply normal derivative
        // sum_ngh G_f * (c_n - c_v) / dn * dA
        
        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let diff_dom = &vars.scldom[self.diff_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_v, -coeff)?;
            }

        }
    
        Ok(())
    }
}
