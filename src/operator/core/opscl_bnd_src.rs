use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndSource {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
}

impl OpSclBndSource {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, src_id: usize) -> Result<OpSclBndSource, VFVChemError> {
        // adds the source term to the control volumes of boundary nodes
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        //                                                                      ^^^

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, src_id)?;

        // create struct
        let mut oper = OpSclBndSource::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        
        Ok(oper)
    }
}

impl OperatorBase for OpSclBndSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol R dV
        // R_c * dV
        //
        // apply centroid skewness correction
        // R_c = R_v + grad(R_v) . skew
        
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let src_dom = &vars.scldom[self.src_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let (sx, sy, sz) = tet.cvol_skew[vid_v];
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let src_v = src_dom.value[vid_v];
            let (src_x_v, src_y_v, src_z_v) = src_dom.grad[vid_v];

            // get source at centroid
            let src_c = src_v + src_x_v * sx + src_y_v * sy + src_z_v * sz;
            
            // load into matrix
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_c * vol)?;
        }
    
        Ok(())
    }
}
