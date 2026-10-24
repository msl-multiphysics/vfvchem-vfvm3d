use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
}

impl OpSclDomSource {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, src_id: usize) -> Result<OpSclDomSource, VFVChemError> {
        // adds the source term to scalar transport equations
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        //                                                                      ^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, src_id)?;

        // create struct
        let mut oper = OpSclDomSource::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;
        
        Ok(oper)
    }
}

impl OperatorBase for OpSclDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol R dV
        // R_c * dV
        //
        // apply centroid skewness correction
        // R_c = R_v + grad(R_v) . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let src_dom = &vars.scldom[self.src_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
