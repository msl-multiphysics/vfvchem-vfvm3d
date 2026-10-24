use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndDirichlet {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub val_id: usize,  // value (c_set)
}

impl OpSclBndDirichlet {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, val_id: usize) -> Result<OpSclBndDirichlet, VFVChemError> {
        // sets the boundary value
        // c_i = c_set

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_sclbnd(vars, bnd_id, val_id)?;
        
        // create struct
        let mut oper = OpSclBndDirichlet::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.val_id = val_id;
        
        Ok(oper)
    }

    pub fn apply_initial(&self, vars: &mut Variables) {
        // loads the boundary value into the unknown
        // c_i = c_set

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // load into unknown
            vars.scldom[self.unk_id].value[vid_v] = vars.sclbnd[self.val_id].value[idx_v];
        }
    }
}

impl OperatorBase for OpSclBndDirichlet {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let val_bnd = &vars.sclbnd[self.val_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;
            
            // get value
            let value = val_bnd.value[idx_v];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, 1.0)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, value)?;
        }

        Ok(())
    }
}
