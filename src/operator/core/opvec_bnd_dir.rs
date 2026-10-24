use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndDirichlet {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown
    pub val_id: usize,  // value (vel_set)
}

impl OpVecBndDirichlet {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, val_id: usize) -> Result<OpVecBndDirichlet, VFVChemError> {
        // sets the boundary value
        // vel = vel_set

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_vecbnd(vars, bnd_id, val_id)?;
        
        // create struct
        let mut oper = OpVecBndDirichlet::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.val_id = val_id;
        
        Ok(oper)
    }

    pub fn apply_initial(&self, vars: &mut Variables) {
        // loads the boundary value into the unknown
        // vel = vel_set

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
            vars.vecdom[self.vel_id].value[vid_v] = vars.vecbnd[self.val_id].value[idx_v];
        }
    }
}

impl OperatorBase for OpVecBndDirichlet {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let val_bnd = &vars.vecbnd[self.val_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;
            
            // get value
            let (u, v, w) = val_bnd.value[idx_v];

            // load into matrix
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, 1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, 1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 2, vid_v, self.vel_id, 2, vid_v, 1.0)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, u)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, v)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 2, vid_v, w)?;
        }

        Ok(())
    }
}
