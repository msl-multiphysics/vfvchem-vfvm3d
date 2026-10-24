use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfContinuity {
    // interface
    pub itf_id: usize,

    // vectors - domain a
    pub vel_a_id: usize,  // velocity (vel^a); unknown

    // vectors - domain b
    pub vel_b_id: usize,  // velocity (vel^b); unknown
}

impl OpVecItfContinuity {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, vel_b_id: usize) -> Result<OpVecItfContinuity, VFVChemError> {
        // sets a continuous interface velocity
        // vel^a = vel^b
        //
        // loads into the rows of domain a

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;

        // create struct
        let mut oper = OpVecItfContinuity::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.vel_b_id = vel_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfContinuity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_a = &vars.dom[vars.vecdom[self.vel_a_id].dom_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_a_v = dom_a.num_pts() + dom_a.nds_nid[&gnd_v];
            let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

            // load into matrix
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 0, vid_a_v, self.vel_a_id, 0, vid_a_v,  1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 0, vid_a_v, self.vel_b_id, 0, vid_b_v, -1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 1, vid_a_v, self.vel_a_id, 1, vid_a_v,  1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 1, vid_a_v, self.vel_b_id, 1, vid_b_v, -1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 2, vid_a_v, self.vel_a_id, 2, vid_a_v,  1.0)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_a_id, 2, vid_a_v, self.vel_b_id, 2, vid_b_v, -1.0)?;
        }

        Ok(())
    }
}
