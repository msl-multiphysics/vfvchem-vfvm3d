use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfContinuity {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,  // unknown scalar (c_i^a)

    // scalars - domain b
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
}

impl OpSclItfContinuity {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize) -> Result<OpSclItfContinuity, VFVChemError> {
        // sets a continuous interface value
        // c_i^a = c_i^b
        //
        // loads into the rows of domain a

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;

        // create struct
        let mut oper = OpSclItfContinuity::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.unk_b_id = unk_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfContinuity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_a = &vars.dom[vars.scldom[self.unk_a_id].dom_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_a_v = dom_a.num_pts() + dom_a.nds_nid[&gnd_v];
            let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.unk_a_id, vid_a_v, self.unk_a_id, vid_a_v,  1.0)?;
            self.add_sclscl_a(vars, a_triplet, self.unk_a_id, vid_a_v, self.unk_b_id, vid_b_v, -1.0)?;
        }

        Ok(())
    }
}
