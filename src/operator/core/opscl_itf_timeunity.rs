use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfTimeUnity {
    // interface
    pub itf_id: usize,

    // scalars
    pub unk_a_id: usize,  // unknown scalar (c_i^a)
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
}

impl OpSclItfTimeUnity {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, unk_b_id: usize) -> Result<OpSclItfTimeUnity, VFVChemError> {
        // adds the time derivative term with unit weight to the combined control volumes of interface nodes
        // dc_i/dt + div(c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        // ^^^^^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;

        // create struct
        let mut oper = OpSclItfTimeUnity::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.unk_b_id = unk_b_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -dc/dt dV
        // -(c_v - c_v_prev) * dV / dt
        //
        // lumped at the vertex

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for unk_id in [self.unk_a_id, self.unk_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[unk_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let unk_dom = &vars.scldom[unk_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let vol = tet.cvol_vol[vid_v];

                // get vertex properties
                let unk_v_prev = unk_dom.value_prev[vid_v];

                // load into matrix and vector
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, unk_id, vid_v, -vol / dt)?;
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -unk_v_prev * vol / dt)?;
            }
        }

        Ok(())
    }
}

impl OperatorBase for OpSclItfTimeUnity {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
