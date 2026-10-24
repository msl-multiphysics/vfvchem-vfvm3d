use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndTimeUnity {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
}

impl OpSclBndTimeUnity {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize) -> Result<OpSclBndTimeUnity, VFVChemError> {
        // adds the time derivative term with unit weight to the control volumes of boundary nodes
        // dc_i/dt + div(c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        // ^^^^^^^

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;

        // create struct
        let mut oper = OpSclBndTimeUnity::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -dc/dt dV
        // -(c_v - c_v_prev) * dV / dt
        //
        // lumped at the vertex

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -vol / dt)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -unk_v_prev * vol / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclBndTimeUnity {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
