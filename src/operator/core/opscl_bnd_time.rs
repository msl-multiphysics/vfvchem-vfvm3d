use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndTime {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub wgt_id: usize,  // weight coefficient (m_i)
}

impl OpSclBndTime {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, wgt_id: usize) -> Result<OpSclBndTime, VFVChemError> {
        // adds the time derivative term to the control volumes of boundary nodes
        // m_i * dc_i/dt + m_i * vel . grad(c_i) = -div( -G_ij * grad(c_j) ) + R_i
        // ^^^^^^^^^^^^^

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, wgt_id)?;

        // create struct
        let mut oper = OpSclBndTime::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -m * dc/dt dV
        // -m_v * (c_v - c_v_prev) * dV / dt

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let wgt_dom = &vars.scldom[self.wgt_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let wgt_v = wgt_dom.value[vid_v];
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -wgt_v * vol / dt)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -wgt_v * unk_v_prev * vol / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclBndTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
