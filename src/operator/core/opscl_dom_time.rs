use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomTime {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub wgt_id: usize,  // weight coefficient (m_i)
}

impl OpSclDomTime {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, wgt_id: usize) -> Result<OpSclDomTime, VFVChemError> {
        // adds the time derivative term to scalar transport equations
        // m_i * dc_i/dt + m_i * vel . grad(c_i) = -div( -G_ij * grad(c_j) ) + R_i
        // ^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, wgt_id)?;

        // create struct
        let mut oper = OpSclDomTime::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -m * dc/dt dV
        // -m_v * (c_v - c_v_prev) * dV / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let wgt_dom = &vars.scldom[self.wgt_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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

impl OperatorBase for OpSclDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
