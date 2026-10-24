use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecBndTime {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub den_id: usize,  // density (rho)
}

impl OpVecBndTime {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize) -> Result<OpVecBndTime, VFVChemError> {
        // adds the time derivative term to the control volumes of boundary nodes
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        // ^^^^^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;

        // create struct
        let mut oper = OpVecBndTime::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -rho * d(vel)/dt dV
        // -rho_v * (vel_v - vel_v_prev) * dV / dt

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let (u_v_prev, v_v_prev, w_v_prev) = vel_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 0, vid_v, self.vel_id, 0, vid_v, -den_v * vol / dt)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 1, vid_v, self.vel_id, 1, vid_v, -den_v * vol / dt)?;
            self.add_vecvec_a(vars, a_triplet, self.vel_id, 2, vid_v, self.vel_id, 2, vid_v, -den_v * vol / dt)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -den_v * u_v_prev * vol / dt)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -den_v * v_v_prev * vol / dt)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 2, vid_v, -den_v * w_v_prev * vol / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpVecBndTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
