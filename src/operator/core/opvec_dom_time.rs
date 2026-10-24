use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomTime {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub den_id: usize,  // density (rho)
}

impl OpVecDomTime {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize) -> Result<OpVecDomTime, VFVChemError> {
        // adds the time derivative term to the Navier-Stokes equation
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        // ^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = OpVecDomTime::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -rho * d(vel)/dt dV
        // -rho_v * (vel_v - vel_v_prev) * dV / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let vel_dom = &vars.vecdom[self.vel_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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

impl OperatorBase for OpVecDomTime {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
