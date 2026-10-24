use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomSource {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown
    pub fce_id: usize,  // body force (fce)
}

impl OpVecDomSource {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, fce_id: usize) -> Result<OpVecDomSource, VFVChemError> {
        // adds the body force term to the Navier-Stokes equation
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                                                                                    ^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_vecdom(vars, dom_id, fce_id)?;

        // create struct
        let mut oper = OpVecDomSource::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.fce_id = fce_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol fce dV
        // fce_c * dV
        //
        // apply centroid skewness correction
        // fce_c = fce_v + grad(fce_v) . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let fce_dom = &vars.vecdom[self.fce_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let (sx, sy, sz) = tet.cvol_skew[vid_v];
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let (fce_x_v, fce_y_v, fce_z_v) = fce_dom.value[vid_v];
            let (fce_x_x_v, fce_x_y_v, fce_x_z_v, fce_y_x_v, fce_y_y_v, fce_y_z_v, fce_z_x_v, fce_z_y_v, fce_z_z_v) = fce_dom.grad[vid_v];

            // get body force at centroid
            let fce_x_c = fce_x_v + fce_x_x_v * sx + fce_x_y_v * sy + fce_x_z_v * sz;
            let fce_y_c = fce_y_v + fce_y_x_v * sx + fce_y_y_v * sy + fce_y_z_v * sz;
            let fce_z_c = fce_z_v + fce_z_x_v * sx + fce_z_y_v * sy + fce_z_z_v * sz;

            // load into vector
            self.add_vec_b(vars, b_vec, self.vel_id, 0, vid_v, -fce_x_c * vol)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 1, vid_v, -fce_y_c * vol)?;
            self.add_vec_b(vars, b_vec, self.vel_id, 2, vid_v, -fce_z_c * vol)?;
        }

        Ok(())
    }
}
