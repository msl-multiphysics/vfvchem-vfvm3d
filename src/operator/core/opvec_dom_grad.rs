use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecDomGradient {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p); unknown
}

impl OpVecDomGradient {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, pres_id: usize) -> Result<OpVecDomGradient, VFVChemError> {
        // adds the pressure gradient term to the Navier-Stokes equation
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                          ^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, pres_id)?;

        // create struct
        let mut oper = OpVecDomGradient::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecDomGradient {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -grad(p) dV
        // int_surf -p n dA
        // sum_ngh -p_f * n * dA
        //
        // apply face skewness correction
        // p_f = (p_v + p_n) / 2 + (grad(p_v) + grad(p_n)) / 2 . skew
        //
        // expand the least-squares gradients; see Tetrahedron
        // sum_m -(c_x, c_y, c_z) * p_m

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];

        // get pressure gradient coefficients
        let grad_vec = &tet.grad_coeff;

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // load into matrix
            for &(vid_m, c_x, c_y, c_z) in &grad_vec[vid_v] {
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 0, vid_v, self.pres_id, vid_m, -c_x)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 1, vid_v, self.pres_id, vid_m, -c_y)?;
                self.add_vecscl_a(vars, a_triplet, self.vel_id, 2, vid_v, self.pres_id, vid_m, -c_z)?;
            }

        }

        Ok(())
    }
}
