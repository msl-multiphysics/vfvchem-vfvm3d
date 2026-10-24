use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfSource {
    // interface
    pub itf_id: usize,

    // vectors - domain a
    pub vel_a_id: usize,  // velocity (vel^a); unknown
    pub fce_a_id: usize,  // body force (fce^a)

    // vectors - domain b
    pub vel_b_id: usize,  // velocity (vel^b); unknown
    pub fce_b_id: usize,  // body force (fce^b)
}

impl OpVecItfSource {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, fce_a_id: usize, vel_b_id: usize, fce_b_id: usize) -> Result<OpVecItfSource, VFVChemError> {
        // adds the body force term to the combined control volumes of interface nodes
        // d(rho * vel)/dt + div(rho * vel * vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                                                                                                                    ^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_vecdom(vars, itf_id, fce_a_id, fce_b_id)?;

        // create struct
        let mut oper = OpVecItfSource::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.fce_a_id = fce_a_id;
        oper.vel_b_id = vel_b_id;
        oper.fce_b_id = fce_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol fce dV
        // fce_c * dV
        //
        // apply centroid skewness correction
        // fce_c = fce_v + grad(fce_v) . skew

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for fce_id in [self.fce_a_id, self.fce_b_id] {
            // get objects on this side
            let dom_id = vars.vecdom[fce_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let fce_dom = &vars.vecdom[fce_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

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
                self.add_vec_b(vars, b_vec, self.vel_b_id, 0, vid_b_v, -fce_x_c * vol)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 1, vid_b_v, -fce_y_c * vol)?;
                self.add_vec_b(vars, b_vec, self.vel_b_id, 2, vid_b_v, -fce_z_c * vol)?;
            }
        }

        Ok(())
    }
}
