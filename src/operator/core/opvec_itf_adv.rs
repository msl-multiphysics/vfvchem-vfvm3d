use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpVecItfAdvection {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub den_a_id: usize,   // density (rho^a)
    pub visc_a_id: usize,  // viscosity (mu^a); for the Peclet number
    pub flw_a_id: usize,   // flow data (domain a); for the Rhie-Chow mass flux

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub den_b_id: usize,   // density (rho^b)
    pub visc_b_id: usize,  // viscosity (mu^b); for the Peclet number
    pub flw_b_id: usize,   // flow data (domain b); for the Rhie-Chow mass flux
}

impl OpVecItfAdvection {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, den_a_id: usize, visc_a_id: usize, flw_a_id: usize, vel_b_id: usize, den_b_id: usize, visc_b_id: usize, flw_b_id: usize) -> Result<OpVecItfAdvection, VFVChemError> {
        // adds the advection term to the combined control volumes of interface nodes
        // rho * d(vel)/dt + rho * vel . grad(vel) = -grad(p) + div( mu * ( grad(vel) + grad(vel)^T ) + lam * div(vel) * I ) + fce
        //                   ^^^^^^^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // the boundary sub-faces are upwinded from the node and contribute nothing

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;
        check_itf_scldom(vars, itf_id, visc_a_id, visc_b_id)?;
        check_dom_sclflw(vars, vars.vecdom[vel_a_id].dom_id, flw_a_id)?;
        check_dom_sclflw(vars, vars.vecdom[vel_b_id].dom_id, flw_b_id)?;

        // create struct
        let mut oper = OpVecItfAdvection::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.den_a_id = den_a_id;
        oper.visc_a_id = visc_a_id;
        oper.flw_a_id = flw_a_id;
        oper.vel_b_id = vel_b_id;
        oper.den_b_id = den_b_id;
        oper.visc_b_id = visc_b_id;
        oper.flw_b_id = flw_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpVecItfAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -rho * vel . grad(vel) dV
        // int_surf -rho * (vel - vel_v) * (vel . n) dA
        // sum_ngh -rho * (vel - vel_v) * (vel . n) dA + sum_bnd -rho * (vel - vel_v) * (vel . n) dA
        //
        // introduce mass flux
        // F = rho * vel * dA
        // sum_ngh -F_f * (vel_f - vel_v) + sum_bnd -F_b * (vel_v - vel_v)
        // sum_ngh -F_f * (vel_f - vel_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (D_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (vel_n - vel_v)

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.vecdom[self.vel_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, visc_id, flw_id) in [(self.vel_a_id, self.visc_a_id, self.flw_a_id), (self.vel_b_id, self.visc_b_id, self.flw_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let visc_dom = &vars.scldom[visc_id];
            let flow_dom = &vars.sclflw[flw_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let list_n = &tet.ngh_vid[vid_v];
                let dist_n = &tet.ngh_dist[vid_v];
                let face_n = &tet.face_area[vid_v];
                let skew_n = &tet.face_skew[vid_v];
                let flux_n = &flow_dom.flux[vid_v];

                // get vertex properties
                let visc_v = visc_dom.value[vid_v];
                let (visc_x_v, visc_y_v, visc_z_v) = visc_dom.grad[vid_v];

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (sx, sy, sz) = skew_n[loc];

                    // get neighbor properties
                    let visc_n = visc_dom.value[vid_n];
                    let (visc_x_n, visc_y_n, visc_z_n) = visc_dom.grad[vid_n];

                    // get face viscosity
                    let visc_m = if visc_v + visc_n != 0.0 { 2.0 * visc_v * visc_n / (visc_v + visc_n) } else { 0.0 };  // harmonic mean
                    let visc_f = visc_m + 0.5 * (sx * (visc_x_v + visc_x_n) + sy * (visc_y_v + visc_y_n) + sz * (visc_z_v + visc_z_n));  // skew correction
                    let visc_f = if visc_f > 0.0 { visc_f } else { visc_m };  // guard against negative viscosity

                    // get power law term
                    let cond = visc_f * face / dist;
                    let flux = flux_n[loc];
                    let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                    let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                    // load into matrix
                    let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 0, vid_b_v, vel_id, 0, vid_v, -coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 1, vid_b_v, vel_id, 1, vid_v, -coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 2, vid_b_v, vel_id, 2, vid_n,  coeff)?;
                    self.add_vecvec_a(vars, a_triplet, self.vel_b_id, 2, vid_b_v, vel_id, 2, vid_v, -coeff)?;
                }

            }
        }

        Ok(())
    }
}
