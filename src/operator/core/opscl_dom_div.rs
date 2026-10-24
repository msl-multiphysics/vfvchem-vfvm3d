use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,  // density (rho)
}

impl OpSclDomDivergence {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize) -> Result<OpSclDomDivergence, VFVChemError> {
        // adds the divergence term to the continuity equation
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^

        // input checking
        check_dom_vecdom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, pres_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = OpSclDomDivergence::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -d(rho)/dt dV
        // -(rho_v - rho_v_prev) * dV / dt

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let den_v_prev = den_dom.value_prev[vid_v];

            // load into vector
            self.add_scl_b(vars, b_vec, self.pres_id, vid_v, (den_v - den_v_prev) * vol / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div(rho * vel) dV
        // int_surf -rho * vel . n dA
        // sum_ngh -rho * vel . n dA
        // sum_ngh -rho_f * vel_f . n * dA
        // 
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let dom = &vars.dom[self.dom_id];
        let tet = &vars.tet[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &tet.ngh_vid[vid_v];
            let face_n = &tet.face_area[vid_v];
            let norm_n = &tet.face_norm[vid_v];
            let skew_n = &tet.face_skew[vid_v];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let (den_x_v, den_y_v, den_z_v) = den_dom.grad[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let face = face_n[loc];
                let (nx, ny, nz) = norm_n[loc];
                let (sx, sy, sz) = skew_n[loc];

                // get neighbor properties
                let den_n = den_dom.value[vid_n];
                let (den_x_n, den_y_n, den_z_n) = den_dom.grad[vid_n];

                // get face density
                let den_m = 0.5 * (den_v + den_n);  // arithmetic mean
                let den_f = den_m + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n) + sz * (den_z_v + den_z_n));  // skewness correction

                // load into matrix
                let coeff = 0.5 * den_f * face;
                for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_n, -coeff * n_c)?;
                    self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_v, -coeff * n_c)?;
                }

                // vertex gradient contribution
                let list_vm = &tet.ngh_vid[vid_v];
                let wgt_vm = &tet.grad_wgt[vid_v];
                for (loc_m, &vid_m) in list_vm.iter().enumerate() {
                    // project gradient weight onto skewness vector
                    let (w_x, w_y, w_z) = wgt_vm[loc_m];
                    let w_skew = w_x * sx + w_y * sy + w_z * sz;

                    // (vel_m - vel_v) terms
                    for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                        let val = -coeff * w_skew * n_c;
                        self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_m,  val)?;
                        self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_v, -val)?;
                    }
                }

                // neighbor gradient contribution
                let list_nm = &tet.ngh_vid[vid_n];
                let wgt_nm = &tet.grad_wgt[vid_n];
                for (loc_m, &vid_m) in list_nm.iter().enumerate() {
                    // project gradient weight onto skewness vector
                    let (w_x, w_y, w_z) = wgt_nm[loc_m];
                    let w_skew = w_x * sx + w_y * sy + w_z * sz;

                    // (vel_m - vel_n) terms
                    for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                        let val = -coeff * w_skew * n_c;
                        self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_m,  val)?;
                        self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_n, -val)?;
                    }
                }
            }

        }

        Ok(())
    }
}
