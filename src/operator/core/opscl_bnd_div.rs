use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndDivergence {
    // boundary
    pub bnd_id: usize,

    // vectors
    pub vel_id: usize,  // velocity (vel); unknown

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,  // density (rho)
}

impl OpSclBndDivergence {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, pres_id: usize, den_id: usize) -> Result<OpSclBndDivergence, VFVChemError> {
        // adds the divergence term to the control volumes of boundary nodes
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^

        // input checking
        check_bnd_vecdom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, pres_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;

        // create struct
        let mut oper = OpSclBndDivergence::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.pres_id = pres_id;
        oper.den_id = den_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -d(rho)/dt dV
        // -(rho_v - rho_v_prev) * dV / dt

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
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
            let den_v_prev = den_dom.value_prev[vid_v];

            // load into vector
            self.add_scl_b(vars, b_vec, self.pres_id, vid_v, (den_v - den_v_prev) * vol / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclBndDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div(rho * vel) dV
        // int_surf -rho * vel . n dA
        // sum_ngh -rho * vel . n dA + sum_bnd -rho * vel . n dA
        // sum_ngh -rho_f * vel_f . n * dA + sum_bnd -(w_v * rho_v * vel_v + w_j * rho_j * vel_j + w_k * rho_k * vel_k) . n * dA
        //
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &tet.ngh_vid[vid_v];
            let face_n = &tet.face_area[vid_v];
            let norm_n = &tet.face_norm[vid_v];
            let skew_n = &tet.face_skew[vid_v];
            let list_b = &tet.bnd_vid[vid_v];
            let area_b = &tet.bnd_area[vid_v];
            let norm_b = &tet.bnd_norm[vid_v];
            let wgt_b = &tet.bnd_wgt[vid_v];

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

            // iterate through boundary sub-faces
            for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                // get geometry
                let area = area_b[blc];
                let (nx, ny, nz) = norm_b[blc];
                let (w_v, w_j, w_k) = wgt_b[blc];

                // get neighbor properties
                let den_j = den_dom.value[vid_j];
                let den_k = den_dom.value[vid_k];

                // load into matrix
                for (vid_m, coeff) in [(vid_v, w_v * den_v * area), (vid_j, w_j * den_j * area), (vid_k, w_k * den_k * area)] {
                    for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                        self.add_sclvec_a(vars, a_triplet, self.pres_id, vid_v, self.vel_id, comp, vid_m, -coeff * n_c)?;
                    }
                }
            }

        }

        Ok(())
    }
}
