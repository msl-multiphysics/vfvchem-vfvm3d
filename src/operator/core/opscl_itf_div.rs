use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfDivergence {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,   // velocity (vel^a); unknown
    pub pres_a_id: usize,  // pressure (p^a)
    pub den_a_id: usize,   // density (rho^a)

    // domain b
    pub vel_b_id: usize,   // velocity (vel^b); unknown
    pub pres_b_id: usize,  // pressure (p^b)
    pub den_b_id: usize,   // density (rho^b)
}

impl OpSclItfDivergence {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, pres_a_id: usize, den_a_id: usize, vel_b_id: usize, pres_b_id: usize, den_b_id: usize) -> Result<OpSclItfDivergence, VFVChemError> {
        // adds the divergence term to the combined control volumes of interface nodes
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b
        // sub-faces on the interface are inside the combined control volume and are skipped

        // input checking
        check_itf_vecdom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, pres_a_id, pres_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;

        // create struct
        let mut oper = OpSclItfDivergence::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.pres_a_id = pres_a_id;
        oper.den_a_id = den_a_id;
        oper.vel_b_id = vel_b_id;
        oper.pres_b_id = pres_b_id;
        oper.den_b_id = den_b_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_vol -d(rho)/dt dV
        // -(rho_v - rho_v_prev) * dV / dt

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.pres_b_id].dom_id];

        // iterate through both sides of the interface
        for den_id in [self.den_a_id, self.den_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[den_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let den_dom = &vars.scldom[den_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let vol = tet.cvol_vol[vid_v];

                // get vertex properties
                let den_v = den_dom.value[vid_v];
                let den_v_prev = den_dom.value_prev[vid_v];

                // load into vector
                self.add_scl_b(vars, b_vec, self.pres_b_id, vid_b_v, (den_v - den_v_prev) * vol / dt)?;
            }
        }

        Ok(())
    }
}

impl OperatorBase for OpSclItfDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div(rho * vel) dV
        // int_surf -rho * vel . n dA
        // sum_ngh -rho * vel . n dA + sum_bnd -rho * vel . n dA
        // sum_ngh -rho_f * vel_f . n * dA + sum_bnd -(w_v * rho_v * vel_v + w_j * rho_j * vel_j + w_k * rho_k * vel_k) . n * dA
        //
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.pres_b_id].dom_id];

        // iterate through both sides of the interface
        for (vel_id, den_id) in [(self.vel_a_id, self.den_a_id), (self.vel_b_id, self.den_b_id)] {
            // get objects on this side
            let dom_id = vars.vecdom[vel_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let den_dom = &vars.scldom[den_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

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
                        self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, comp, vid_n, -coeff * n_c)?;
                        self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, comp, vid_v, -coeff * n_c)?;
                    }

                    // vertex and neighbor gradient contributions; (vel_m - vel_s) terms for s = v, n
                    for vid_s in [vid_v, vid_n] {
                        for (loc_m, &vid_m) in tet.ngh_vid[vid_s].iter().enumerate() {
                            // project gradient weight onto skewness vector
                            let (w_x, w_y, w_z) = tet.grad_wgt[vid_s][loc_m];
                            let w_skew = w_x * sx + w_y * sy + w_z * sz;

                            // load into matrix
                            for (comp, n_c) in [(0, nx), (1, ny), (2, nz)] {
                                let val = -coeff * w_skew * n_c;
                                self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, comp, vid_m,  val)?;
                                self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, comp, vid_s, -val)?;
                            }
                        }
                    }
                }

                // iterate through boundary sub-faces off the interface
                for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                    let on_itf_j = itf.nds_nid.contains_key(&dom.nds_gnd[vid_j - dom.num_pts()]);
                    let on_itf_k = itf.nds_nid.contains_key(&dom.nds_gnd[vid_k - dom.num_pts()]);
                    if on_itf_j && on_itf_k { continue; }

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
                            self.add_sclvec_a(vars, a_triplet, self.pres_b_id, vid_b_v, vel_id, comp, vid_m, -coeff * n_c)?;
                        }
                    }
                }

            }
        }

        Ok(())
    }
}
