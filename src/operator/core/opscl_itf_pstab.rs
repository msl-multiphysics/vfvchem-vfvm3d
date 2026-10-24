use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfPressureStabilizer {
    // interface
    pub itf_id: usize,

    // domain a
    pub pres_a_id: usize,  // pressure (p^a); unknown
    pub den_a_id: usize,   // density (rho^a)
    pub fce_a_id: usize,   // body force (fce^a)
    pub flw_a_id: usize,   // flow data (domain a); for the Rhie-Chow coefficient

    // domain b
    pub pres_b_id: usize,  // pressure (p^b); unknown
    pub den_b_id: usize,   // density (rho^b)
    pub fce_b_id: usize,   // body force (fce^b)
    pub flw_b_id: usize,   // flow data (domain b); for the Rhie-Chow coefficient
}

impl OpSclItfPressureStabilizer {
    pub fn new(vars: &Variables, itf_id: usize, pres_a_id: usize, den_a_id: usize, fce_a_id: usize, flw_a_id: usize, pres_b_id: usize, den_b_id: usize, fce_b_id: usize, flw_b_id: usize) -> Result<OpSclItfPressureStabilizer, VFVChemError> {
        // adds the Rhie-Chow correction to the face velocities of the combined control volumes of interface nodes
        // d(rho)/dt = -div(rho * vel)
        //             ^^^^^^^^^^^^^^^
        // vel_f . n = (vel_v + vel_n) / 2 . n - d_f * ((p_n - p_v) / dn - fce_f . n) + (d_v * (grad(p_v) - fce_v) + d_n * (grad(p_n) - fce_n)) / 2 . n
        // see OpSclDomPressureStabilizer
        // each domain contributes the faces of its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_scldom(vars, itf_id, pres_a_id, pres_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;
        check_itf_vecdom(vars, itf_id, fce_a_id, fce_b_id)?;
        check_dom_sclflw(vars, vars.scldom[pres_a_id].dom_id, flw_a_id)?;
        check_dom_sclflw(vars, vars.scldom[pres_b_id].dom_id, flw_b_id)?;

        // create struct
        let mut oper = OpSclItfPressureStabilizer::default();
        oper.itf_id = itf_id;
        oper.pres_a_id = pres_a_id;
        oper.den_a_id = den_a_id;
        oper.fce_a_id = fce_a_id;
        oper.flw_a_id = flw_a_id;
        oper.pres_b_id = pres_b_id;
        oper.den_b_id = den_b_id;
        oper.fce_b_id = fce_b_id;
        oper.flw_b_id = flw_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfPressureStabilizer {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol -div(rho * vel) dV
        // int_surf -rho * vel . n dA
        // sum_ngh -rho_f * vel_f . n * dA
        //
        // apply Rhie-Chow correction
        // sum_ngh rho_f * (d_f * (p_n - p_v) / dn - (d_v * grad(p_v) + d_n * grad(p_n)) / 2 . n) * dA
        // sum_ngh -rho_f * (d_f * fce_f - (d_v * fce_v + d_n * fce_n) / 2) . n * dA
        //
        // apply face skewness correction
        // rho_f = (rho_v + rho_n) / 2 + (grad(rho_v) + grad(rho_n)) / 2 . skew
        //
        // use the momentum pressure gradient and body force
        // grad(p_v) = sum_m (c_x, c_y, c_z) * p_m / V_v; see Tetrahedron
        // fce_v = fce at the control volume centroid; see OpVecItfSource
        // fce_f = (fce_v + fce_n) / 2 at the edge midpoint, consistent with (p_n - p_v) / dn

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.pres_b_id].dom_id];

        // iterate through both sides of the interface
        for (pres_id, den_id, fce_id, flw_id) in [(self.pres_a_id, self.den_a_id, self.fce_a_id, self.flw_a_id), (self.pres_b_id, self.den_b_id, self.fce_b_id, self.flw_b_id)] {
            // get objects on this side
            let dom_id = vars.scldom[pres_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let den_dom = &vars.scldom[den_id];
            let fce_dom = &vars.vecdom[fce_id];
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
                let norm_n = &tet.face_norm[vid_v];
                let skew_n = &tet.face_skew[vid_v];
                let (cx_v, cy_v, cz_v) = tet.cvol_skew[vid_v];
                let vol_v = tet.cvol_vol[vid_v];

                // get vertex properties
                let den_v = den_dom.value[vid_v];
                let (den_x_v, den_y_v, den_z_v) = den_dom.grad[vid_v];
                let (fce_x_v, fce_y_v, fce_z_v) = fce_dom.value[vid_v];
                let (fce_x_x_v, fce_x_y_v, fce_x_z_v, fce_y_x_v, fce_y_y_v, fce_y_z_v, fce_z_x_v, fce_z_y_v, fce_z_z_v) = fce_dom.grad[vid_v];
                let fce_x_c_v = fce_x_v + fce_x_x_v * cx_v + fce_x_y_v * cy_v + fce_x_z_v * cz_v;
                let fce_y_c_v = fce_y_v + fce_y_x_v * cx_v + fce_y_y_v * cy_v + fce_y_z_v * cz_v;
                let fce_z_c_v = fce_z_v + fce_z_x_v * cx_v + fce_z_y_v * cy_v + fce_z_z_v * cz_v;
                let rc_v = flow_dom.rc_coeff[vid_v];

                // iterate through faces
                for (loc, &vid_n) in list_n.iter().enumerate() {
                    // get geometry
                    let face = face_n[loc];
                    let dist = dist_n[loc];
                    let (nx, ny, nz) = norm_n[loc];
                    let (sx, sy, sz) = skew_n[loc];
                    let (cx_n, cy_n, cz_n) = tet.cvol_skew[vid_n];
                    let vol_n = tet.cvol_vol[vid_n];

                    // get neighbor properties
                    let den_n = den_dom.value[vid_n];
                    let (den_x_n, den_y_n, den_z_n) = den_dom.grad[vid_n];
                    let (fce_x_n, fce_y_n, fce_z_n) = fce_dom.value[vid_n];
                    let (fce_x_x_n, fce_x_y_n, fce_x_z_n, fce_y_x_n, fce_y_y_n, fce_y_z_n, fce_z_x_n, fce_z_y_n, fce_z_z_n) = fce_dom.grad[vid_n];
                    let fce_x_c_n = fce_x_n + fce_x_x_n * cx_n + fce_x_y_n * cy_n + fce_x_z_n * cz_n;
                    let fce_y_c_n = fce_y_n + fce_y_x_n * cx_n + fce_y_y_n * cy_n + fce_y_z_n * cz_n;
                    let fce_z_c_n = fce_z_n + fce_z_x_n * cx_n + fce_z_y_n * cy_n + fce_z_z_n * cz_n;
                    let rc_n = flow_dom.rc_coeff[vid_n];

                    // get face density
                    let den_m = 0.5 * (den_v + den_n);  // arithmetic mean
                    let den_f = den_m + 0.5 * (sx * (den_x_v + den_x_n) + sy * (den_y_v + den_y_n) + sz * (den_z_v + den_z_n));  // skewness correction

                    // load into matrix
                    let coeff = den_f * 0.5 * (rc_v + rc_n) * face;
                    self.add_sclscl_a(vars, a_triplet, self.pres_b_id, vid_b_v, pres_id, vid_n,  coeff / dist)?;
                    self.add_sclscl_a(vars, a_triplet, self.pres_b_id, vid_b_v, pres_id, vid_v, -coeff / dist)?;

                    // vertex gradient contribution
                    let coeff_v = den_f * 0.5 * rc_v / vol_v * face;
                    for &(vid_m, c_x, c_y, c_z) in &tet.grad_coeff[vid_v] {
                        self.add_sclscl_a(vars, a_triplet, self.pres_b_id, vid_b_v, pres_id, vid_m, -coeff_v * (c_x * nx + c_y * ny + c_z * nz))?;
                    }

                    // neighbor gradient contribution
                    let coeff_n = den_f * 0.5 * rc_n / vol_n * face;
                    for &(vid_m, c_x, c_y, c_z) in &tet.grad_coeff[vid_n] {
                        self.add_sclscl_a(vars, a_triplet, self.pres_b_id, vid_b_v, pres_id, vid_m, -coeff_n * (c_x * nx + c_y * ny + c_z * nz))?;
                    }

                    // get body force correction
                    let fce_x_rc = 0.5 * (rc_v + rc_n) * 0.5 * (fce_x_v + fce_x_n) - 0.5 * (rc_v * fce_x_c_v + rc_n * fce_x_c_n);
                    let fce_y_rc = 0.5 * (rc_v + rc_n) * 0.5 * (fce_y_v + fce_y_n) - 0.5 * (rc_v * fce_y_c_v + rc_n * fce_y_c_n);
                    let fce_z_rc = 0.5 * (rc_v + rc_n) * 0.5 * (fce_z_v + fce_z_n) - 0.5 * (rc_v * fce_z_c_v + rc_n * fce_z_c_n);

                    // load into vector
                    self.add_scl_b(vars, b_vec, self.pres_b_id, vid_b_v, den_f * (fce_x_rc * nx + fce_y_rc * ny + fce_z_rc * nz) * face)?;
                }

            }
        }

        Ok(())
    }
}
