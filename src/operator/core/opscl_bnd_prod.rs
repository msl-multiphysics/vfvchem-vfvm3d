use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndProduction {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub prd_id: usize,  // production coefficient (s_i)
}

impl OpSclBndProduction {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, prd_id: usize) -> Result<OpSclBndProduction, VFVChemError> {
        // adds the production term to the control volumes of boundary nodes
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i + s_i * c_i
        //                                                                            ^^^^^^^^^
        // s_i is positive if production and negative if consumption
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, prd_id)?;

        // create struct
        let mut oper = OpSclBndProduction::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.prd_id = prd_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol s * c dV
        // s_c * c_c * dV
        //
        // apply centroid skewness correction
        // s_c * c_c = s_c * c_v + s_v * grad(c_v) . skew (to first order)
        // s_c = s_v + grad(s_v) . skew
        // the first term is implicit; the second term is lagged

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];
        let prd_dom = &vars.scldom[self.prd_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let (sx, sy, sz) = tet.cvol_skew[vid_v];
            let vol = tet.cvol_vol[vid_v];

            // get vertex properties
            let prd_v = prd_dom.value[vid_v];
            let (prd_x_v, prd_y_v, prd_z_v) = prd_dom.grad[vid_v];
            let (unk_x_v, unk_y_v, unk_z_v) = unk_dom.grad[vid_v];

            // get production coefficient at centroid and skewness correction of the unknown
            let prd_c = prd_v + prd_x_v * sx + prd_y_v * sy + prd_z_v * sz;
            let unk_skew = unk_x_v * sx + unk_y_v * sy + unk_z_v * sz;

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, prd_c * vol)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -prd_v * unk_skew * vol)?;
        }

        Ok(())
    }
}
