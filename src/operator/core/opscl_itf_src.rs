use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfSource {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,  // unknown scalar (c_i^a)
    pub src_a_id: usize,  // source (R_i^a)

    // scalars - domain b
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
    pub src_b_id: usize,  // source (R_i^b)
}

impl OpSclItfSource {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, src_a_id: usize, unk_b_id: usize, src_b_id: usize) -> Result<OpSclItfSource, VFVChemError> {
        // adds the source term to the combined control volumes of interface nodes
        // d(m_i * c_i)/dt + div(m_i * c_i * vel) = -div( -G_ij * grad(c_j) ) + R_i
        //                                                                      ^^^
        // each domain contributes its half of the control volume; loads into the rows of domain b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, src_a_id, src_b_id)?;

        // create struct
        let mut oper = OpSclItfSource::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.src_a_id = src_a_id;
        oper.unk_b_id = unk_b_id;
        oper.src_b_id = src_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_vol R dV
        // R_c * dV
        //
        // apply centroid skewness correction
        // R_c = R_v + grad(R_v) . skew

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for src_id in [self.src_a_id, self.src_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[src_id].dom_id;
            let dom = &vars.dom[dom_id];
            let tet = &vars.tet[dom_id];
            let src_dom = &vars.scldom[src_id];

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
                let src_v = src_dom.value[vid_v];
                let (src_x_v, src_y_v, src_z_v) = src_dom.grad[vid_v];

                // get source at centroid
                let src_c = src_v + src_x_v * sx + src_y_v * sy + src_z_v * sz;

                // load into vector
                self.add_scl_b(vars, b_vec, self.unk_b_id, vid_b_v, -src_c * vol)?;
            }
        }

        Ok(())
    }
}
