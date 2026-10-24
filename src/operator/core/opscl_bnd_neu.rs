use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndNeumann {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub flux_id: usize,  // flux (N_set)
}

impl OpSclBndNeumann {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, flux_id: usize) -> Result<OpSclBndNeumann, VFVChemError> {
        // sets the outward boundary flux
        // N_i . n = N_set

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_sclbnd(vars, bnd_id, flux_id)?;
        
        // create struct
        let mut oper = OpSclBndNeumann::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.flux_id = flux_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndNeumann {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_bnd N dA
        // sum_bnd (w_v * N_v + w_j * N_j + w_k * N_k) * dA
        
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let tet = &vars.tet[bnd.dom_id];
        let flux_bnd = &vars.sclbnd[self.flux_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_b = &tet.bnd_vid[vid_v];
            let area_b = &tet.bnd_area[vid_v];
            let wgt_b = &tet.bnd_wgt[vid_v];

            // get vertex properties
            let flux_v = flux_bnd.value[idx_v];

            // iterate through boundary sub-faces
            for (blc, &(vid_j, vid_k)) in list_b.iter().enumerate() {
                // get indices
                let gnd_j = dom.nds_gnd[vid_j - dom.num_pts()];
                let gnd_k = dom.nds_gnd[vid_k - dom.num_pts()];

                // get geometry
                let area = area_b[blc];
                let (w_v, w_j, w_k) = wgt_b[blc];

                // get neighbor properties
                // use vertex flux if neighbor is on another boundary
                let flux_j = match bnd.nds_idx.get(&gnd_j) {
                    Some(&idx_j) => flux_bnd.value[idx_j],
                    None => flux_v,
                };
                let flux_k = match bnd.nds_idx.get(&gnd_k) {
                    Some(&idx_k) => flux_bnd.value[idx_k],
                    None => flux_v,
                };

                // load into matrix
                self.add_scl_b(vars, b_vec, self.unk_id, vid_v, area * (w_v * flux_v + w_j * flux_j + w_k * flux_k))?;
            }
        }

        Ok(())
    }
}
