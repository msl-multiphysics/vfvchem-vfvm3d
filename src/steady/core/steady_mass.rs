use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::steady::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct SteadyMass {
    // internal data
    // indexed by component, then by domain
    pub num_comp: usize,                                       // number of components
    pub itr_conc: Vec<HashMap<usize, usize>>,                  // concentration
    pub itr_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // diffusion coefficient; component j -> D_ij
    pub itr_mprd: Vec<HashMap<usize, usize>>,                  // mass production coefficient
    pub itr_msrc: Vec<HashMap<usize, usize>>,                  // mass source

    // boundary data
    pub conc_bnd: Vec<Vec<usize>>,              // boundary with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>,  // concentration
    pub mflx_bnd: Vec<Vec<usize>>,              // boundary with diffusive mass flux BC
    pub mflx_mflx: Vec<HashMap<usize, usize>>,  // diffusive mass flux

    // interface data
    pub mcnt_itf: Vec<Vec<usize>>,              // interface with continuity IC

    // operators
    pub oper_itr_mass: Vec<(Vec<OpSclDomDiffusion>, OpSclDomProduction, OpSclDomSource)>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<(Vec<OpSclBndDiffusion>, OpSclBndProduction, OpSclBndSource, OpSclBndNeumann)>,
    pub oper_itf_mcnt: Vec<(OpSclItfContinuity, Vec<OpSclItfDiffusion>, OpSclItfProduction, OpSclItfSource)>,
}

impl SteadyMass {
    pub fn new(num_comp: usize) -> Result<SteadyMass, VFVChemError> {
        // input checking
        if num_comp == 0 { return Err(VFVChemError::MassNonPositiveNumComponents { num_comp }) }

        // create struct with data for each component
        let mut phys = SteadyMass::default();
        phys.num_comp = num_comp;
        phys.itr_conc = vec![HashMap::new(); num_comp];
        phys.itr_diff = vec![HashMap::new(); num_comp];
        phys.itr_mprd = vec![HashMap::new(); num_comp];
        phys.itr_msrc = vec![HashMap::new(); num_comp];
        phys.conc_bnd = vec![Vec::new(); num_comp];
        phys.conc_conc = vec![HashMap::new(); num_comp];
        phys.mflx_bnd = vec![Vec::new(); num_comp];
        phys.mflx_mflx = vec![HashMap::new(); num_comp];
        phys.mcnt_itf = vec![Vec::new(); num_comp];
        Ok(phys)
    }

    pub fn add_mass_dom(&mut self, comp_id: usize, dom_id: usize, conc_id: usize, diff_ids: HashMap<usize, usize>, mprd_id: usize, msrc_id: usize) -> Result<(), VFVChemError> {
        for &comp in diff_ids.keys().chain([&comp_id]) {
            if comp >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp, num_comp: self.num_comp }) }
        }
        self.itr_conc[comp_id].insert(dom_id, conc_id);
        self.itr_diff[comp_id].insert(dom_id, diff_ids);
        self.itr_mprd[comp_id].insert(dom_id, mprd_id);
        self.itr_msrc[comp_id].insert(dom_id, msrc_id);
        Ok(())
    }

    pub fn add_conc_bnd(&mut self, comp_id: usize, bnd_id: usize, conc_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.conc_bnd[comp_id].push(bnd_id);
        self.conc_conc[comp_id].insert(bnd_id, conc_id);
        Ok(())
    }

    pub fn add_mflx_bnd(&mut self, comp_id: usize, bnd_id: usize, mflx_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mflx_bnd[comp_id].push(bnd_id);
        self.mflx_mflx[comp_id].insert(bnd_id, mflx_id);
        Ok(())
    }

    pub fn add_mcnt_itf(&mut self, comp_id: usize, itf_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mcnt_itf[comp_id].push(itf_id);
        Ok(())
    }

}

impl SteadyBase for SteadyMass {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown scalars
        let mut xid = 0;  // index in solution vector
        for itr_conc in &self.itr_conc {
        for (&dom_id, &conc_id) in itr_conc {
            vars.scldom[conc_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }}

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr_mass.clear();
        self.oper_bnd_conc.clear();
        self.oper_bnd_mflx.clear();
        self.oper_itf_mcnt.clear();

        for comp_i in 0..self.num_comp {
            // internal
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let mprd_id = self.itr_mprd[comp_i][&dom_id];
                let msrc_id = self.itr_msrc[comp_i][&dom_id];
                let mut oper_diff = Vec::new();
                for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclDomDiffusion::new(vars, dom_id, conc_id, diff_id, conc_j_id)?);
                }
                let oper_prd = OpSclDomProduction::new(vars, dom_id, conc_id, mprd_id)?;
                let oper_src = OpSclDomSource::new(vars, dom_id, conc_id, msrc_id)?;
                self.oper_itr_mass.push((oper_diff, oper_prd, oper_src));
            }

            // boundary
            for &bnd_id in &self.conc_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_conc_id = self.conc_conc[comp_i][&bnd_id];
                let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_conc_id, bnd_conc_id)?;
                self.oper_bnd_conc.push(oper_dir);
            }
            for &bnd_id in &self.mflx_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let dom_mprd_id = self.itr_mprd[comp_i][&dom_id];
                let dom_msrc_id = self.itr_msrc[comp_i][&dom_id];
                let bnd_mflx_id = self.mflx_mflx[comp_i][&bnd_id];
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let dom_conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclBndDiffusion::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_conc_j_id)?);
                }
                let oper_prd = OpSclBndProduction::new(vars, bnd_id, dom_conc_id, dom_mprd_id)?;
                let oper_src = OpSclBndSource::new(vars, bnd_id, dom_conc_id, dom_msrc_id)?;
                let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_conc_id, bnd_mflx_id)?;
                self.oper_bnd_mflx.push((oper_diff, oper_prd, oper_src, oper_neu));
            }

            // interface
            for &itf_id in &self.mcnt_itf[comp_i] {
                let dom_a_id = vars.itf[itf_id].dom_a_id;
                let dom_b_id = vars.itf[itf_id].dom_b_id;
                let dom_conc_a_id = self.itr_conc[comp_i][&dom_a_id];
                let dom_conc_b_id = self.itr_conc[comp_i][&dom_b_id];
                let dom_mprd_a_id = self.itr_mprd[comp_i][&dom_a_id];
                let dom_mprd_b_id = self.itr_mprd[comp_i][&dom_b_id];
                let dom_msrc_a_id = self.itr_msrc[comp_i][&dom_a_id];
                let dom_msrc_b_id = self.itr_msrc[comp_i][&dom_b_id];
                let oper_cont = OpSclItfContinuity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id)?;
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_a_id) in &self.itr_diff[comp_i][&dom_a_id] {
                    let dom_diff_b_id = self.itr_diff[comp_i][&dom_b_id][&comp_j];
                    let dom_conc_j_a_id = self.itr_conc[comp_j][&dom_a_id];
                    let dom_conc_j_b_id = self.itr_conc[comp_j][&dom_b_id];
                    oper_diff.push(OpSclItfDiffusion::new(vars, itf_id, dom_conc_a_id, dom_diff_a_id, dom_conc_j_a_id, dom_conc_b_id, dom_diff_b_id, dom_conc_j_b_id)?);
                }
                let oper_prd = OpSclItfProduction::new(vars, itf_id, dom_conc_a_id, dom_mprd_a_id, dom_conc_b_id, dom_mprd_b_id)?;
                let oper_src = OpSclItfSource::new(vars, itf_id, dom_conc_a_id, dom_msrc_a_id, dom_conc_b_id, dom_msrc_b_id)?;
                self.oper_itf_mcnt.push((oper_cont, oper_diff, oper_prd, oper_src));
            }
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // internal
        for (oper_diff, oper_prd, oper_src) in &self.oper_itr_mass {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        // boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_prd, oper_src, oper_neu) in &self.oper_bnd_mflx {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }

        // interface
        for (oper_cont, oper_diff, oper_prd, oper_src) in &self.oper_itf_mcnt {
            oper_cont.apply(vars, a_triplet, b_vec)?;
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply_initial(vars);
        }
    }

    fn update_flow(&self, _vars: &mut Variables) {
        // no flow data
    }
}
