use crate::base::prelude::*;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientHeat {
    // internal data
    pub itr_dom: Vec<usize>,              // domain
    pub itr_temp: HashMap<usize, usize>,  // temperature
    pub itr_vlcp: HashMap<usize, usize>,  // volumetric heat capacity
    pub itr_cond: HashMap<usize, usize>,  // thermal conductivity
    pub itr_hprd: HashMap<usize, usize>,  // heat production coefficient
    pub itr_hsrc: HashMap<usize, usize>,  // heat source

    // boundary data
    pub temp_bnd: Vec<usize>,              // boundary with temperature BC
    pub temp_temp: HashMap<usize, usize>,  // temperature
    pub hflx_bnd: Vec<usize>,              // boundary with heat flux BC
    pub hflx_hflx: HashMap<usize, usize>,  // heat flux

    // interface data
    pub hcnt_itf: Vec<usize>,              // interface with continuity IC

    // operators
    pub oper_itr: Vec<(OpSclDomTime, OpSclDomDiffusion, OpSclDomProduction, OpSclDomSource)>,
    pub oper_bnd_temp: Vec<OpSclBndDirichlet>,
    pub oper_bnd_hflx: Vec<(OpSclBndTime, OpSclBndDiffusion, OpSclBndProduction, OpSclBndSource, OpSclBndNeumann)>,
    pub oper_itf_hcnt: Vec<(OpSclItfTime, OpSclItfContinuity, OpSclItfDiffusion, OpSclItfProduction, OpSclItfSource)>,
}

impl TransientHeat {
    pub fn new() -> Result<TransientHeat, VFVChemError> {
        let phys = TransientHeat::default();
        Ok(phys)
    }

    pub fn add_heat_dom(&mut self, dom_id: usize, temp_id: usize, vlcp_id: usize, cond_id: usize, hprd_id: usize, hsrc_id: usize) -> Result<(), VFVChemError> {
        self.itr_dom.push(dom_id);
        self.itr_temp.insert(dom_id, temp_id);
        self.itr_vlcp.insert(dom_id, vlcp_id);
        self.itr_cond.insert(dom_id, cond_id);
        self.itr_hprd.insert(dom_id, hprd_id);
        self.itr_hsrc.insert(dom_id, hsrc_id);
        Ok(())
    }

    pub fn add_temp_bnd(&mut self, bnd_id: usize, temp_id: usize) -> Result<(), VFVChemError> {
        self.temp_bnd.push(bnd_id);
        self.temp_temp.insert(bnd_id, temp_id);
        Ok(())
    }

    pub fn add_hflx_bnd(&mut self, bnd_id: usize, hflx_id: usize) -> Result<(), VFVChemError> {
        self.hflx_bnd.push(bnd_id);
        self.hflx_hflx.insert(bnd_id, hflx_id);
        Ok(())
    }

    pub fn add_hcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.hcnt_itf.push(itf_id);
        Ok(())
    }

}

impl TransientBase for TransientHeat {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &temp_id) in &self.itr_temp {
            vars.scldom[temp_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr.clear();
        self.oper_bnd_temp.clear();
        self.oper_bnd_hflx.clear();
        self.oper_itf_hcnt.clear();

        // internal
        for &dom_id in &self.itr_dom {
            let temp_id = self.itr_temp[&dom_id];
            let vlcp_id = self.itr_vlcp[&dom_id];
            let cond_id = self.itr_cond[&dom_id];
            let hprd_id = self.itr_hprd[&dom_id];
            let hsrc_id = self.itr_hsrc[&dom_id];
            let oper_time = OpSclDomTime::new(vars, dom_id, temp_id, vlcp_id)?;
            let oper_diff = OpSclDomDiffusion::new(vars, dom_id, temp_id, cond_id, temp_id)?;
            let oper_prd = OpSclDomProduction::new(vars, dom_id, temp_id, hprd_id)?;
            let oper_src = OpSclDomSource::new(vars, dom_id, temp_id, hsrc_id)?;
            self.oper_itr.push((oper_time, oper_diff, oper_prd, oper_src));
        }

        // boundary
        for &bnd_id in &self.temp_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let bnd_temp_id = self.temp_temp[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_temp_id, bnd_temp_id)?;
            self.oper_bnd_temp.push(oper_dir);
        }
        for &bnd_id in &self.hflx_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_temp_id = self.itr_temp[&dom_id];
            let dom_vlcp_id = self.itr_vlcp[&dom_id];
            let dom_cond_id = self.itr_cond[&dom_id];
            let dom_hprd_id = self.itr_hprd[&dom_id];
            let dom_hsrc_id = self.itr_hsrc[&dom_id];
            let bnd_hflx_id = self.hflx_hflx[&bnd_id];
            let oper_time = OpSclBndTime::new(vars, bnd_id, dom_temp_id, dom_vlcp_id)?;
            let oper_diff = OpSclBndDiffusion::new(vars, bnd_id, dom_temp_id, dom_cond_id, dom_temp_id)?;
            let oper_prd = OpSclBndProduction::new(vars, bnd_id, dom_temp_id, dom_hprd_id)?;
            let oper_src = OpSclBndSource::new(vars, bnd_id, dom_temp_id, dom_hsrc_id)?;
            let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_temp_id, bnd_hflx_id)?;
            self.oper_bnd_hflx.push((oper_time, oper_diff, oper_prd, oper_src, oper_neu));
        }

        // interface
        for &itf_id in &self.hcnt_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let dom_temp_a_id = self.itr_temp[&dom_a_id];
            let dom_temp_b_id = self.itr_temp[&dom_b_id];
            let dom_vlcp_a_id = self.itr_vlcp[&dom_a_id];
            let dom_vlcp_b_id = self.itr_vlcp[&dom_b_id];
            let dom_cond_a_id = self.itr_cond[&dom_a_id];
            let dom_cond_b_id = self.itr_cond[&dom_b_id];
            let dom_hprd_a_id = self.itr_hprd[&dom_a_id];
            let dom_hprd_b_id = self.itr_hprd[&dom_b_id];
            let dom_hsrc_a_id = self.itr_hsrc[&dom_a_id];
            let dom_hsrc_b_id = self.itr_hsrc[&dom_b_id];
            let oper_time = OpSclItfTime::new(vars, itf_id, dom_temp_a_id, dom_vlcp_a_id, dom_temp_b_id, dom_vlcp_b_id)?;
            let oper_cont = OpSclItfContinuity::new(vars, itf_id, dom_temp_a_id, dom_temp_b_id)?;
            let oper_diff = OpSclItfDiffusion::new(vars, itf_id, dom_temp_a_id, dom_cond_a_id, dom_temp_a_id, dom_temp_b_id, dom_cond_b_id, dom_temp_b_id)?;
            let oper_prd = OpSclItfProduction::new(vars, itf_id, dom_temp_a_id, dom_hprd_a_id, dom_temp_b_id, dom_hprd_b_id)?;
            let oper_src = OpSclItfSource::new(vars, itf_id, dom_temp_a_id, dom_hsrc_a_id, dom_temp_b_id, dom_hsrc_b_id)?;
            self.oper_itf_hcnt.push((oper_time, oper_cont, oper_diff, oper_prd, oper_src));
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // internal
        for (oper_time, oper_diff, oper_prd, oper_src) in &self.oper_itr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        // boundary
        for oper_dir in &self.oper_bnd_temp {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_time, oper_diff, oper_prd, oper_src, oper_neu) in &self.oper_bnd_hflx {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }

        // interface
        for (oper_time, oper_cont, oper_diff, oper_prd, oper_src) in &self.oper_itf_hcnt {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_cont.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_temp {
            oper_dir.apply_initial(vars);
        }
    }

    fn update_flow(&self, _vars: &mut Variables) {
        // no flow data
    }
}
