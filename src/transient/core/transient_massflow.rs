use crate::base::prelude::*;
use crate::operator::core::flow_helper;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientMassFlow {
    // mass internal data
    // indexed by component, then by domain
    pub num_comp: usize,                                       // number of components
    pub itr_conc: Vec<HashMap<usize, usize>>,                  // concentration
    pub itr_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // diffusion coefficient; component j -> D_ij
    pub itr_mprd: Vec<HashMap<usize, usize>>,                  // mass production coefficient
    pub itr_msrc: Vec<HashMap<usize, usize>>,                  // mass source

    // mass boundary data
    pub conc_bnd: Vec<Vec<usize>>,              // boundary with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>,  // concentration
    pub mflx_bnd: Vec<Vec<usize>>,              // boundary with diffusive mass flux BC
    pub mflx_mflx: Vec<HashMap<usize, usize>>,  // diffusive mass flux
    pub mout_bnd: Vec<Vec<usize>>,              // boundary with mass outflow BC

    // mass interface data
    pub mcnt_itf: Vec<Vec<usize>>,              // interface with continuity IC

    // flow internal data
    pub itr_fdom: Vec<usize>,             // domain
    pub itr_vel: HashMap<usize, usize>,   // velocity
    pub itr_pres: HashMap<usize, usize>,  // pressure
    pub itr_den: HashMap<usize, usize>,   // density
    pub itr_visc: HashMap<usize, usize>,  // viscosity
    pub itr_fce: HashMap<usize, usize>,   // body force
    pub itr_flw: HashMap<usize, usize>,   // flow data; created in initial_start

    // flow boundary data
    pub vel_bnd: Vec<usize>,              // boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,   // velocity
    pub pres_bnd: Vec<usize>,             // boundary with pressure BC
    pub pres_pres: HashMap<usize, usize>, // pressure
    pub trac_bnd: Vec<usize>,             // boundary with traction BC
    pub trac_pres: HashMap<usize, usize>, // pressure

    // flow interface data
    pub fcnt_itf: Vec<usize>,             // interface with continuity IC

    // reference pressure
    pub pref_dom: Vec<usize>,             // domain with reference pressure
    pub pref_nid: HashMap<usize, usize>,  // node with reference pressure
    pub pref_pres: HashMap<usize, f64>,   // reference pressure

    // mass operators
    pub oper_itr_mass: Vec<(Vec<OpSclDomDiffusion>, OpSclDomProduction, OpSclDomSource)>,
    pub oper_itr_mtim: Vec<OpSclDomTimeUnity>,
    pub oper_itr_mtpr: Vec<OpSclDomTimePrevious>,
    pub oper_itr_madv: Vec<OpSclDomAdvectionUnity>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<(Vec<OpSclBndDiffusion>, OpSclBndProduction, OpSclBndSource, OpSclBndNeumann)>,
    pub oper_bnd_mtim: Vec<OpSclBndTimeUnity>,
    pub oper_bnd_mtpr: Vec<OpSclBndTimePrevious>,
    pub oper_bnd_madv: Vec<OpSclBndAdvectionUnity>,
    pub oper_bnd_mout: Vec<(OpSclBndTimePrevious, Vec<OpSclBndDiffusion>, OpSclBndAdvectionUnity, OpSclBndProduction, OpSclBndSource)>,
    pub oper_itf_mcnt: Vec<(OpSclItfContinuity, Vec<OpSclItfDiffusion>, OpSclItfProduction, OpSclItfSource)>,
    pub oper_itf_mtim: Vec<OpSclItfTimeUnity>,
    pub oper_itf_mtpr: Vec<OpSclItfTimePrevious>,
    pub oper_itf_madv: Vec<OpSclItfAdvectionUnity>,

    // flow operators
    pub oper_itr_flow: Vec<(OpVecDomTime, OpVecDomAdvection, OpVecDomGradient, OpVecDomDiffusion, OpVecDomTranspose, OpVecDomDilatatory, OpVecDomSource, OpSclDomDivergence, OpSclDomPressureStabilizer)>,
    pub oper_bnd_vel: Vec<(OpVecBndDirichlet, OpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_bnd_pres: Vec<(OpSclBndDirichlet, OpVecBndTime, OpVecBndAdvection, OpVecBndGradient, OpVecBndDiffusion, OpVecBndTranspose, OpVecBndDilatatory, OpVecBndSource)>,
    pub oper_bnd_trac: Vec<(OpVecBndTraction, OpVecBndTime, OpVecBndAdvection, OpVecBndDiffusion, OpVecBndTranspose, OpVecBndDilatatory, OpVecBndSource, OpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_itf_fcnt: Vec<(OpVecItfContinuity, OpSclItfContinuity, OpVecItfTime, OpVecItfAdvection, OpVecItfGradient, OpVecItfDiffusion, OpVecItfTranspose, OpVecItfDilatatory, OpVecItfSource, OpSclItfDivergence, OpSclItfPressureStabilizer)>,
    pub oper_pref: Vec<OpSclDomPressureReference>,
}

impl TransientMassFlow {
    pub fn new(num_comp: usize) -> Result<TransientMassFlow, VFVChemError> {
        // mass transfer is solved in conservative form minus c_i times continuity
        // rho * dc_i/dt + rho * vel . grad(c_i) = -div( -D_ij * grad(c_j) ) + R_i
        // on domains with flow, c_i is per unit mass and D_ij includes the density (rho * D_ij)
        // uniform c_i is preserved even if discrete continuity is not satisfied; c_i is conserved once it is
        // on domains without flow, c_i is per unit volume and the time derivative has unit weight

        // input checking
        if num_comp == 0 { return Err(VFVChemError::MassNonPositiveNumComponents { num_comp }) }

        // create struct with data for each component
        let mut phys = TransientMassFlow::default();
        phys.num_comp = num_comp;
        phys.itr_conc = vec![HashMap::new(); num_comp];
        phys.itr_diff = vec![HashMap::new(); num_comp];
        phys.itr_mprd = vec![HashMap::new(); num_comp];
        phys.itr_msrc = vec![HashMap::new(); num_comp];
        phys.conc_bnd = vec![Vec::new(); num_comp];
        phys.conc_conc = vec![HashMap::new(); num_comp];
        phys.mflx_bnd = vec![Vec::new(); num_comp];
        phys.mflx_mflx = vec![HashMap::new(); num_comp];
        phys.mout_bnd = vec![Vec::new(); num_comp];
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

    pub fn add_mout_bnd(&mut self, comp_id: usize, bnd_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mout_bnd[comp_id].push(bnd_id);
        Ok(())
    }

    pub fn add_mcnt_itf(&mut self, comp_id: usize, itf_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mcnt_itf[comp_id].push(itf_id);
        Ok(())
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<(), VFVChemError> {
        self.itr_fdom.push(dom_id);
        self.itr_vel.insert(dom_id, vel_id);
        self.itr_pres.insert(dom_id, pres_id);
        self.itr_den.insert(dom_id, den_id);
        self.itr_visc.insert(dom_id, visc_id);
        self.itr_fce.insert(dom_id, fce_id);
        Ok(())
    }

    pub fn add_vel_bnd(&mut self, bnd_id: usize, vel_id: usize) -> Result<(), VFVChemError> {
        self.vel_bnd.push(bnd_id);
        self.vel_vel.insert(bnd_id, vel_id);
        Ok(())
    }

    pub fn add_pres_bnd(&mut self, bnd_id: usize, pres_id: usize) -> Result<(), VFVChemError> {
        self.pres_bnd.push(bnd_id);
        self.pres_pres.insert(bnd_id, pres_id);
        Ok(())
    }

    pub fn add_trac_bnd(&mut self, bnd_id: usize, pres_id: usize) -> Result<(), VFVChemError> {
        self.trac_bnd.push(bnd_id);
        self.trac_pres.insert(bnd_id, pres_id);
        Ok(())
    }

    pub fn add_fcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.fcnt_itf.push(itf_id);
        Ok(())
    }

    pub fn set_pres_ref(&mut self, dom_id: usize, nid: usize, pres: f64) -> Result<(), VFVChemError> {
        self.pref_dom.push(dom_id);
        self.pref_nid.insert(dom_id, nid);
        self.pref_pres.insert(dom_id, pres);
        Ok(())
    }

}

impl TransientBase for TransientMassFlow {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown vectors and scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in &self.itr_vel {
            vars.vecdom[vel_id].vecdom_type = VecDomType::Unknown { start: xid };
            xid += 3 * vars.dom[dom_id].num_tot();
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            vars.scldom[pres_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }
        for itr_conc in &self.itr_conc {
        for (&dom_id, &conc_id) in itr_conc {
            vars.scldom[conc_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }}

        // create flow data
        // replace in place on a re-solve to avoid orphaned flow data
        for &dom_id in &self.itr_fdom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];
            match self.itr_flw.get(&dom_id) {
                Some(&flw_id) => { vars.sclflw[flw_id] = ScalarFlow::new(flw_id, &vars.dom[dom_id], vel_id, pres_id, den_id, visc_id, fce_id)?; }
                None => { self.itr_flw.insert(dom_id, vars.add_sclflw(dom_id, vel_id, pres_id, den_id, visc_id, fce_id)?); }
            }
        }

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr_mass.clear();
        self.oper_itr_mtim.clear();
        self.oper_itr_mtpr.clear();
        self.oper_itr_madv.clear();
        self.oper_bnd_conc.clear();
        self.oper_bnd_mflx.clear();
        self.oper_bnd_mtim.clear();
        self.oper_bnd_mtpr.clear();
        self.oper_bnd_madv.clear();
        self.oper_bnd_mout.clear();
        self.oper_itf_mcnt.clear();
        self.oper_itf_mtim.clear();
        self.oper_itf_mtpr.clear();
        self.oper_itf_madv.clear();
        self.oper_itr_flow.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_pres.clear();
        self.oper_bnd_trac.clear();
        self.oper_itf_fcnt.clear();
        self.oper_pref.clear();

        for comp_i in 0..self.num_comp {
            // mass internal
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

                // advection only on domains with flow; the Peclet number uses the self-diffusivity
                // the time derivative is weighted by the previous density on domains with flow
                match self.itr_flw.get(&dom_id) {
                    Some(&flw_id) => {
                        let diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let den_id = self.itr_den[&dom_id];
                        let oper_time = OpSclDomTimePrevious::new(vars, dom_id, conc_id, den_id)?;
                        let oper_adv = OpSclDomAdvectionUnity::new(vars, dom_id, conc_id, diff_id, flw_id)?;
                        self.oper_itr_mtpr.push(oper_time);
                        self.oper_itr_madv.push(oper_adv);
                    }
                    None => {
                        let oper_time = OpSclDomTimeUnity::new(vars, dom_id, conc_id)?;
                        self.oper_itr_mtim.push(oper_time);
                    }
                }
            }

            // mass boundary
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

                // advection only on domains with flow; the Neumann value is the diffusive flux
                // the time derivative is weighted by the previous density on domains with flow
                match self.itr_flw.get(&dom_id) {
                    Some(&dom_flw_id) => {
                        let dom_diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let dom_vel_id = self.itr_vel[&dom_id];
                        let dom_den_id = self.itr_den[&dom_id];
                        let oper_time = OpSclBndTimePrevious::new(vars, bnd_id, dom_conc_id, dom_den_id)?;
                        let oper_adv = OpSclBndAdvectionUnity::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_vel_id, dom_flw_id)?;
                        self.oper_bnd_mtpr.push(oper_time);
                        self.oper_bnd_madv.push(oper_adv);
                    }
                    None => {
                        let oper_time = OpSclBndTimeUnity::new(vars, bnd_id, dom_conc_id)?;
                        self.oper_bnd_mtim.push(oper_time);
                    }
                }
            }
            for &bnd_id in &self.mout_bnd[comp_i] {
                // zero diffusive flux; mass leaves by advection only
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let dom_mprd_id = self.itr_mprd[comp_i][&dom_id];
                let dom_msrc_id = self.itr_msrc[comp_i][&dom_id];
                let dom_vel_id = self.itr_vel[&dom_id];
                let dom_den_id = self.itr_den[&dom_id];
                let dom_flw_id = self.itr_flw[&dom_id];
                let oper_time = OpSclBndTimePrevious::new(vars, bnd_id, dom_conc_id, dom_den_id)?;
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let dom_conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclBndDiffusion::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_conc_j_id)?);
                }
                let dom_diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                let oper_adv = OpSclBndAdvectionUnity::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_vel_id, dom_flw_id)?;
                let oper_prd = OpSclBndProduction::new(vars, bnd_id, dom_conc_id, dom_mprd_id)?;
                let oper_src = OpSclBndSource::new(vars, bnd_id, dom_conc_id, dom_msrc_id)?;
                self.oper_bnd_mout.push((oper_time, oper_diff, oper_adv, oper_prd, oper_src));
            }

            // mass interface
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

                // the time derivative is weighted by the previous density if both sides have flow
                // c is per unit mass with flow and per unit volume without, so both sides must match
                match (self.itr_den.get(&dom_a_id), self.itr_den.get(&dom_b_id)) {
                    (Some(&den_a_id), Some(&den_b_id)) => {
                        let oper_time = OpSclItfTimePrevious::new(vars, itf_id, dom_conc_a_id, den_a_id, dom_conc_b_id, den_b_id)?;
                        self.oper_itf_mtpr.push(oper_time);
                    }
                    (None, None) => {
                        let oper_time = OpSclItfTimeUnity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id)?;
                        self.oper_itf_mtim.push(oper_time);
                    }
                    _ => return Err(VFVChemError::MassInterfaceMixedFlow { itf_id }),
                }

                // advection on each side with flow
                for dom_id in [dom_a_id, dom_b_id] {
                    if let Some(&flw_id) = self.itr_flw.get(&dom_id) {
                        let diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let vel_id = self.itr_vel[&dom_id];
                        let oper_adv = OpSclItfAdvectionUnity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id, diff_id, vel_id, flw_id)?;
                        self.oper_itf_madv.push(oper_adv);
                    }
                }
            }
        }

        // flow internal
        for &dom_id in &self.itr_fdom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];
            let flw_id = self.itr_flw[&dom_id];
            let oper_time = OpVecDomTime::new(vars, dom_id, vel_id, den_id)?;
            let oper_adv = OpVecDomAdvection::new(vars, dom_id, vel_id, den_id, visc_id, flw_id)?;
            let oper_grad = OpVecDomGradient::new(vars, dom_id, vel_id, pres_id)?;
            let oper_diff = OpVecDomDiffusion::new(vars, dom_id, vel_id, visc_id)?;
            let oper_dtrn = OpVecDomTranspose::new(vars, dom_id, vel_id, visc_id)?;
            let oper_ddil = OpVecDomDilatatory::new(vars, dom_id, vel_id, visc_id)?;
            let oper_src = OpVecDomSource::new(vars, dom_id, vel_id, fce_id)?;
            let oper_div = OpSclDomDivergence::new(vars, dom_id, vel_id, pres_id, den_id)?;
            let oper_pres = OpSclDomPressureStabilizer::new(vars, dom_id, pres_id, den_id, fce_id, flw_id)?;
            self.oper_itr_flow.push((oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
        }

        // flow boundary
        for &bnd_id in &self.vel_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = OpVecBndDirichlet::new(vars, bnd_id, dom_vel_id, bnd_vel_id)?;
            let oper_div = OpSclBndDivergence::new(vars, bnd_id, dom_vel_id, dom_pres_id, dom_den_id)?;
            let oper_pres = OpSclBndPressureStabilizer::new(vars, bnd_id, dom_pres_id, dom_den_id, dom_fce_id, dom_flw_id)?;
            self.oper_bnd_vel.push((oper_dir, oper_div, oper_pres));
        }
        for &bnd_id in &self.pres_bnd {
            // momentum over the whole node control volumes; mu * grad(vel) . n is do-nothing on the boundary half-edges
            // continuity is replaced with the boundary pressure
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_visc_id = self.itr_visc[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_pres_id = self.pres_pres[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_pres_id, bnd_pres_id)?;
            let oper_time = OpVecBndTime::new(vars, bnd_id, dom_vel_id, dom_den_id)?;
            let oper_adv = OpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_flw_id)?;
            let oper_grad = OpVecBndGradient::new(vars, bnd_id, dom_vel_id, dom_pres_id)?;
            let oper_diff = OpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_dtrn = OpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_ddil = OpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            self.oper_bnd_pres.push((oper_dir, oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src));
        }
        for &bnd_id in &self.trac_bnd {
            // momentum over the whole node control volumes; the boundary pressure sets the traction on the boundary half-edges
            // continuity is kept, so pressure stays unknown on the boundary
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_pres_id = self.itr_pres[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let dom_visc_id = self.itr_visc[&dom_id];
            let dom_fce_id = self.itr_fce[&dom_id];
            let dom_flw_id = self.itr_flw[&dom_id];
            let bnd_pres_id = self.trac_pres[&bnd_id];
            let oper_trac = OpVecBndTraction::new(vars, bnd_id, dom_vel_id, dom_pres_id, bnd_pres_id)?;
            let oper_time = OpVecBndTime::new(vars, bnd_id, dom_vel_id, dom_den_id)?;
            let oper_adv = OpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_flw_id)?;
            let oper_diff = OpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_dtrn = OpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_ddil = OpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            let oper_div = OpSclBndDivergence::new(vars, bnd_id, dom_vel_id, dom_pres_id, dom_den_id)?;
            let oper_pres = OpSclBndPressureStabilizer::new(vars, bnd_id, dom_pres_id, dom_den_id, dom_fce_id, dom_flw_id)?;
            self.oper_bnd_trac.push((oper_trac, oper_time, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
        }

        // flow interface
        // rows of domain a set continuous velocity and pressure; rows of domain b balance momentum and mass over both halves
        for &itf_id in &self.fcnt_itf {
            let dom_a_id = vars.itf[itf_id].dom_a_id;
            let dom_b_id = vars.itf[itf_id].dom_b_id;
            let (vel_a_id, vel_b_id) = (self.itr_vel[&dom_a_id], self.itr_vel[&dom_b_id]);
            let (pres_a_id, pres_b_id) = (self.itr_pres[&dom_a_id], self.itr_pres[&dom_b_id]);
            let (den_a_id, den_b_id) = (self.itr_den[&dom_a_id], self.itr_den[&dom_b_id]);
            let (visc_a_id, visc_b_id) = (self.itr_visc[&dom_a_id], self.itr_visc[&dom_b_id]);
            let (fce_a_id, fce_b_id) = (self.itr_fce[&dom_a_id], self.itr_fce[&dom_b_id]);
            let (flw_a_id, flw_b_id) = (self.itr_flw[&dom_a_id], self.itr_flw[&dom_b_id]);
            let oper_vcnt = OpVecItfContinuity::new(vars, itf_id, vel_a_id, vel_b_id)?;
            let oper_pcnt = OpSclItfContinuity::new(vars, itf_id, pres_a_id, pres_b_id)?;
            let oper_time = OpVecItfTime::new(vars, itf_id, vel_a_id, den_a_id, vel_b_id, den_b_id)?;
            let oper_adv = OpVecItfAdvection::new(vars, itf_id, vel_a_id, den_a_id, visc_a_id, flw_a_id, vel_b_id, den_b_id, visc_b_id, flw_b_id)?;
            let oper_grad = OpVecItfGradient::new(vars, itf_id, vel_a_id, pres_a_id, vel_b_id, pres_b_id)?;
            let oper_diff = OpVecItfDiffusion::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_dtrn = OpVecItfTranspose::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_ddil = OpVecItfDilatatory::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_src = OpVecItfSource::new(vars, itf_id, vel_a_id, fce_a_id, vel_b_id, fce_b_id)?;
            let oper_div = OpSclItfDivergence::new(vars, itf_id, vel_a_id, pres_a_id, den_a_id, vel_b_id, pres_b_id, den_b_id)?;
            let oper_pres = OpSclItfPressureStabilizer::new(vars, itf_id, pres_a_id, den_a_id, fce_a_id, flw_a_id, pres_b_id, den_b_id, fce_b_id, flw_b_id)?;
            self.oper_itf_fcnt.push((oper_vcnt, oper_pcnt, oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
        }

        // reference pressure
        for &dom_id in &self.pref_dom {
            let pres_id = self.itr_pres[&dom_id];
            let nid = self.pref_nid[&dom_id];
            let pres = self.pref_pres[&dom_id];
            let oper_pref = OpSclDomPressureReference::new(vars, dom_id, pres_id, nid, pres)?;
            self.oper_pref.push(oper_pref);
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // flow internal
        for (oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_itr_flow {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // flow boundary
        for (oper_dir, oper_div, oper_pres) in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_dir, oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src) in &self.oper_bnd_pres {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_trac, oper_time, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_bnd_trac {
            oper_trac.apply(vars, a_triplet, b_vec)?;
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // flow interface
        for (oper_vcnt, oper_pcnt, oper_time, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_itf_fcnt {
            oper_vcnt.apply(vars, a_triplet, b_vec)?;
            oper_pcnt.apply(vars, a_triplet, b_vec)?;
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // mass internal
        for (oper_diff, oper_prd, oper_src) in &self.oper_itr_mass {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_itr_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_itr_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_itr_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }

        // mass boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_prd, oper_src, oper_neu) in &self.oper_bnd_mflx {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_bnd_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_bnd_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_bnd_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_time, oper_diff, oper_adv, oper_prd, oper_src) in &self.oper_bnd_mout {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        // mass interface
        for (oper_cont, oper_diff, oper_prd, oper_src) in &self.oper_itf_mcnt {
            oper_cont.apply(vars, a_triplet, b_vec)?;
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_itf_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_itf_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_itf_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }

        // reference pressure; replaces equations, so it is applied last
        for oper_pref in &self.oper_pref {
            oper_pref.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply_initial(vars);
        }
        for (oper_dir, ..) in &self.oper_bnd_vel {
            oper_dir.apply_initial(vars);
        }
        for (oper_dir, ..) in &self.oper_bnd_pres {
            oper_dir.apply_initial(vars);
        }
    }

    fn update_flow(&self, vars: &mut Variables) {
        // Rhie-Chow face mass flux of a free fluid
        for &dom_id in &self.itr_fdom {
            flow_helper::update_flow(vars, self.itr_flw[&dom_id]);
        }
    }
}
