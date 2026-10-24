use crate::base::prelude::*;
use crate::operator::core::flow_helper;
use crate::operator::prelude::*;
use crate::steady::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct SteadyFlow {
    // internal data
    pub itr_dom: Vec<usize>,              // domain
    pub itr_vel: HashMap<usize, usize>,   // velocity
    pub itr_pres: HashMap<usize, usize>,  // pressure
    pub itr_den: HashMap<usize, usize>,   // density
    pub itr_visc: HashMap<usize, usize>,  // viscosity
    pub itr_fce: HashMap<usize, usize>,   // body force
    pub itr_flw: HashMap<usize, usize>,   // flow data; created in initial_start

    // boundary data
    pub vel_bnd: Vec<usize>,              // boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,   // velocity
    pub pres_bnd: Vec<usize>,             // boundary with pressure BC
    pub pres_pres: HashMap<usize, usize>, // pressure
    pub trac_bnd: Vec<usize>,             // boundary with traction BC
    pub trac_pres: HashMap<usize, usize>, // pressure

    // interface data
    pub fcnt_itf: Vec<usize>,             // interface with continuity IC

    // reference pressure
    pub pref_dom: Vec<usize>,             // domain with reference pressure
    pub pref_nid: HashMap<usize, usize>,  // node with reference pressure
    pub pref_pres: HashMap<usize, f64>,   // reference pressure

    // operators
    pub oper_itr: Vec<(OpVecDomAdvection, OpVecDomGradient, OpVecDomDiffusion, OpVecDomTranspose, OpVecDomDilatatory, OpVecDomSource, OpSclDomDivergence, OpSclDomPressureStabilizer)>,
    pub oper_bnd_vel: Vec<(OpVecBndDirichlet, OpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_bnd_pres: Vec<(OpSclBndDirichlet, OpVecBndAdvection, OpVecBndGradient, OpVecBndDiffusion, OpVecBndTranspose, OpVecBndDilatatory, OpVecBndSource)>,
    pub oper_bnd_trac: Vec<(OpVecBndTraction, OpVecBndAdvection, OpVecBndDiffusion, OpVecBndTranspose, OpVecBndDilatatory, OpVecBndSource, OpSclBndDivergence, OpSclBndPressureStabilizer)>,
    pub oper_itf_fcnt: Vec<(OpVecItfContinuity, OpSclItfContinuity, OpVecItfAdvection, OpVecItfGradient, OpVecItfDiffusion, OpVecItfTranspose, OpVecItfDilatatory, OpVecItfSource, OpSclItfDivergence, OpSclItfPressureStabilizer)>,
    pub oper_pref: Vec<OpSclDomPressureReference>,
}

impl SteadyFlow {
    pub fn new() -> Result<SteadyFlow, VFVChemError> {
        let phys = SteadyFlow::default();
        Ok(phys)
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<(), VFVChemError> {
        self.itr_dom.push(dom_id);
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

impl SteadyBase for SteadyFlow {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in &self.itr_vel {
            vars.vecdom[vel_id].vecdom_type = VecDomType::Unknown { start: xid };
            xid += 3 * vars.dom[dom_id].num_tot();
        }
        for (&dom_id, &pres_id) in &self.itr_pres {
            vars.scldom[pres_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }

        // create flow data
        // replace in place on a re-solve to avoid orphaned flow data
        for &dom_id in &self.itr_dom {
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
        self.oper_itr.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_pres.clear();
        self.oper_bnd_trac.clear();
        self.oper_itf_fcnt.clear();
        self.oper_pref.clear();

        // internal
        for &dom_id in &self.itr_dom {
            let vel_id = self.itr_vel[&dom_id];
            let pres_id = self.itr_pres[&dom_id];
            let den_id = self.itr_den[&dom_id];
            let visc_id = self.itr_visc[&dom_id];
            let fce_id = self.itr_fce[&dom_id];
            let flw_id = self.itr_flw[&dom_id];
            let oper_adv = OpVecDomAdvection::new(vars, dom_id, vel_id, den_id, visc_id, flw_id)?;
            let oper_grad = OpVecDomGradient::new(vars, dom_id, vel_id, pres_id)?;
            let oper_diff = OpVecDomDiffusion::new(vars, dom_id, vel_id, visc_id)?;
            let oper_dtrn = OpVecDomTranspose::new(vars, dom_id, vel_id, visc_id)?;
            let oper_ddil = OpVecDomDilatatory::new(vars, dom_id, vel_id, visc_id)?;
            let oper_src = OpVecDomSource::new(vars, dom_id, vel_id, fce_id)?;
            let oper_div = OpSclDomDivergence::new(vars, dom_id, vel_id, pres_id, den_id)?;
            let oper_pres = OpSclDomPressureStabilizer::new(vars, dom_id, pres_id, den_id, fce_id, flw_id)?;
            self.oper_itr.push((oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
        }

        // boundary
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
            let oper_adv = OpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_flw_id)?;
            let oper_grad = OpVecBndGradient::new(vars, bnd_id, dom_vel_id, dom_pres_id)?;
            let oper_diff = OpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_dtrn = OpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_ddil = OpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            self.oper_bnd_pres.push((oper_dir, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src));
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
            let oper_adv = OpVecBndAdvection::new(vars, bnd_id, dom_vel_id, dom_den_id, dom_visc_id, dom_flw_id)?;
            let oper_diff = OpVecBndDiffusion::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_dtrn = OpVecBndTranspose::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_ddil = OpVecBndDilatatory::new(vars, bnd_id, dom_vel_id, dom_visc_id)?;
            let oper_src = OpVecBndSource::new(vars, bnd_id, dom_vel_id, dom_fce_id)?;
            let oper_div = OpSclBndDivergence::new(vars, bnd_id, dom_vel_id, dom_pres_id, dom_den_id)?;
            let oper_pres = OpSclBndPressureStabilizer::new(vars, bnd_id, dom_pres_id, dom_den_id, dom_fce_id, dom_flw_id)?;
            self.oper_bnd_trac.push((oper_trac, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
        }

        // interface
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
            let oper_adv = OpVecItfAdvection::new(vars, itf_id, vel_a_id, den_a_id, visc_a_id, flw_a_id, vel_b_id, den_b_id, visc_b_id, flw_b_id)?;
            let oper_grad = OpVecItfGradient::new(vars, itf_id, vel_a_id, pres_a_id, vel_b_id, pres_b_id)?;
            let oper_diff = OpVecItfDiffusion::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_dtrn = OpVecItfTranspose::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_ddil = OpVecItfDilatatory::new(vars, itf_id, vel_a_id, visc_a_id, vel_b_id, visc_b_id)?;
            let oper_src = OpVecItfSource::new(vars, itf_id, vel_a_id, fce_a_id, vel_b_id, fce_b_id)?;
            let oper_div = OpSclItfDivergence::new(vars, itf_id, vel_a_id, pres_a_id, den_a_id, vel_b_id, pres_b_id, den_b_id)?;
            let oper_pres = OpSclItfPressureStabilizer::new(vars, itf_id, pres_a_id, den_a_id, fce_a_id, flw_a_id, pres_b_id, den_b_id, fce_b_id, flw_b_id)?;
            self.oper_itf_fcnt.push((oper_vcnt, oper_pcnt, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres));
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

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // internal
        for (oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_itr {
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // boundary
        for (oper_dir, oper_div, oper_pres) in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_dir, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src) in &self.oper_bnd_pres {
            oper_dir.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_trac, oper_adv, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_bnd_trac {
            oper_trac.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // interface
        for (oper_vcnt, oper_pcnt, oper_adv, oper_grad, oper_diff, oper_dtrn, oper_ddil, oper_src, oper_div, oper_pres) in &self.oper_itf_fcnt {
            oper_vcnt.apply(vars, a_triplet, b_vec)?;
            oper_pcnt.apply(vars, a_triplet, b_vec)?;
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_grad.apply(vars, a_triplet, b_vec)?;
            oper_diff.apply(vars, a_triplet, b_vec)?;
            oper_dtrn.apply(vars, a_triplet, b_vec)?;
            oper_ddil.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_pres.apply(vars, a_triplet, b_vec)?;
        }

        // reference pressure; replaces equations, so it is applied last
        for oper_pref in &self.oper_pref {
            oper_pref.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for (oper_dir, ..) in &self.oper_bnd_vel {
            oper_dir.apply_initial(vars);
        }
        for (oper_dir, ..) in &self.oper_bnd_pres {
            oper_dir.apply_initial(vars);
        }
    }

    fn update_flow(&self, vars: &mut Variables) {
        // Rhie-Chow face mass flux of a free fluid
        for &dom_id in &self.itr_dom {
            flow_helper::update_flow(vars, self.itr_flw[&dom_id]);
        }
    }
}
