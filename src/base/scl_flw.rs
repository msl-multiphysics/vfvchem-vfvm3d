use crate::base::prelude::*;

#[derive(Default)]
pub struct ScalarFlow {
    // ids
    pub sclflw_id: usize,
    pub dom_id: usize,  // domain this is attached to

    // vectors
    pub vel_id: usize,  // velocity (vel)
    pub fce_id: usize,  // body force (fce)

    // scalars
    pub pres_id: usize,  // pressure (p)
    pub den_id: usize,   // density (rho)
    pub visc_id: usize,  // viscosity (mu)

    // flow data
    // shared by the momentum and continuity operators so that both use the same face mass flux
    pub rc_coeff: Vec<f64>,   // vid -> V / a_P for the Rhie-Chow correction
    pub flux: Vec<Vec<f64>>,  // vid, loc -> face mass flux
}

impl ScalarFlow {
    pub fn new(sclflw_id: usize, dom: &Domain, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<ScalarFlow, VFVChemError> {
        // default class
        let mut sclflw = ScalarFlow::default();
        sclflw.sclflw_id = sclflw_id;
        sclflw.dom_id = dom.dom_id;

        // field ids
        sclflw.vel_id = vel_id;
        sclflw.fce_id = fce_id;
        sclflw.pres_id = pres_id;
        sclflw.den_id = den_id;
        sclflw.visc_id = visc_id;

        // flow data is computed by the solver after the tetrahedralization; see flow_helper
        sclflw.rc_coeff = Vec::new();
        sclflw.flux = Vec::new();

        Ok(sclflw)
    }

}
