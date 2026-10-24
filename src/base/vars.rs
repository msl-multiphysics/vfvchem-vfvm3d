use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Variables {
    // geometry
    pub pts: Points,
    pub nds: Nodes,
    pub tet: Vec<Tetrahedron>,

    // grouping
    pub dom: Vec<Domain>,
    pub bnd: Vec<Boundary>,
    pub itf: Vec<Interface>,

    // scalar
    pub scldom: Vec<ScalarDomain>,
    pub sclbnd: Vec<ScalarBoundary>,
    pub sclitf: Vec<ScalarInterface>,
    pub sclflw: Vec<ScalarFlow>,

    // vector
    pub vecdom: Vec<VectorDomain>,
    pub vecbnd: Vec<VectorBoundary>,
    pub vecitf: Vec<VectorInterface>,

    // convergence scaling
    // reference scales of unknown fields used in the Picard residual and the adaptive error norm (see scaling.rs)
    pub scale_scldom: HashMap<usize, f64>,  // scldom id -> reference scale (for pressure: reference pressure difference)
    pub scale_vecdom: HashMap<usize, f64>,  // vecdom id -> reference scale (velocity magnitude)
    pub scale_off: bool,                    // if true, the unscaled residual |x_new - x_old| / |x_old| is used
}

impl Variables {
    // geometry
    // boundary nodes are spaced by h_min
    // point spacing grows linearly from h_min at the boundary to h_max at a distance ref_len
    pub fn new_obj(file_path: &str, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Result<Variables, VFVChemError> {
        let mut vars = Variables::default();
        vars.nds = Nodes::new_from_obj(file_path, h_min)?;
        vars.pts = Points::new(&vars.nds, h_min, h_max, ref_len, max_try, seed)?;
        Ok(vars)
    }

    pub fn new_bounds(x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Result<Variables, VFVChemError> {
        let mut vars = Variables::default();
        vars.nds = Nodes::new_from_bounds(x_min, y_min, z_min, x_max, y_max, z_max, h_min)?;
        vars.pts = Points::new(&vars.nds, h_min, h_max, ref_len, max_try, seed)?;
        Ok(vars)
    }

    // groupings
    pub fn add_dom(&mut self, grp3d_id: &str) -> Result<usize, VFVChemError> {
        let dom_id = self.dom.len();
        let dom = Domain::new(dom_id, &self.pts, &self.nds, grp3d_id)?;
        let tet = Tetrahedron::new(dom_id, &dom)?;  // tet_id = dom_id;
        self.dom.push(dom);
        self.tet.push(tet);
        Ok(dom_id)
    }

    pub fn add_bnd(&mut self, dom_id: usize, grp2d_id: &str) -> Result<usize, VFVChemError> {
        let bnd_id = self.bnd.len();
        let dom = &self.dom[dom_id];
        let bnd = Boundary::new(bnd_id, &self.nds, dom, grp2d_id)?;
        self.bnd.push(bnd);
        Ok(bnd_id)
    }

    pub fn add_itf(&mut self, dom_a_id: usize, dom_b_id: usize, grp2d_id: &str) -> Result<usize, VFVChemError> {
        let itf_id = self.itf.len();
        let dom_a = &self.dom[dom_a_id];
        let dom_b = &self.dom[dom_b_id];
        let itf = Interface::new(itf_id, &self.nds, dom_a, dom_b, grp2d_id)?;
        self.itf.push(itf);
        Ok(itf_id)
    }

    // scalar
    pub fn add_scldom_cnst(&mut self, dom_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_cnst(scldom_id, dom, value, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_func(&mut self, dom_id: usize, func: SclDomFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_func(scldom_id, dom, &self.scldom, func, unkn_id, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_unkn(&mut self, dom_id: usize, value_init: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_unkn(scldom_id, dom, value_init, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_sclbnd_cnst(&mut self, bnd_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_cnst(sclbnd_id, bnd, value, file_path)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclbnd_func(&mut self, bnd_id: usize, func: SclBndFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_func(sclbnd_id, bnd, &self.scldom, func, unkn_id, file_path)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclitf_cnst(&mut self, itf_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_cnst(sclitf_id, itf, value, file_path)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_sclitf_func(&mut self, itf_id: usize, func: SclItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_func(sclitf_id, itf, &self.scldom, func, unkn_a_id, unkn_b_id, file_path)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_sclflw(&mut self, dom_id: usize, vel_id: usize, pres_id: usize, den_id: usize, visc_id: usize, fce_id: usize) -> Result<usize, VFVChemError> {
        let sclflw_id = self.sclflw.len();
        let dom = &self.dom[dom_id];
        let sclflw = ScalarFlow::new(sclflw_id, dom, vel_id, pres_id, den_id, visc_id, fce_id)?;
        self.sclflw.push(sclflw);
        Ok(sclflw_id)
    }

    // vector
    pub fn add_vecdom_cnst(&mut self, dom_id: usize, value: (f64, f64, f64), file_path: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_cnst(vecdom_id, dom, value, file_path)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_func(&mut self, dom_id: usize, func: VecDomFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_func(vecdom_id, dom, &self.scldom, func, unkn_id, file_path)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecdom_unkn(&mut self, dom_id: usize, value_init: (f64, f64, f64), file_path: &str) -> Result<usize, VFVChemError> {
        let vecdom_id = self.vecdom.len();
        let dom = &self.dom[dom_id];
        let vecdom = VectorDomain::new_unkn(vecdom_id, dom, value_init, file_path)?;
        self.vecdom.push(vecdom);
        Ok(vecdom_id)
    }

    pub fn add_vecbnd_cnst(&mut self, bnd_id: usize, value: (f64, f64, f64), file_path: &str) -> Result<usize, VFVChemError> {
        let vecbnd_id = self.vecbnd.len();
        let bnd = &self.bnd[bnd_id];
        let vecbnd = VectorBoundary::new_cnst(vecbnd_id, bnd, value, file_path)?;
        self.vecbnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn add_vecbnd_func(&mut self, bnd_id: usize, func: VecBndFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let vecbnd_id = self.vecbnd.len();
        let bnd = &self.bnd[bnd_id];
        let vecbnd = VectorBoundary::new_func(vecbnd_id, bnd, &self.scldom, func, unkn_id, file_path)?;
        self.vecbnd.push(vecbnd);
        Ok(vecbnd_id)
    }

    pub fn add_vecitf_cnst(&mut self, itf_id: usize, value: (f64, f64, f64), file_path: &str) -> Result<usize, VFVChemError> {
        let vecitf_id = self.vecitf.len();
        let itf = &self.itf[itf_id];
        let vecitf = VectorInterface::new_cnst(vecitf_id, itf, value, file_path)?;
        self.vecitf.push(vecitf);
        Ok(vecitf_id)
    }

    pub fn add_vecitf_func(&mut self, itf_id: usize, func: VecItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let vecitf_id = self.vecitf.len();
        let itf = &self.itf[itf_id];
        let vecitf = VectorInterface::new_func(vecitf_id, itf, &self.scldom, func, unkn_a_id, unkn_b_id, file_path)?;
        self.vecitf.push(vecitf);
        Ok(vecitf_id)
    }


    // convergence scaling
    pub fn set_scldom_scale(&mut self, scldom_id: usize, scale: f64) -> Result<(), VFVChemError> {
        if scale.is_nan() || scale <= 0.0 { return Err(VFVChemError::VariablesNonPositiveScale { scale }); }
        self.scale_scldom.insert(scldom_id, scale);
        Ok(())
    }

    pub fn set_vecdom_scale(&mut self, vecdom_id: usize, scale: f64) -> Result<(), VFVChemError> {
        if scale.is_nan() || scale <= 0.0 { return Err(VFVChemError::VariablesNonPositiveScale { scale }); }
        self.scale_vecdom.insert(vecdom_id, scale);
        Ok(())
    }

}
