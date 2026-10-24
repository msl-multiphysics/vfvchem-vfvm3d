use crate::base::prelude::*;

// input: time, (x, y, z), slice of unknown values at a point
// output: function value at the point
pub type VecDomFunc = Box<dyn Fn(f64, (f64, f64, f64), &[f64]) -> (f64, f64, f64) + Send + Sync>;

#[derive(Default)]
pub enum VecDomType {
    #[default] Constant,
    Function {
        func: VecDomFunc,  
        unkn_id: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

#[derive(Default)]
pub struct VectorDomain {
    // ids
    pub vecdom_id: usize,
    pub dom_id: usize,  // domain this is attached to

    // point and node data
    pub vecdom_type: VecDomType,
    pub value: Vec<(f64, f64, f64)>,  // vid -> (u, v, w)
    pub value_prev: Vec<(f64, f64, f64)>,

    // gradient data
    pub grad: Vec<(f64, f64, f64, f64, f64, f64, f64, f64, f64)>,  // vid -> (du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx, dw/dy, dw/dz)
    pub grad_prev: Vec<(f64, f64, f64, f64, f64, f64, f64, f64, f64)>,

    // output file
    pub file_name: String,  // file path without extension
    pub file_type: String,  // file extension
}

impl VectorDomain {
    pub fn new_cnst(vecdom_id: usize, dom: &Domain, value: (f64, f64, f64), file_path: &str) -> Result<VectorDomain, VFVChemError> {
        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point and node data
        vecdom.vecdom_type = VecDomType::Constant;
        vecdom.value = vec![value; dom.num_tot()];
        vecdom.value_prev = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); dom.num_tot()];
        vecdom.grad_prev = vecdom.grad.clone();

        // output file
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::VectorUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            vecdom.file_name = file_name.to_string();
            vecdom.file_type = file_type.to_string();
        }

        Ok(vecdom)
    }

    pub fn new_func(vecdom_id: usize, dom: &Domain, scldom_all: &Vec<ScalarDomain>, func: VecDomFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<VectorDomain, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::VectorFunctionDependsOnNonUnknown { geom: "Domain", scldom_id: i }) }
            }
        }

        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point and node data
        vecdom.vecdom_type = VecDomType::Function { func, unkn_id };
        vecdom.value = vec![(f64::NAN, f64::NAN, f64::NAN); dom.num_tot()];  // placeholder
        vecdom.value_prev = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN); dom.num_tot()];  // placeholder
        vecdom.grad_prev = vecdom.grad.clone();

        // output file
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::VectorUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            vecdom.file_name = file_name.to_string();
            vecdom.file_type = file_type.to_string();
        }

        Ok(vecdom)
    }

    pub fn new_unkn(vecdom_id: usize, dom: &Domain, value_init: (f64, f64, f64), file_path: &str) -> Result<VectorDomain, VFVChemError> {
        // default class
        let mut vecdom = VectorDomain::default();
        vecdom.vecdom_id = vecdom_id;
        vecdom.dom_id = dom.dom_id;

        // point and node data
        vecdom.vecdom_type = VecDomType::Unknown { start: usize::MAX };  // placeholder
        vecdom.value = vec![value_init; dom.num_tot()];
        vecdom.value_prev = vecdom.value.clone();

        // gradient data
        vecdom.grad = vec![(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); dom.num_tot()];
        vecdom.grad_prev = vecdom.grad.clone();

        // output file
        if file_path == "" {
            vecdom.file_name = String::new();
            vecdom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::VectorUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            vecdom.file_name = file_name.to_string();
            vecdom.file_type = file_type.to_string();
        }

        Ok(vecdom)
    }

    pub fn write(&self, pts: &Points, nds: &Nodes, dom: &Domain, tet: &Tetrahedron, ts: usize) -> Result<(), VFVChemError> {
        match self.file_type.as_str() {
            "" => Ok(()),  // no output if blank
            "csv" => write_vecdom_csv(&self, pts, nds, dom, ts),
            "vtu" => write_vecdom_vtu(&self, pts, nds, dom, tet, ts),
            _ => Err(VFVChemError::VectorUnsupportedOutputType { geom: "Domain", file_type: self.file_type.clone() })
        }
    }

    pub fn update_prev(&mut self) {
        self.value_prev = self.value.clone();
        self.grad_prev = self.grad.clone();
    }

    pub fn update_grad(&mut self, tet: &Tetrahedron) {
        // initialize
        self.grad.clear();

        // iterate through neighbors
        for (vid_c, list_n) in tet.ngh_vid.iter().enumerate() {
            let mut grad_ux = 0.0;
            let mut grad_uy = 0.0;
            let mut grad_uz = 0.0;
            let mut grad_vx = 0.0;
            let mut grad_vy = 0.0;
            let mut grad_vz = 0.0;
            let mut grad_wx = 0.0;
            let mut grad_wy = 0.0;
            let mut grad_wz = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                let (w_x, w_y, w_z) = tet.grad_wgt[vid_c][loc];
                let delta_u = self.value[vid_n].0 - self.value[vid_c].0;
                let delta_v = self.value[vid_n].1 - self.value[vid_c].1;
                let delta_w = self.value[vid_n].2 - self.value[vid_c].2;
                grad_ux += w_x * delta_u;
                grad_uy += w_y * delta_u;
                grad_uz += w_z * delta_u;
                grad_vx += w_x * delta_v;
                grad_vy += w_y * delta_v;
                grad_vz += w_z * delta_v;
                grad_wx += w_x * delta_w;
                grad_wy += w_y * delta_w;
                grad_wz += w_z * delta_w;
            }
            self.grad.push((grad_ux, grad_uy, grad_uz, grad_vx, grad_vy, grad_vz, grad_wx, grad_wy, grad_wz));
        }
    }

    pub fn update_func(vars: &mut Variables, vecdom_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as vecdom_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type vecdom
        
        // compute updated values
        let vecdom = &vars.vecdom[vecdom_id];
        match &vecdom.vecdom_type {
            VecDomType::Function { func, unkn_id } => {
                // get number of points
                let num_tot = vecdom.value.len();
                let dom = &vars.dom[vecdom.dom_id];

                // collect positions
                let mut pos = Vec::with_capacity(num_tot);  // vid -> (x, y, z)
                for &gpt in &dom.pts_gpt {
                    pos.push(vars.pts.pts_pos[gpt]);
                }
                for &gnd in &dom.nds_gnd {
                    pos.push(vars.nds.nds_pos[gnd]);
                }
                
                // collect values from unknown scalars
                let mut value_u = Vec::with_capacity(num_tot);  // vid -> values of unknowns
                for vid in 0..num_tot {
                    let mut value_sub = Vec::with_capacity(unkn_id.len());  // values of unknowns
                    for &i in unkn_id {
                        let scldom_u = &vars.scldom[i];
                        value_sub.push(scldom_u.value[vid]);
                    }
                    value_u.push(value_sub);
                }

                // compute function value
                let mut value_f = Vec::with_capacity(num_tot);
                for vid in 0..num_tot {
                    let value_sub = func(t, pos[vid], value_u[vid].as_slice());
                    value_f.push(value_sub)
                }

                // modify value array in vecdom
                let vecdom = &mut vars.vecdom[vecdom_id];
                vecdom.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
