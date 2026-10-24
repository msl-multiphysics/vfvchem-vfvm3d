use crate::base::prelude::*;

// input: time, (x, y, z), slice of unknown values at a point
// output: function value at the point
pub type SclDomFunc = Box<dyn Fn(f64, (f64, f64, f64), &[f64]) -> f64 + Send + Sync>;

#[derive(Default)]
pub enum SclDomType {
    #[default] Constant,
    Function {
        func: SclDomFunc,  
        unkn_id: Vec<usize>,
    },
    Unknown {
        start: usize,
    },
}

#[derive(Default)]
pub struct ScalarDomain {
    // ids
    pub scldom_id: usize,
    pub dom_id: usize,  // domain this is attached to

    // value data
    pub scldom_type: SclDomType,
    pub value: Vec<f64>,  // vid -> value (v)
    pub value_prev: Vec<f64>,

    // gradient data
    pub grad: Vec<(f64, f64, f64)>,  // vid -> (dv/dx, dv/dy, dv/dz)
    pub grad_prev: Vec<(f64, f64, f64)>,

    // output file
    pub file_name: String,  // file path without extension
    pub file_type: String,  // file extension
}

impl ScalarDomain {
    pub fn new_cnst(scldom_id: usize, dom: &Domain, value: f64, file_path: &str) -> Result<ScalarDomain, VFVChemError> {
        // default class
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point and node data
        scldom.scldom_type = SclDomType::Constant;
        scldom.value = vec![value; dom.num_tot()];
        scldom.value_prev = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(0.0, 0.0, 0.0); dom.num_tot()];
        scldom.grad_prev = scldom.grad.clone();

        // output file
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            scldom.file_name = file_name.to_string();
            scldom.file_type = file_type.to_string();
        }

        Ok(scldom)
    }

    pub fn new_func(scldom_id: usize, dom: &Domain, scldom_all: &Vec<ScalarDomain>, func: SclDomFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<ScalarDomain, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::ScalarFunctionDependsOnNonUnknown { geom: "Domain", scldom_id: i }) }
            }
        }

        // default class
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point data
        scldom.scldom_type = SclDomType::Function { func, unkn_id };
        scldom.value = vec![f64::NAN; dom.num_tot()];  // placeholder
        scldom.value_prev = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(f64::NAN, f64::NAN, f64::NAN); dom.num_tot()];  // placeholder
        scldom.grad_prev = scldom.grad.clone();

        // output file
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            scldom.file_name = file_name.to_string();
            scldom.file_type = file_type.to_string();
        }

        Ok(scldom)
    }

    pub fn new_unkn(scldom_id: usize, dom: &Domain, value_init: f64, file_path: &str) -> Result<ScalarDomain, VFVChemError> {
        // default class
        let mut scldom = ScalarDomain::default();
        scldom.scldom_id = scldom_id;
        scldom.dom_id = dom.dom_id;

        // point data
        scldom.scldom_type = SclDomType::Unknown { start: usize::MAX };  // placeholder
        scldom.value = vec![value_init; dom.num_tot()];
        scldom.value_prev = scldom.value.clone();

        // gradient data
        scldom.grad = vec![(0.0, 0.0, 0.0); dom.num_tot()];
        scldom.grad_prev = scldom.grad.clone();

        // output file
        if file_path == "" {
            scldom.file_name = String::new();
            scldom.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Domain", file_path: file_path.to_string() })?;
            scldom.file_name = file_name.to_string();
            scldom.file_type = file_type.to_string();
        }

        Ok(scldom)
    }

    pub fn write(&self, pts: &Points, nds: &Nodes, dom: &Domain, tet: &Tetrahedron, ts: usize) -> Result<(), VFVChemError> {
        match self.file_type.as_str() {
            "" => Ok(()),  // no output if blank
            "csv" => write_scldom_csv(&self, pts, nds, dom, ts),
            "vtu" => write_scldom_vtu(&self, pts, nds, dom, tet, ts),
            _ => Err(VFVChemError::ScalarUnsupportedOutputType { geom: "Domain", file_type: self.file_type.clone() })
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
            let mut grad_x = 0.0;
            let mut grad_y = 0.0;
            let mut grad_z = 0.0;
            for (loc, &vid_n) in list_n.iter().enumerate() {
                let (w_x, w_y, w_z) = tet.grad_wgt[vid_c][loc];
                let delta_value = self.value[vid_n] - self.value[vid_c];
                grad_x += w_x * delta_value;
                grad_y += w_y * delta_value;
                grad_z += w_z * delta_value;
            }
            self.grad.push((grad_x, grad_y, grad_z));
        }
    }

    pub fn update_func(vars: &mut Variables, scldom_id: usize, t: f64) {
        // does not use &mut self to prevent borrowing problems
        // assume that scldoms in unkn_id have the same dom as scldom_id

        // pseudocode:
        // collect values from unknown-type scldom (scldom_u)
        // update the selected function-type scldom (scldom_f)

        // compute updated values
        let scldom_f = &vars.scldom[scldom_id];
        match &scldom_f.scldom_type {
            SclDomType::Function { func, unkn_id } => {
                // get number of points
                let num_tot = scldom_f.value.len();
                let dom = &vars.dom[scldom_f.dom_id];

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

                // modify scldom values
                let scldom_f = &mut vars.scldom[scldom_id];
                scldom_f.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
