use crate::base::prelude::*;

// input: time, (x, y, z), slice of unknown values at a point
// output: function value at the point
pub type SclBndFunc = Box<dyn Fn(f64, (f64, f64, f64), &[f64]) -> f64 + Send + Sync>;

#[derive(Default)]
pub enum SclBndType {
    #[default] Constant,
    Function {
        func: SclBndFunc,  
        unkn_id: Vec<usize>,
    },
}

#[derive(Default)]
pub struct ScalarBoundary {
    // ids
    pub sclbnd_id: usize,
    pub bnd_id: usize,  // boundary this is attached to

    // node data
    pub sclbnd_type: SclBndType,
    pub value: Vec<f64>,  // idx -> value
    pub value_prev: Vec<f64>,

    // output file
    pub file_name: String,  // file path without extension
    pub file_type: String,  // file extension
}

impl ScalarBoundary {
    pub fn new_cnst(sclbnd_id: usize, bnd: &Boundary, value: f64, file_path: &str) -> Result<ScalarBoundary, VFVChemError> {
        // default class
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // node data
        sclbnd.sclbnd_type = SclBndType::Constant;
        sclbnd.value = vec![value; bnd.num_nds()];
        sclbnd.value_prev = sclbnd.value.clone();
        
        // output file
        if file_path == "" {
            sclbnd.file_name = String::new();
            sclbnd.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Boundary", file_path: file_path.to_string() })?;
            sclbnd.file_name = file_name.to_string();
            sclbnd.file_type = file_type.to_string();
        }

        Ok(sclbnd)
    }

    pub fn new_func(sclbnd_id: usize, bnd: &Boundary, scldom_all: &Vec<ScalarDomain>, func: SclBndFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<ScalarBoundary, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_id.iter() {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::ScalarFunctionDependsOnNonUnknown { geom: "Boundary", scldom_id: i }) }
            }
        }

        // default class
        let mut sclbnd = ScalarBoundary::default();
        sclbnd.sclbnd_id = sclbnd_id;
        sclbnd.bnd_id = bnd.bnd_id;

        // node data
        sclbnd.sclbnd_type = SclBndType::Function { func, unkn_id };
        sclbnd.value = vec![f64::NAN; bnd.num_nds()];  // placeholder
        sclbnd.value_prev = sclbnd.value.clone();
        
        // output file
        if file_path == "" {
            sclbnd.file_name = String::new();
            sclbnd.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Boundary", file_path: file_path.to_string() })?;
            sclbnd.file_name = file_name.to_string();
            sclbnd.file_type = file_type.to_string();
        }

        Ok(sclbnd)
    }

    pub fn write(&self, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
        match self.file_type.as_str() {
            "" => Ok(()),  // no output if blank
            "csv" => write_sclbnd_csv(&self, nds, bnd, ts),
            "vtu" => write_sclbnd_vtu(&self, nds, bnd, ts),
            _ => Err(VFVChemError::ScalarUnsupportedOutputType { geom: "Boundary", file_type: self.file_type.clone() })
        }
    }

    pub fn update_prev(&mut self) {
        self.value_prev = self.value.clone();
    }

    pub fn update_func(vars: &mut Variables, sclbnd_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as sclbnd_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type sclbnd

        // compute updated values
        let sclbnd = &vars.sclbnd[sclbnd_id];
        match &sclbnd.sclbnd_type {
            SclBndType::Function { func, unkn_id } => {
                // get objects
                let bnd = &vars.bnd[sclbnd.bnd_id];
                let dom = &vars.dom[bnd.dom_id];

                // collect positions and values from unknown scalars
                let mut pos = Vec::with_capacity(bnd.num_nds());  // idx -> (x, y, z)
                let mut value_u = Vec::with_capacity(bnd.num_nds());  // idx -> values of unknowns
                for idx in 0..bnd.num_nds() {
                    // get domain index
                    let gnd = bnd.nds_gnd[idx];
                    let nid = dom.nds_nid[&gnd];
                    let vid = dom.num_pts() + nid;

                    // collect position and values
                    pos.push(vars.nds.nds_pos[gnd]);
                    let mut value_sub = Vec::with_capacity(unkn_id.len());  // values of unknowns
                    for &i in unkn_id {
                        let scldom = &vars.scldom[i];
                        value_sub.push(scldom.value[vid]);
                    }
                    value_u.push(value_sub);
                }

                // compute function value
                let mut value_f = Vec::with_capacity(bnd.num_nds());
                for idx in 0..bnd.num_nds() {
                    let value_sub = func(t, pos[idx], value_u[idx].as_slice());
                    value_f.push(value_sub)
                }

                // modify scldom values
                let sclbnd = &mut vars.sclbnd[sclbnd_id];
                sclbnd.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
