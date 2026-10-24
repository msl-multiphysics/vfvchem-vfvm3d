use crate::base::prelude::*;

// input: time, slice of unknown values at a point
// output: function value at the point
pub type SclItfFunc = Box<dyn Fn(f64, &[f64], &[f64]) -> f64 + Send + Sync>;

#[derive(Default)]
pub enum SclItfType {
    #[default] Constant,
    Function {
        func: SclItfFunc,  
        unkn_a_id: Vec<usize>,
        unkn_b_id: Vec<usize>,
    },
}

#[derive(Default)]
pub struct ScalarInterface {
    // ids
    pub sclitf_id: usize,
    pub itf_id: usize,  // interface this is attached to

    // node data
    pub sclitf_type: SclItfType,
    pub value: Vec<f64>,  // idx -> value
    pub value_prev: Vec<f64>,

    // output file
    pub file_name: String,  // file path without extension
    pub file_type: String,  // file extension
}

impl ScalarInterface {
    pub fn new_cnst(sclitf_id: usize, itf: &Interface, value: f64, file_path: &str) -> Result<ScalarInterface, VFVChemError> {
        // default class
        let mut sclitf = ScalarInterface::default();
        sclitf.sclitf_id = sclitf_id;
        sclitf.itf_id = itf.itf_id;

        // node data
        sclitf.sclitf_type = SclItfType::Constant;
        sclitf.value = vec![value; itf.num_nds()];
        sclitf.value_prev = sclitf.value.clone();
        
        // output file
        if file_path == "" {
            sclitf.file_name = String::new();
            sclitf.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Interface", file_path: file_path.to_string() })?;
            sclitf.file_name = file_name.to_string();
            sclitf.file_type = file_type.to_string();
        }

        Ok(sclitf)
    }

    pub fn new_func(sclitf_id: usize, itf: &Interface, scldom_all: &Vec<ScalarDomain>, func: SclItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>, file_path: &str) -> Result<ScalarInterface, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_a_id.iter().chain(unkn_b_id.iter()) {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::ScalarFunctionDependsOnNonUnknown { geom: "Interface", scldom_id: i }) }
            }
        }

        // default class
        let mut sclitf = ScalarInterface::default();
        sclitf.sclitf_id = sclitf_id;
        sclitf.itf_id = itf.itf_id;

        // node data
        sclitf.sclitf_type = SclItfType::Function { func, unkn_a_id, unkn_b_id };
        sclitf.value = vec![f64::NAN; itf.num_nds()];  // placeholder
        sclitf.value_prev = sclitf.value.clone();
        
        // output file
        if file_path == "" {
            sclitf.file_name = String::new();
            sclitf.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::ScalarUnknownOutputType { geom: "Interface", file_path: file_path.to_string() })?;
            sclitf.file_name = file_name.to_string();
            sclitf.file_type = file_type.to_string();
        }

        Ok(sclitf)
    }

    pub fn write(&self, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
        match self.file_type.as_str() {
            "" => Ok(()),  // no output if blank
            "csv" => write_sclitf_csv(&self, nds, itf, ts),
            "vtu" => write_sclitf_vtu(&self, nds, itf, ts),
            _ => Err(VFVChemError::ScalarUnsupportedOutputType { geom: "Interface", file_type: self.file_type.clone() })
        }
    }

    pub fn update_prev(&mut self) {
        self.value_prev = self.value.clone();
    }

    pub fn update_func(vars: &mut Variables, sclitf_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as sclitf_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type sclitf

        // compute updated values
        let sclitf = &vars.sclitf[sclitf_id];
        match &sclitf.sclitf_type {
            SclItfType::Function { func, unkn_a_id, unkn_b_id } => {
                // get objects
                let itf = &vars.itf[sclitf.itf_id];
                let dom_a = &vars.dom[itf.dom_a_id];
                let dom_b = &vars.dom[itf.dom_b_id];

                // collect values from unknown scalars
                let mut value_a = Vec::with_capacity(itf.num_nds());  // idx -> values of unknowns
                let mut value_b = Vec::with_capacity(itf.num_nds());
                for idx in 0..itf.num_nds() {
                    // get domain index
                    let gnd = itf.nds_gnd[idx];
                    let nid_a = dom_a.nds_nid[&gnd];
                    let vid_a = dom_a.num_pts() + nid_a;
                    let nid_b = dom_b.nds_nid[&gnd];
                    let vid_b = dom_b.num_pts() + nid_b;

                    // collect values
                    let mut value_a_sub = Vec::with_capacity(unkn_a_id.len());  // values of unknowns
                    for &i in unkn_a_id {
                        let scldom_a = &vars.scldom[i];
                        value_a_sub.push(scldom_a.value[vid_a]);
                    }
                    value_a.push(value_a_sub);
                    let mut value_b_sub = Vec::with_capacity(unkn_b_id.len());  // values of unknowns
                    for &i in unkn_b_id {
                        let scldom_b = &vars.scldom[i];
                        value_b_sub.push(scldom_b.value[vid_b]);
                    }
                    value_b.push(value_b_sub);
                }

                // compute function value
                let mut value_f = Vec::with_capacity(itf.num_nds());
                for idx in 0..itf.num_nds() {
                    let value_sub = func(t, value_a[idx].as_slice(), value_b[idx].as_slice());
                    value_f.push(value_sub)
                }

                // modify scldom values
                let sclitf = &mut vars.sclitf[sclitf_id];
                sclitf.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }

}
