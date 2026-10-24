use crate::base::prelude::*;

// input: time, slice of unknown values at a point
// output: function value at the point
pub type VecItfFunc = Box<dyn Fn(f64, &[f64], &[f64]) -> (f64, f64, f64) + Send + Sync>;

#[derive(Default)]
pub enum VecItfType {
    #[default] Constant,
    Function {
        func: VecItfFunc,  
        unkn_a_id: Vec<usize>,
        unkn_b_id: Vec<usize>,
    },
}

#[derive(Default)]
pub struct VectorInterface {
    // ids
    pub vecitf_id: usize,
    pub itf_id: usize,  // interface this is attached to

    // data
    pub vecitf_type: VecItfType,
    pub value: Vec<(f64, f64, f64)>,  // idx -> value
    pub value_prev: Vec<(f64, f64, f64)>,
    
    // output file
    pub file_name: String,  // file path without extension
    pub file_type: String,  // file extension
}

impl VectorInterface {
    pub fn new_cnst(vecitf_id: usize, itf: &Interface, value: (f64, f64, f64), file_path: &str) -> Result<VectorInterface, VFVChemError> {
        // default class
        let mut vecitf = VectorInterface::default();
        vecitf.vecitf_id = vecitf_id;
        vecitf.itf_id = itf.itf_id;

        // data
        vecitf.vecitf_type = VecItfType::Constant;
        vecitf.value = vec![value; itf.num_nds()];
        vecitf.value_prev = vecitf.value.clone();

        // output file
        if file_path == "" {
            vecitf.file_name = String::new();
            vecitf.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::VectorUnknownOutputType { geom: "Interface", file_path: file_path.to_string() })?;
            vecitf.file_name = file_name.to_string();
            vecitf.file_type = file_type.to_string();
        }

        Ok(vecitf)
    }

    pub fn new_func(vecitf_id: usize, itf: &Interface, scldom_all: &Vec<ScalarDomain>, func: VecItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>, file_path: &str) -> Result<VectorInterface, VFVChemError> {
        // check that unkn_ids refer to unknown scalars
        for &i in unkn_a_id.iter().chain(unkn_b_id.iter()) {
            match scldom_all[i].scldom_type {
                SclDomType::Unknown { start: _ } => { continue; }
                _ => { return Err(VFVChemError::VectorFunctionDependsOnNonUnknown { geom: "Interface", scldom_id: i }) }
            }
        }

        // default class
        let mut vecitf = VectorInterface::default();
        vecitf.vecitf_id = vecitf_id;
        vecitf.itf_id = itf.itf_id;

        // data
        vecitf.vecitf_type = VecItfType::Function { func, unkn_a_id, unkn_b_id };
        vecitf.value = vec![(f64::NAN, f64::NAN, f64::NAN); itf.num_nds()];  // placeholder
        vecitf.value_prev = vecitf.value.clone();

        // output file
        if file_path == "" {
            vecitf.file_name = String::new();
            vecitf.file_type = String::new();
        } else {
            let (file_name, file_type) = file_path.rsplit_once(".").ok_or(VFVChemError::VectorUnknownOutputType { geom: "Interface", file_path: file_path.to_string() })?;
            vecitf.file_name = file_name.to_string();
            vecitf.file_type = file_type.to_string();
        }

        Ok(vecitf)
    }

    pub fn write(&self, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
        match self.file_type.as_str() {
            "" => Ok(()),  // no output if blank
            "csv" => write_vecitf_csv(&self, nds, itf, ts),
            "vtu" => write_vecitf_vtu(&self, nds, itf, ts),
            _ => Err(VFVChemError::VectorUnsupportedOutputType { geom: "Interface", file_type: self.file_type.clone() })
        }
    }

    pub fn update_prev(&mut self) {
        self.value_prev = self.value.clone();
    }

    pub fn update_func(vars: &mut Variables, vecitf_id: usize, t: f64) {
        // does not use &mut self for similarity with scldom's update_func
        // assume that scldoms in unkn_id have the same dom as vecitf_id
        
        // pseudocode:
        // collect values from unknown-type scldom
        // update the selected function-type vecitf

        // compute updated values
        let vecitf = &vars.vecitf[vecitf_id];
        match &vecitf.vecitf_type {
            VecItfType::Function { func, unkn_a_id, unkn_b_id } => {
                // get objects
                let itf = &vars.itf[vecitf.itf_id];
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
                let vecitf = &mut vars.vecitf[vecitf_id];
                vecitf.value = value_f;

            }
            _ => { return; }  // skip if not function-type
        };

    }


}
