use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Interface {
    // ids
    pub itf_id: usize,
    pub dom_a_id: usize,  // domains this interface is attached to
    pub dom_b_id: usize,
    pub grp2d_id: String,
    
    // nodes on interface
    pub nds_gnd: Vec<usize>,  // idx (index) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> idx; global to local
    pub nds_norm_ab: Vec<(f64, f64, f64)>,  // idx -> (nx, ny, nz) from dom a to b
}

impl Interface {
    pub fn new(itf_id: usize, nds: &Nodes, dom_a: &Domain, dom_b: &Domain, grp2d_id: &str) -> Result<Interface, VFVChemError> {
        // input checking
        if !nds.grp2d_gnd.contains_key(grp2d_id) { return Err(VFVChemError::InterfaceGroupIDNotFound { grp2d_id: grp2d_id.to_string() })}
        
        // check that grp2d_id nodes are a subset of dom_a and dom_b
        let grp2d_nds = &nds.grp2d_gnd[grp2d_id];
        let not_subset = !grp2d_nds.iter().all(|gnd| {dom_a.nds_nid.contains_key(gnd) && dom_b.nds_nid.contains_key(gnd)});
        if not_subset { return Err(VFVChemError::InterfaceNoCommonInterface { dom_a: dom_a.dom_id, dom_b: dom_b.dom_id, grp2d_id: grp2d_id.to_string() })}

        // default struct
        let mut itf = Interface::default();
        itf.itf_id = itf_id;
        itf.dom_a_id = dom_a.dom_id;
        itf.dom_b_id = dom_b.dom_id;
        itf.grp2d_id = grp2d_id.to_string();

        // get nodes
        itf.nds_gnd = nds.grp2d_gnd[grp2d_id].clone();
        itf.nds_nid = dom_a.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();

        // get outward normals
        let nds_norm_ab = &nds.grp3d_norm[&dom_a.grp3d_id];
        itf.nds_norm_ab = itf.nds_gnd.iter().map(|&gnd| nds_norm_ab[itf.nds_nid[&gnd]]).collect();

        Ok(itf)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }
}
