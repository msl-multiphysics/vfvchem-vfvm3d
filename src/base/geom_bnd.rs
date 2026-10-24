use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Boundary {
    // ids
    pub bnd_id: usize,
    pub dom_id: usize,  // domain this boundary is attached to
    pub grp2d_id: String,
    
    // nodes on boundary
    pub nds_gnd: Vec<usize>,  // idx (index) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> idx; global to local
    pub nds_idx: HashMap<usize, usize>,  // gnd -> idx; only nodes on this boundary
    pub nds_norm: Vec<(f64, f64, f64)>,  // idx -> (nx, ny, nz)
}

impl Boundary {
    pub fn new(bnd_id: usize, nds: &Nodes, dom: &Domain, grp2d_id: &str) -> Result<Boundary, VFVChemError> {
        // input checking
        if !nds.grp2d_gnd.contains_key(grp2d_id) { return Err(VFVChemError::BoundaryGroupIDNotFound { grp2d_id: grp2d_id.to_string() })}

        // default struct
        let mut bnd = Boundary::default();
        bnd.bnd_id = bnd_id;
        bnd.dom_id = dom.dom_id;
        bnd.grp2d_id = grp2d_id.to_string();

        // get nodes
        bnd.nds_gnd = nds.grp2d_gnd[grp2d_id].clone();
        bnd.nds_nid = dom.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();
        bnd.nds_idx = bnd.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();

        // get outward normals
        let nds_norm = &nds.grp3d_norm[&dom.grp3d_id];
        bnd.nds_norm = bnd.nds_gnd.iter().map(|&gnd| nds_norm[bnd.nds_nid[&gnd]]).collect();

        Ok(bnd)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }
}
