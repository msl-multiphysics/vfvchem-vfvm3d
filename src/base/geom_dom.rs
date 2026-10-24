use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Domain {
    // ids
    pub dom_id: usize,
    pub grp3d_id: String,

    // points in domain
    pub pts_gpt: Vec<usize>,  // pid (point id) -> gpt; local to global
    pub pts_pid: HashMap<usize, usize>,  // gpt -> pid; global to local

    // nodes around domain
    pub nds_gnd: Vec<usize>,  // nid (node id) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> nid; global to local
}

impl Domain {
    pub fn new(dom_id: usize, pts: &Points, nds: &Nodes, grp3d_id: &str) -> Result<Domain, VFVChemError> {
        // input checking
        if !pts.grp3d_gpt.contains_key(grp3d_id) { return Err(VFVChemError::DomainGroupIDNotFound { grp3d_id: grp3d_id.to_string() })}
        if !nds.grp3d_gnd.contains_key(grp3d_id) { return Err(VFVChemError::DomainGroupIDNotFound { grp3d_id: grp3d_id.to_string() })}
        
        // default struct
        let mut dom = Domain::default();
        dom.dom_id = dom_id;
        dom.grp3d_id = grp3d_id.to_string();

        // get points
        dom.pts_gpt = pts.grp3d_gpt[grp3d_id].clone();
        dom.pts_pid = dom.pts_gpt.iter().enumerate().map(|(pid, &gpt)| (gpt, pid)).collect();

        // get nodes
        dom.nds_gnd = nds.grp3d_gnd[grp3d_id].clone();
        dom.nds_nid = dom.nds_gnd.iter().enumerate().map(|(nid, &gnd)| (gnd, nid)).collect();

        Ok(dom)
    }

    #[inline]
    pub fn num_pts(&self) -> usize {
        self.pts_gpt.len()
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }

    #[inline]
    pub fn num_tot(&self) -> usize {
        self.pts_gpt.len() + self.nds_gnd.len()
    }
}
