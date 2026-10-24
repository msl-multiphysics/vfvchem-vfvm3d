use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomPressureReference {
    // domain
    pub dom_id: usize,

    // scalars
    pub pres_id: usize,  // pressure (p); unknown

    // pressure reference
    pub nid: usize,  // node where the pressure is set
    pub pres_ref: f64,   // reference pressure (p_ref)
}

impl OpSclDomPressureReference {
    pub fn new(vars: &Variables, dom_id: usize, pres_id: usize, nid: usize, pres_ref: f64) -> Result<OpSclDomPressureReference, VFVChemError> {
        // sets the pressure at one node
        // p = p_ref
        //
        // must be applied after all other operators
        // interface nodes carry interface conditions and cannot be modified

        // input checking
        check_dom_scldom(vars, dom_id, pres_id)?;
        if nid >= vars.dom[dom_id].num_nds() { return Err(VFVChemError::OperatorNodeNotFound { dom_id, nid }) }
        let gnd = vars.dom[dom_id].nds_gnd[nid];
        for itf in &vars.itf {
            let on_dom = itf.dom_a_id == dom_id || itf.dom_b_id == dom_id;
            if on_dom && itf.nds_gnd.contains(&gnd) { return Err(VFVChemError::OperatorNodeOnInterface { dom_id, nid, itf_id: itf.itf_id }) }
        }

        // create struct
        let mut oper = OpSclDomPressureReference::default();
        oper.dom_id = dom_id;
        oper.pres_id = pres_id;
        oper.nid = nid;
        oper.pres_ref = pres_ref;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomPressureReference {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get row in the global system
        let dom = &vars.dom[self.dom_id];
        let start = match vars.scldom[self.pres_id].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let xid = start + dom.num_pts() + self.nid;

        // replace the existing equation -> load into matrix
        a_triplet.retain(|t| t.row != xid);
        a_triplet.push(Triplet::new(xid, xid, 1.0));
        b_vec[xid] = self.pres_ref;

        Ok(())
    }
}
