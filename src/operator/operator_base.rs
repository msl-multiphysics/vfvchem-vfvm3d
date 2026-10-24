use crate::base::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

pub trait OperatorBase {
    // to implement in specific solvers
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError>;

    // matrix insertion convenience functions

    // scalar-scalar blocks
    fn add_sclscl_a(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, row_scl: usize, row_idx: usize, col_scl: usize, col_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.scldom[row_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let col_start = match vars.scldom[col_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let row_xid = row_start + row_idx;
        let col_xid = col_start + col_idx;
        a_triplet.push(Triplet::new(row_xid, col_xid, value));
        Ok(())
    }
    fn add_scl_b(&self, vars: &Variables, b_vec: &mut Col<f64>, row_scl: usize, row_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.scldom[row_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let row_xid = row_start + row_idx;
        b_vec[row_xid] += value;
        Ok(())
    }

    // vector-vector blocks
    fn add_vecvec_a(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, row_vec: usize, row_comp: usize, row_idx: usize, col_vec: usize, col_comp: usize, col_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.vecdom[row_vec].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let col_start = match vars.vecdom[col_vec].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let row_vert = vars.vecdom[row_vec].value.len();
        let col_vert = vars.vecdom[col_vec].value.len();
        let row_xid = row_start + row_comp * row_vert + row_idx;
        let col_xid = col_start + col_comp * col_vert + col_idx;
        a_triplet.push(Triplet::new(row_xid, col_xid, value));
        Ok(())
    }
    fn add_vec_b(&self, vars: &Variables, b_vec: &mut Col<f64>, row_vec: usize, row_comp: usize, row_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.vecdom[row_vec].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let num_pts = vars.vecdom[row_vec].value.len();
        let row_xid = row_start + row_comp * num_pts + row_idx;
        b_vec[row_xid] += value;
        Ok(())
    }

    // scalar-vector blocks
    fn add_sclvec_a(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, row_scl: usize, row_idx: usize, col_vec: usize, col_comp: usize, col_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.scldom[row_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let col_start = match vars.vecdom[col_vec].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let num_pts = vars.vecdom[col_vec].value.len();
        let row_xid = row_start + row_idx;
        let col_xid = col_start + col_comp * num_pts + col_idx;
        a_triplet.push(Triplet::new(row_xid, col_xid, value));
        Ok(())
    }
    fn add_vecscl_a(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, row_vec: usize, row_comp: usize, row_idx: usize, col_scl: usize, col_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.vecdom[row_vec].vecdom_type {
            VecDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let col_start = match vars.scldom[col_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let num_pts = vars.vecdom[row_vec].value.len();
        let row_xid = row_start + row_comp * num_pts + row_idx;
        let col_xid = col_start + col_idx;
        a_triplet.push(Triplet::new(row_xid, col_xid, value));
        Ok(())
    }

}

// input checking convenience functions

// domain
pub fn check_dom_scldom(vars: &Variables, dom_id: usize, scldom_id: usize) -> Result<(), VFVChemError> {
    let scldom = &vars.scldom[scldom_id];
    if scldom.dom_id != dom_id { return Err(VFVChemError::OperatorScalarDomainNotOnDomain { scldom_id, dom_id }) }
    Ok(())
}

// boundary
pub fn check_bnd_scldom(vars: &Variables, bnd_id: usize, scldom_id: usize) -> Result<(), VFVChemError> {
    let scldom = &vars.scldom[scldom_id];
    let bnd = &vars.bnd[bnd_id];
    if scldom.dom_id != bnd.dom_id { return Err(VFVChemError::OperatorScalarDomainNotOnBoundary { scldom_id, bnd_id }) }
    Ok(())
}
pub fn check_bnd_sclbnd(vars: &Variables, bnd_id: usize, sclbnd_id: usize) -> Result<(), VFVChemError> {
    let sclbnd = &vars.sclbnd[sclbnd_id];
    if sclbnd.bnd_id != bnd_id { return Err(VFVChemError::OperatorScalarBoundaryNotOnBoundary { sclbnd_id, bnd_id }) }
    Ok(())
}

// interface
pub fn check_itf_scldom(vars: &Variables, itf_id: usize, scldom_a_id: usize, scldom_b_id: usize) -> Result<(), VFVChemError> {
    let scldom_a = &vars.scldom[scldom_a_id];
    let scldom_b = &vars.scldom[scldom_b_id];
    let itf = &vars.itf[itf_id];
    let order_same = (itf.dom_a_id == scldom_a.dom_id) && (itf.dom_b_id == scldom_b.dom_id);
    let order_swap = (itf.dom_a_id == scldom_b.dom_id) && (itf.dom_b_id == scldom_a.dom_id);
    if !(order_same || order_swap) { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id, scldom_b_id, itf_id }) }
    Ok(())
}
pub fn check_itf_sclitf(vars: &Variables, itf_id: usize, sclitf_id: usize) -> Result<(), VFVChemError> {
    let sclitf = &vars.sclitf[sclitf_id];
    if sclitf.itf_id != itf_id { return Err(VFVChemError::OperatorScalarInterfaceNotOnInterface { sclitf_id, itf_id }) }
    Ok(())
}

// domain
pub fn check_dom_vecdom(vars: &Variables, dom_id: usize, vecdom_id: usize) -> Result<(), VFVChemError> {
    let vecdom = &vars.vecdom[vecdom_id];
    if vecdom.dom_id != dom_id { return Err(VFVChemError::OperatorVectorDomainNotOnDomain { vecdom_id, dom_id }) }
    Ok(())
}

// boundary
pub fn check_bnd_vecdom(vars: &Variables, bnd_id: usize, vecdom_id: usize) -> Result<(), VFVChemError> {
    let vecdom = &vars.vecdom[vecdom_id];
    let bnd = &vars.bnd[bnd_id];
    if vecdom.dom_id != bnd.dom_id { return Err(VFVChemError::OperatorVectorDomainNotOnBoundary { vecdom_id, bnd_id }) }
    Ok(())
}
pub fn check_bnd_vecbnd(vars: &Variables, bnd_id: usize, vecbnd_id: usize) -> Result<(), VFVChemError> {
    let vecbnd = &vars.vecbnd[vecbnd_id];
    if vecbnd.bnd_id != bnd_id { return Err(VFVChemError::OperatorVectorBoundaryNotOnBoundary { vecbnd_id, bnd_id }) }
    Ok(())
}

// interface
pub fn check_itf_vecdom(vars: &Variables, itf_id: usize, vecdom_a_id: usize, vecdom_b_id: usize) -> Result<(), VFVChemError> {
    let vecdom_a = &vars.vecdom[vecdom_a_id];
    let vecdom_b = &vars.vecdom[vecdom_b_id];
    let itf = &vars.itf[itf_id];
    let order_same = (itf.dom_a_id == vecdom_a.dom_id) && (itf.dom_b_id == vecdom_b.dom_id);
    let order_swap = (itf.dom_a_id == vecdom_b.dom_id) && (itf.dom_b_id == vecdom_a.dom_id);
    if !(order_same || order_swap) { return Err(VFVChemError::OperatorVectorDomainNotOnInterface { vecdom_a_id, vecdom_b_id, itf_id }) }
    Ok(())
}
pub fn check_itf_vecitf(vars: &Variables, itf_id: usize, vecitf_id: usize) -> Result<(), VFVChemError> {
    let vecitf = &vars.vecitf[vecitf_id];
    if vecitf.itf_id != itf_id { return Err(VFVChemError::OperatorVectorInterfaceNotOnInterface { vecitf_id, itf_id }) }
    Ok(())
}

// domain
pub fn check_dom_sclflw(vars: &Variables, dom_id: usize, sclflw_id: usize) -> Result<(), VFVChemError> {
    let sclflw = &vars.sclflw[sclflw_id];
    if sclflw.dom_id != dom_id { return Err(VFVChemError::OperatorScalarFlowNotOnDomain { sclflw_id, dom_id }) }
    Ok(())
}

// boundary
pub fn check_bnd_sclflw(vars: &Variables, bnd_id: usize, sclflw_id: usize) -> Result<(), VFVChemError> {
    let sclflw = &vars.sclflw[sclflw_id];
    let bnd = &vars.bnd[bnd_id];
    if sclflw.dom_id != bnd.dom_id { return Err(VFVChemError::OperatorScalarFlowNotOnBoundary { sclflw_id, bnd_id }) }
    Ok(())
}
