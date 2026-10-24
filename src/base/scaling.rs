use crate::base::prelude::*;
use faer::Col;

// Per-field scaling of convergence measures (Picard residual and adaptive error norm).
//
// Every unknown field (scalar domain or vector domain) is measured against its own reference scale, and the
// measure is the maximum over fields, so a field of large magnitude (e.g., T ~ 300 K) does not hide the others.
//   Picard residual: max_f rms_f(x_new - x_old) / s_f
// Reference scale s_f (largest of the candidates; a value set with set_scldom_scale / set_vecdom_scale is a lower bound):
//   scalar   : max |v| over both states
//   vector   : max |(v_x, v_y, v_z)| over both states (all components share one scale)
//   pressure : max(range of p, rho_max * U_max^2) over both states; the absolute value is not used, so a large
//              constant offset (absolute pressure) does not loosen the measure; at near-zero pressure differences
//              the user reference pressure difference (set_scldom_scale) or the dynamic pressure is the floor
// A field is a pressure if it is the pressure of a flow (sclflw); its velocity and density come from the same flow.
// A tiny absolute floor (SCALE_FLOOR) avoids division by zero for fields that are identically zero.
// Pressure round-off: with a large offset (absolute pressure), the linear solve resolves p only to about
// PRES_ROUNDOFF * max|p| (observed ~3e-13 relative), which can exceed tol * range(p). Pressure changes below
// this level are not counted (dead band) in the Picard residual and are added to atol in the error norm, so the
// iteration does not stagnate on round-off; at 1e5 Pa this is 1e-6 Pa, far below any resolved pressure difference.

pub const SCALE_FLOOR: f64 = 1e-30;
pub const PRES_ROUNDOFF: f64 = 1e-11;

pub fn pres_roundoff(field: &FieldNorm, x_a: &Col<f64>, x_b: &Col<f64>) -> f64 {
    // absolute pressure round-off level PRES_ROUNDOFF * max|p| (0 for other fields)
    match field.kind {
        FieldKind::Pressure { .. } => {
            let mut pabs: f64 = 0.0;
            for x in [x_a, x_b] { for i in field.start..field.start + field.len { pabs = pabs.max(x[i].abs()); } }
            PRES_ROUNDOFF * pabs
        }
        _ => 0.0,
    }
}

pub enum FieldKind {
    Scalar,
    Vector,
    Pressure { vel_id: usize, den_id: usize },
}

pub struct FieldNorm {
    pub name: String,   // e.g. "scldom 0" or "vecdom 0", for output
    pub id: usize,      // scldom or vecdom id
    pub start: usize,   // start in solution vector
    pub len: usize,     // number of points (per component for vectors)
    pub kind: FieldKind,
}

pub fn field_norms(vars: &Variables, exclude: &[usize]) -> Vec<FieldNorm> {
    // unknown fields in the solution vector; scldom ids in exclude are left out
    let mut fields = Vec::new();
    for (id, scldom) in vars.scldom.iter().enumerate() {
        if exclude.contains(&id) { continue; }
        if let SclDomType::Unknown { start } = scldom.scldom_type {
            let kind = match vars.sclflw.iter().find(|f| f.pres_id == id) {
                Some(f) => FieldKind::Pressure { vel_id: f.vel_id, den_id: f.den_id },
                None => FieldKind::Scalar,
            };
            fields.push(FieldNorm { name: format!("scldom {id}"), id, start, len: scldom.value.len(), kind });
        }
    }
    for (id, vecdom) in vars.vecdom.iter().enumerate() {
        if let VecDomType::Unknown { start } = vecdom.vecdom_type {
            fields.push(FieldNorm { name: format!("vecdom {id}"), id, start, len: vecdom.value.len(), kind: FieldKind::Vector });
        }
    }
    fields
}

fn _vec_max(vars: &Variables, vecdom_id: usize, x_a: &Col<f64>, x_b: &Col<f64>) -> f64 {
    // max |(v_x, v_y, v_z)| over both states
    let vecdom = &vars.vecdom[vecdom_id];
    let mut vmax: f64 = 0.0;
    match vecdom.vecdom_type {
        VecDomType::Unknown { start } => {
            let n = vecdom.value.len();
            for x in [x_a, x_b] {
                for i in 0..n { vmax = vmax.max((x[start + i].powi(2) + x[start + n + i].powi(2) + x[start + 2 * n + i].powi(2)).sqrt()); }
            }
        }
        _ => { for v in &vecdom.value { vmax = vmax.max((v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()); } }
    }
    vmax
}

pub fn field_scale(vars: &Variables, field: &FieldNorm, x_a: &Col<f64>, x_b: &Col<f64>) -> f64 {
    // reference scale of a field from two states (see the header)
    let (s, n) = (field.start, field.len);
    let scale = match field.kind {
        FieldKind::Scalar => {
            let mut vmax: f64 = 0.0;
            for x in [x_a, x_b] { for i in s..s + n { vmax = vmax.max(x[i].abs()); } }
            vmax.max(vars.scale_scldom.get(&field.id).copied().unwrap_or(0.0))
        }
        FieldKind::Vector => {
            _vec_max(vars, field.id, x_a, x_b).max(vars.scale_vecdom.get(&field.id).copied().unwrap_or(0.0))
        }
        FieldKind::Pressure { vel_id, den_id } => {
            let mut pmin = f64::MAX;
            let mut pmax = f64::MIN;
            for x in [x_a, x_b] { for i in s..s + n { pmin = pmin.min(x[i]); pmax = pmax.max(x[i]); } }
            let range = if n > 0 { pmax - pmin } else { 0.0 };
            let rho = vars.scldom[den_id].value.iter().fold(0.0_f64, |a, &v| a.max(v.abs()));
            let vel = _vec_max(vars, vel_id, x_a, x_b);
            range.max(rho * vel * vel).max(vars.scale_scldom.get(&field.id).copied().unwrap_or(0.0))
        }
    };
    scale.max(SCALE_FLOOR)
}

pub fn scaled_residual(vars: &Variables, fields: &[FieldNorm], x_old: &Col<f64>, x_new: &Col<f64>) -> (f64, usize) {
    // Picard residual: max over fields of rms(x_new - x_old) / s_f; returns the residual and the index of the largest field
    // with scale_off, the unscaled residual |x_new - x_old| / |x_old| (field index 0)
    // NaN if x_new is not finite; f64::max drops NaN, so the per-field sums below would hide it
    if (0..x_new.nrows()).any(|i| !x_new[i].is_finite()) {
        return (f64::NAN, 0);
    }
    if vars.scale_off || fields.is_empty() {
        return ((x_new - x_old).norm_l2() / (x_old.norm_l2() + 1e-10), 0);
    }
    let mut res: f64 = 0.0;
    let mut fid = 0;
    for (k, field) in fields.iter().enumerate() {
        let (s, n) = (field.start, field.len);
        let ncomp = if let FieldKind::Vector = field.kind { 3 } else { 1 };
        let dead = pres_roundoff(field, x_old, x_new);  // pressure round-off dead band (0 for other fields)
        let mut sum = 0.0;
        for i in s..s + ncomp * n { sum += ((x_new[i] - x_old[i]).abs() - dead).max(0.0).powi(2); }
        let rms = (sum / (n.max(1) as f64)).sqrt();
        let r = rms / field_scale(vars, field, x_old, x_new);
        if r > res || r.is_nan() { res = r; fid = k; }
    }
    (res, fid)
}
