use vfvchem_vfvm3d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Production coefficient versus mass source for a fast first-order consumption of species 0.
/// Run with: `cargo run --release --example feat_prod`
///
/// The consumption is given either as a production coefficient (s_0 = -k) or as a mass source (R_0 = -k c_0).
/// The production coefficient is implicit (s_0 c_0 is loaded into the matrix), so it converges in a few iterations.
/// The mass source is lagged by one iteration, so the iteration diverges at large Damkohler numbers.
///
/// Geometry:
/// Cube domain defined by lengths along each axis.
/// Damkohler number Da = k L^2 / D = 20; 1D analytic profile c_0(x) = cosh(sqrt(Da) (1 - x / L)) / cosh(sqrt(Da)).
///
///             bnd_t   bnd_k
///            +---------+
///           /         /|
///          +---------+ |
///    bnd_l |         | + bnd_r
///          |   dom   |/
///          +---------+
///         bnd_f   bnd_b
///
/// Domains:
/// dom - main domain (1.0e-3 m x 1.0e-3 m x 1.0e-3 m)
///
/// Boundaries:
/// bnd_l - left boundary (x = 0.0 m)
/// bnd_r - right boundary (x = 1.0e-3 m)
/// bnd_f - front boundary (y = 0.0 m)
/// bnd_k - back boundary (y = 1.0e-3 m)
/// bnd_b - bottom boundary (z = 0.0 m)
/// bnd_t - top boundary (z = 1.0e-3 m)
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - diffusion coefficient (D_00 = 1.0e-9 m2 s-1)
/// dom - production coefficient (s_0 = -0.02 s-1 if production; 0.0 s-1 if source)
/// dom - reaction rate (R_0 = 0.0 mol m-3 s-1 if production; -0.02 * c_0 mol m-3 s-1 if source)
///
/// Boundary conditions:
/// bnd_l - concentration (c_0 = 1.0 mol m-3)
/// bnd_r - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_f - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_k - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_b - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_t - no flux (N_0 = 0.0 mol m-2 s-1)
///
fn run(use_prd: bool) -> Result<(), VFVChemError> {
    // output directory
    let out = if use_prd { "examples/output_feat_prod/production" } else { "examples/output_feat_prod/source" };
    create_dir_all(out).unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, z_min, x_max, y_max, z_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.0, 1e-3, 1e-3, 1e-3, 5e-5, 5e-5, 5e-5, 100, 42)?;

    // geometry
    // if using new_bounds, the domain is "dom"
    // the left, right, front, back, bottom, and top boundaries are "bnd_l", "bnd_r", "bnd_f", "bnd_k", "bnd_b", and "bnd_t"
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let bnd_f = vars.add_bnd(dom, "bnd_f")?;
    let bnd_k = vars.add_bnd(dom, "bnd_k")?;
    let bnd_b = vars.add_bnd(dom, "bnd_b")?;
    let bnd_t = vars.add_bnd(dom, "bnd_t")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., conc_0_0000.vtu)
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, &format!("{out}/conc_0.vtu"))?;

    // constant domain values
    // arguments: domain, value, output_file
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-9, "")?;  // D_00

    // non-constant domain values
    // arguments: domain, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y, z), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // species 0 is consumed with the rate constant k = Da D / L^2 = 0.02 s-1
    // the consumption is either a production coefficient (s_0 = -k) or a mass source (R_0 = -k c_0)
    let (mprd_0, msrc_0) = if use_prd {
        let mprd_0 = vars.add_scldom_cnst(dom, -0.02, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
        let msrc_0 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)
        (mprd_0, msrc_0)
    } else {
        let msrc_0_func: SclDomFunc = Box::new(|_t, (_x, _y, _z), unkn| { -0.02 * unkn[0] });  // unkn[0] is c_0
        let mprd_0 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
        let msrc_0 = vars.add_scldom_func(dom, msrc_0_func, vec![conc_0], "")?;  // mass source (R_0; positive if source; negative if sink)
        (mprd_0, msrc_0)
    };

    // constant boundary values
    // arguments: boundary, value, output_file
    let conc_0_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;  // concentration
    let mflx_0_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // molar flux (positive if outward; negative if inward)
    let mflx_0_f = vars.add_sclbnd_cnst(bnd_f, 0.0, "")?;  // molar flux
    let mflx_0_k = vars.add_sclbnd_cnst(bnd_k, 0.0, "")?;  // molar flux
    let mflx_0_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux
    let mflx_0_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(1)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00)]), mprd_0, msrc_0)?;  // arguments: component, domain, c, D, s, R
    phys.add_conc_bnd(0, bnd_l, conc_0_l)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(0, bnd_r, mflx_0_r)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_f, mflx_0_f)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_k, mflx_0_k)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_b, mflx_0_b)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_t, mflx_0_t)?;  // arguments: component, boundary, N

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // new_with_options - arguments: num_thread, tol, restart, max_outer; tol must be well below the tolerance of the iteration
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    let mut lin = LinearGmres::new_with_options(1, 1e-10, 50, 200)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-8, 1.0)?;

    Ok(())
}

fn main() -> Result<(), VFVChemError> {
    // production coefficient; converges in a few iterations
    run(true)?;

    // mass source; the iteration is expected to diverge at this Damkohler number
    // the linear solver may stop first as the iteration diverges
    match run(false) {
        Err(VFVChemError::SteadyDidNotConverge { .. }) | Err(VFVChemError::LinearGmresDidNotConverge { .. }) => {},
        res => res?,
    }

    Ok(())
}
