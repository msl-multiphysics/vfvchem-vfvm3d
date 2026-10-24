use vfvchem_vfvm3d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state diffusion with first-order reaction 0 -> 1 in a water-filled porous slab.
/// Run with: `cargo run --release --example mass_analytic`
///
/// Geometry:
/// Cube domain defined by lengths along each axis.
/// No flux except on the left, so the solution is the 1D profile c_0(x) = cosh(sqrt(Da) (1 - x / L)) / cosh(sqrt(Da)).
/// Damkohler number Da = k L^2 / D = 0.5.
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
/// dom - diffusion coefficients (m2 s-1)
/// [D_00 D_01] = [1.0e-9 0.0]
/// [D_10 D_11]   [0.0    1.0e-9]
/// dom - production coefficients (s-1; positive if production)
/// [s_0] = [-5.0e-4]
/// [s_1]   [0.0]
/// dom - reaction rates (mol m-3 s-1; positive if production)
/// [R_0] = [0.0]
/// [R_1]   [+5.0e-4 * c_0]
///
/// Boundary conditions:
/// bnd_l - concentration (mol m-3)
/// [c_0] = [1.0]
/// [c_1]   [0.0]
/// bnd_r - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_f - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_k - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_b - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
/// bnd_t - no flux (mol m-2 s-1)
/// [N_0] = [0.0]
/// [N_1]   [0.0]
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_mass_analytic").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, z_min, x_max, y_max, z_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.0, 1e-3, 1e-3, 1e-3, 4e-5, 4e-5, 1e-4, 100, 42)?;

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
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_analytic/conc_0.vtu")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_analytic/conc_1.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    // the consumption of species 0 uses the production coefficient (s_0 = -k)
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-9, "")?;  // D_00
    let diff_11 = vars.add_scldom_cnst(dom, 1.0e-9, "")?;  // D_11
    let mprd_0 = vars.add_scldom_cnst(dom, -5.0e-4, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_1 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass production coefficient (s_1; positive if source; negative if sink)
    let msrc_0 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)

    // non-constant domain values
    // arguments: domain, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y, z), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // species 0 converts to species 1, first order in the reactant
    let msrc_1_func: SclDomFunc = Box::new(|_t, (_x, _y, _z), unkn| { 5.0e-4 * unkn[0] });  // unkn[0] is c_0
    let msrc_1 = vars.add_scldom_func(dom, msrc_1_func, vec![conc_0], "")?;  // mass source (R_1; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let conc_0_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;  // concentration
    let mflx_0_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // molar flux (positive if outward; negative if inward)
    let mflx_0_f = vars.add_sclbnd_cnst(bnd_f, 0.0, "")?;  // molar flux
    let mflx_0_k = vars.add_sclbnd_cnst(bnd_k, 0.0, "")?;  // molar flux
    let mflx_0_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux
    let mflx_0_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // molar flux
    let conc_1_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;  // concentration
    let mflx_1_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // molar flux
    let mflx_1_f = vars.add_sclbnd_cnst(bnd_f, 0.0, "")?;  // molar flux
    let mflx_1_k = vars.add_sclbnd_cnst(bnd_k, 0.0, "")?;  // molar flux
    let mflx_1_b = vars.add_sclbnd_cnst(bnd_b, 0.0, "")?;  // molar flux
    let mflx_1_t = vars.add_sclbnd_cnst(bnd_t, 0.0, "")?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00)]), mprd_0, msrc_0)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), mprd_1, msrc_1)?;  // arguments: component, domain, c, D, s, R
    phys.add_conc_bnd(0, bnd_l, conc_0_l)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(0, bnd_r, mflx_0_r)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_f, mflx_0_f)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_k, mflx_0_k)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_b, mflx_0_b)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_t, mflx_0_t)?;  // arguments: component, boundary, N
    phys.add_conc_bnd(1, bnd_l, conc_1_l)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(1, bnd_r, mflx_1_r)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_f, mflx_1_f)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_k, mflx_1_k)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_b, mflx_1_b)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(1, bnd_t, mflx_1_t)?;  // arguments: component, boundary, N

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // new_with_options - arguments: num_thread, tol, restart, max_outer; tol must be well below the tolerance of the iteration
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearGmres::new_with_options(1, 1e-10, 50, 200)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-8, 1.0)?;
    Ok(())
}
