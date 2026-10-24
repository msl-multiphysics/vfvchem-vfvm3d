use vfvchem_vfvm3d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Heat and mass transfer through three regions with different properties joined by interfaces.
/// Run with: `cargo run --release --example feat_itf`
///
/// The interface conditions make the unknown and its flux continuous across the interface.
/// Heat: water in the corner regions and glass in the middle, joined with add_hcnt_itf.
/// Mass: two diffusion coefficients, joined with add_mcnt_itf.
///
/// Geometry:
/// Cube with a lower-front-left and an upper-back-right corner cuboid; see examples/obj/obj_threereg.py.
///
///                  +------+------+
///                 /       /dom_up|
///                +-------+------+|
///               /              / |
///              +--------------+  |
///              |    dom_m     |  +
///              |+------+      | /
///              ||dom_lo|      |/
///              ++------+------+
///
/// Domains:
/// dom_lo - lower-front-left region (0.75 m x 0.25 m x 0.25 m)
/// dom_m - middle region (1.0 m x 1.0 m x 1.0 m cube minus the corner regions)
/// dom_up - upper-back-right region (0.75 m x 0.25 m x 0.25 m)
///
/// Boundaries:
/// bnd_lo - outer faces of dom_lo (left, front, and bottom)
/// bnd_mb - outer faces of dom_m (bottom, front, and right)
/// bnd_mt - outer faces of dom_m (top, back, and left)
/// bnd_up - outer faces of dom_up (right, back, and top)
///
/// Interfaces:
/// itf_lo - between dom_lo and dom_m
/// itf_up - between dom_up and dom_m
///
/// Heat Transfer
///
/// Properties:
/// dom_lo, dom_up - thermal conductivity (k = 0.6 W m-1 K-1)
/// dom_m - thermal conductivity (k = 1.4 W m-1 K-1)
/// dom_lo, dom_m, dom_up - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_lo, dom_m, dom_up - heat source (Q = 0.0 W m-3)
///
/// Boundary conditions:
/// bnd_lo - temperature (T = 350 K)
/// bnd_up - temperature (T = 300 K)
/// bnd_mb - no flux (q = 0.0 W m-2)
/// bnd_mt - no flux (q = 0.0 W m-2)
///
/// Interface conditions:
/// itf_lo - temperature and flux continuity
/// itf_up - temperature and flux continuity
///
/// Mass Transfer
///
/// Properties:
/// dom_lo, dom_up - diffusion coefficient (D_00 = 1.0e-9 m2 s-1)
/// dom_m - diffusion coefficient (D_00 = 4.0e-9 m2 s-1)
/// dom_lo, dom_m, dom_up - production coefficient (s_0 = 0.0 s-1)
/// dom_lo, dom_m, dom_up - reaction rate (R_0 = 0.0 mol m-3 s-1)
///
/// Boundary conditions:
/// bnd_lo - concentration (c_0 = 1.0 mol m-3)
/// bnd_up - concentration (c_0 = 0.0 mol m-3)
/// bnd_mb - no flux (N_0 = 0.0 mol m-2 s-1)
/// bnd_mt - no flux (N_0 = 0.0 mol m-2 s-1)
///
/// Interface conditions:
/// itf_lo - concentration and flux continuity
/// itf_up - concentration and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_feat_itf").unwrap();

    // heat transfer

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain surfaces from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_threereg.obj", 0.05, 0.1, 0.2, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    // add_itf - arguments: domain_a, domain_b, group_name
    let dom_lo = vars.add_dom("dom_lo")?;  // lower-front-left
    let dom_m = vars.add_dom("dom_m")?;  // middle
    let dom_up = vars.add_dom("dom_up")?;  // upper-back-right
    let bnd_lo = vars.add_bnd(dom_lo, "bnd_lo")?;  // left, front, and bottom of dom_lo
    let bnd_mb = vars.add_bnd(dom_m, "bnd_mb")?;  // bottom, front, and right of dom_m
    let bnd_mt = vars.add_bnd(dom_m, "bnd_mt")?;  // top, back, and left of dom_m
    let bnd_up = vars.add_bnd(dom_up, "bnd_up")?;  // right, back, and top of dom_up
    let itf_lo = vars.add_itf(dom_lo, dom_m, "itf_lo")?;  // lower interface
    let itf_up = vars.add_itf(dom_up, dom_m, "itf_up")?;  // upper interface

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_lo_0000.vtu)
    let temp_lo = vars.add_scldom_unkn(dom_lo, 300.0, "examples/output_feat_itf/temp_lo.vtu")?;
    let temp_m = vars.add_scldom_unkn(dom_m, 300.0, "examples/output_feat_itf/temp_m.vtu")?;
    let temp_up = vars.add_scldom_unkn(dom_up, 300.0, "examples/output_feat_itf/temp_up.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let cond_lo = vars.add_scldom_cnst(dom_lo, 0.6, "")?;  // thermal conductivity (water)
    let cond_m = vars.add_scldom_cnst(dom_m, 1.4, "")?;  // thermal conductivity (glass)
    let cond_up = vars.add_scldom_cnst(dom_up, 0.6, "")?;  // thermal conductivity (water)
    let hprd_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_m = vars.add_scldom_cnst(dom_m, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_m = vars.add_scldom_cnst(dom_m, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_bl = vars.add_sclbnd_cnst(bnd_lo, 350.0, "")?;  // temperature
    let temp_bu = vars.add_sclbnd_cnst(bnd_up, 300.0, "")?;  // temperature
    let hflx_mb = vars.add_sclbnd_cnst(bnd_mb, 0.0, "")?;  // heat flux (positive if outward; negative if inward)
    let hflx_mt = vars.add_sclbnd_cnst(bnd_mt, 0.0, "")?;  // heat flux

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hflx_bnd - register boundary with heat flux
    // add_hcnt_itf - register interface with temperature and flux continuity
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_lo, temp_lo, cond_lo, hprd_lo, hsrc_lo)?;  // arguments: domain, T, k, s, Q
    phys.add_heat_dom(dom_m, temp_m, cond_m, hprd_m, hsrc_m)?;  // arguments: domain, T, k, s, Q
    phys.add_heat_dom(dom_up, temp_up, cond_up, hprd_up, hsrc_up)?;  // arguments: domain, T, k, s, Q
    phys.add_temp_bnd(bnd_lo, temp_bl)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_up, temp_bu)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_mb, hflx_mb)?;  // arguments: boundary, q
    phys.add_hflx_bnd(bnd_mt, hflx_mt)?;  // arguments: boundary, q
    phys.add_hcnt_itf(itf_lo)?;  // arguments: interface
    phys.add_hcnt_itf(itf_up)?;  // arguments: interface

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // new_with_options - arguments: num_thread, tol, restart, max_outer; tol must be well below the tolerance of the iteration
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    let mut lin = LinearGmres::new_with_options(1, 1e-10, 50, 200)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    // mass transfer
    // each problem uses its own variables

    // nodes, points, and variables
    let mut vars = Variables::new_obj("examples/obj/obj_threereg.obj", 0.05, 0.1, 0.2, 100, 42)?;

    // geometry
    let dom_lo = vars.add_dom("dom_lo")?;  // lower-front-left
    let dom_m = vars.add_dom("dom_m")?;  // middle
    let dom_up = vars.add_dom("dom_up")?;  // upper-back-right
    let bnd_lo = vars.add_bnd(dom_lo, "bnd_lo")?;  // left, front, and bottom of dom_lo
    let bnd_mb = vars.add_bnd(dom_m, "bnd_mb")?;  // bottom, front, and right of dom_m
    let bnd_mt = vars.add_bnd(dom_m, "bnd_mt")?;  // top, back, and left of dom_m
    let bnd_up = vars.add_bnd(dom_up, "bnd_up")?;  // right, back, and top of dom_up
    let itf_lo = vars.add_itf(dom_lo, dom_m, "itf_lo")?;  // lower interface
    let itf_up = vars.add_itf(dom_up, dom_m, "itf_up")?;  // upper interface

    // unknown domain values
    let conc_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_feat_itf/conc_lo.vtu")?;
    let conc_m = vars.add_scldom_unkn(dom_m, 0.0, "examples/output_feat_itf/conc_m.vtu")?;
    let conc_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_feat_itf/conc_up.vtu")?;

    // constant domain values
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_lo = vars.add_scldom_cnst(dom_lo, 1.0e-9, "")?;  // D_00
    let diff_m = vars.add_scldom_cnst(dom_m, 4.0e-9, "")?;  // D_00
    let diff_up = vars.add_scldom_cnst(dom_up, 1.0e-9, "")?;  // D_00
    let mprd_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_m = vars.add_scldom_cnst(dom_m, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let msrc_lo = vars.add_scldom_cnst(dom_lo, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_m = vars.add_scldom_cnst(dom_m, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_up = vars.add_scldom_cnst(dom_up, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)

    // constant boundary values
    let conc_bl = vars.add_sclbnd_cnst(bnd_lo, 1.0, "")?;  // concentration
    let conc_bu = vars.add_sclbnd_cnst(bnd_up, 0.0, "")?;  // concentration
    let mflx_mb = vars.add_sclbnd_cnst(bnd_mb, 0.0, "")?;  // molar flux (positive if outward; negative if inward)
    let mflx_mt = vars.add_sclbnd_cnst(bnd_mt, 0.0, "")?;  // molar flux

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // add_mcnt_itf - register interface with concentration and flux continuity
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(1)?;
    phys.add_mass_dom(0, dom_lo, conc_lo, HashMap::from([(0, diff_lo)]), mprd_lo, msrc_lo)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(0, dom_m, conc_m, HashMap::from([(0, diff_m)]), mprd_m, msrc_m)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(0, dom_up, conc_up, HashMap::from([(0, diff_up)]), mprd_up, msrc_up)?;  // arguments: component, domain, c, D, s, R
    phys.add_conc_bnd(0, bnd_lo, conc_bl)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(0, bnd_up, conc_bu)?;  // arguments: component, boundary, c
    phys.add_mflx_bnd(0, bnd_mb, mflx_mb)?;  // arguments: component, boundary, N
    phys.add_mflx_bnd(0, bnd_mt, mflx_mt)?;  // arguments: component, boundary, N
    phys.add_mcnt_itf(0, itf_lo)?;  // arguments: component, interface
    phys.add_mcnt_itf(0, itf_up)?;  // arguments: component, interface

    // physics solver
    let mut lin = LinearGmres::new_with_options(1, 1e-10, 50, 200)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    Ok(())
}
