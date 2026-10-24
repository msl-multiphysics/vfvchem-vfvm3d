use vfvchem_vfvm3d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation with lid-driven cavity flow in three regions joined by interfaces.
/// Run with: `cargo run --release --example flow_threereg`
///
/// Geometry:
/// Cube with a lower-front-left and an upper-back-right corner cuboid; see examples/obj/obj_threereg.py.
/// The top face (z = 1.0 m) is the lid; it is shared by bnd_mt and bnd_up.
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
/// Properties:
/// dom_lo, dom_m, dom_up - density (rho = 1000.0 kg m-3)
/// dom_lo, dom_m, dom_up - viscosity (mu = 0.001 Pa s)
/// dom_lo, dom_m, dom_up - body force (f = <0.0, 0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_lo - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_mb - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_mt - lid velocity on the top face; no-slip on the back and left faces
/// bnd_up - lid velocity on the top face; no-slip on the right and back faces
/// v = <4.0e-4, 0.0, 0.0> if z = 1.0; else <0.0, 0.0, 0.0>
/// Zero reference pressure at node 0 of dom_m
///
/// Interface conditions:
/// itf_lo - velocity and pressure continuity
/// itf_up - velocity and pressure continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_threereg").unwrap();

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
    // vector values are given as (x, y, z) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_lo_0000.vtu)
    let vel_lo = vars.add_vecdom_unkn(dom_lo, (0.0, 0.0, 0.0), "examples/output_flow_threereg/vel_lo.vtu")?;
    let vel_m = vars.add_vecdom_unkn(dom_m, (0.0, 0.0, 0.0), "examples/output_flow_threereg/vel_m.vtu")?;
    let vel_up = vars.add_vecdom_unkn(dom_up, (0.0, 0.0, 0.0), "examples/output_flow_threereg/vel_up.vtu")?;
    let pres_lo = vars.add_scldom_unkn(dom_lo, 0.0, "examples/output_flow_threereg/pres_lo.vtu")?;
    let pres_m = vars.add_scldom_unkn(dom_m, 0.0, "examples/output_flow_threereg/pres_m.vtu")?;
    let pres_up = vars.add_scldom_unkn(dom_up, 0.0, "examples/output_flow_threereg/pres_up.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    let den_lo = vars.add_scldom_cnst(dom_lo, 1000.0, "")?;  // density
    let den_m = vars.add_scldom_cnst(dom_m, 1000.0, "")?;  // density
    let den_up = vars.add_scldom_cnst(dom_up, 1000.0, "")?;  // density
    let visc_lo = vars.add_scldom_cnst(dom_lo, 0.001, "")?;  // viscosity
    let visc_m = vars.add_scldom_cnst(dom_m, 0.001, "")?;  // viscosity
    let visc_up = vars.add_scldom_cnst(dom_up, 0.001, "")?;  // viscosity
    let fce_lo = vars.add_vecdom_cnst(dom_lo, (0.0, 0.0, 0.0), "")?;  // body force (den * g; not acceleration g)
    let fce_m = vars.add_vecdom_cnst(dom_m, (0.0, 0.0, 0.0), "")?;  // body force (den * g; not acceleration g)
    let fce_up = vars.add_vecdom_cnst(dom_up, (0.0, 0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y, z), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // the lid moves along x; the other faces of the same boundaries are no-slip
    let lid_func = |_t: f64, (_x, _y, z): (f64, f64, f64), _unkn: &[f64]| if z > 1.0 - 1e-9 { (4e-4, 0.0, 0.0) } else { (0.0, 0.0, 0.0) };
    let vel_mt = vars.add_vecbnd_func(bnd_mt, Box::new(lid_func), vec![], "")?;  // velocity
    let vel_bu = vars.add_vecbnd_func(bnd_up, Box::new(lid_func), vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let vel_bl = vars.add_vecbnd_cnst(bnd_lo, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_mb = vars.add_vecbnd_cnst(bnd_mb, (0.0, 0.0, 0.0), "")?;  // velocity

    // steady-state momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_fcnt_itf - register interface with velocity and pressure continuity
    // set_pres_ref - set reference pressure at a node; needed if no boundary sets the pressure
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_lo, vel_lo, pres_lo, den_lo, visc_lo, fce_lo)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_flow_dom(dom_m, vel_m, pres_m, den_m, visc_m, fce_m)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_flow_dom(dom_up, vel_up, pres_up, den_up, visc_up, fce_up)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_lo, vel_bl)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_mb, vel_mb)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_mt, vel_mt)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_up, vel_bu)?;  // arguments: boundary, v
    phys.add_fcnt_itf(itf_lo)?;  // arguments: interface
    phys.add_fcnt_itf(itf_up)?;  // arguments: interface
    phys.set_pres_ref(dom_m, 0, 0.0)?;  // arguments: domain, node_index, p

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // new_with_options - arguments: num_thread, tol, restart, max_outer; a larger restart is needed across the interfaces
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearGmres::new_with_options(1, 1e-6, 100, 400)?;
    phys.solve(&mut vars, &mut lin, 100, 1e-3, 0.6)?;

    Ok(())
}
