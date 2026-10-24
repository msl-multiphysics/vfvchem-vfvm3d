use vfvchem_vfvm3d::*;
use std::fs::create_dir_all;

/// Transient Navier-Stokes equation with backward-facing step flow.
/// Run with: `cargo run --release --example flow_step`
///
/// Geometry:
/// Channel (2.4 m x 0.4 m x 0.4 m) with a step (0.4 m x 0.4 m x 0.2 m) at the bottom left corner.
/// Inlet above the step, outlet at the right end, and walls elsewhere.
/// Side view in the x-z plane; the channel is extruded in y from 0.0 m to 0.4 m.
///
///                        bnd_w
///       +-------------------------------------+
/// bnd_i |                                     |
///       +----+             dom_f              | bnd_o
///            |                                |
///            +--------------------------------+
///                        bnd_w
///
/// Domains:
/// dom_f - channel
///
/// Boundaries:
/// bnd_i - inlet (x = 0.0 m; 0.0 m < y < 0.4 m; 0.2 m < z < 0.4 m)
/// bnd_o - outlet (x = 2.4 m; 0.0 m < y < 0.4 m; 0.0 m < z < 0.4 m)
/// bnd_w - walls (bottom, step, top, front at y = 0.0 m, and back at y = 0.4 m)
///
/// Properties:
/// dom_f - density (rho = 1.2 kg m-3)
/// dom_f - viscosity (mu = 1.8e-5 Pa s)
/// dom_f - body force (f = <0.0, 0.0, 0.0> N m-3)
///
/// Initial conditions:
/// dom_f - velocity (v = <0.015, 0.0, 0.0> m s-1)
/// dom_f - pressure (p = 0.0 Pa)
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_i - parabolic velocity in y and z (v = <u, 0.0, 0.0> m s-1)
/// u = 16 * 0.0225 * y * (0.4 - y) * (z - 0.2) * (0.4 - z) / (0.4^2 * 0.2^2)
/// bnd_o - pressure (p = 0.0 Pa)
/// bnd_w - no-slip (v = <0.0, 0.0, 0.0> m s-1)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_step").unwrap();

    // nodes, points, and variables
    // new_obj - import boundary, interface, and domain surfaces from an obj file
    // arguments: input_file, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    let mut vars = Variables::new_obj("examples/obj/obj_step.obj", 0.04, 0.08, 0.2, 100, 42)?;

    // geometry
    // domains, boundaries, and interfaces are selected by their group names in the obj file
    // add_bnd - arguments: domain, group_name
    let dom_f = vars.add_dom("dom_f")?;  // channel
    let bnd_i = vars.add_bnd(dom_f, "bnd_i")?;  // inlet
    let bnd_o = vars.add_bnd(dom_f, "bnd_o")?;  // outlet
    let bnd_w = vars.add_bnd(dom_f, "bnd_w")?;  // walls

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // vector values are given as (x, y, z) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_0000.vtu)
    let vel = vars.add_vecdom_unkn(dom_f, (0.015, 0.0, 0.0), "examples/output_flow_step/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom_f, 0.0, "examples/output_flow_step/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    let den = vars.add_scldom_cnst(dom_f, 1.2, "")?;  // density
    let visc = vars.add_scldom_cnst(dom_f, 1.8e-5, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom_f, (0.0, 0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position (x, y, z), and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    // parabolic inlet on y in [0.0, 0.4] and z in [0.2, 0.4], maximum 0.0225
    let vel_func: VecBndFunc = Box::new(|_t, (_x, y, z), _unkn| {
        let u = 16.0 * 0.0225 * y * (0.4 - y) * (z - 0.2) * (0.4 - z) / (0.4 * 0.4 * 0.2 * 0.2);
        (u, 0.0, 0.0)
    });
    let vel_i = vars.add_vecbnd_func(bnd_i, vel_func, vec![], "")?;  // velocity

    // constant boundary values
    // arguments: boundary, value, output_file
    let pres_o = vars.add_sclbnd_cnst(bnd_o, 0.0, "")?;  // pressure
    let vel_w = vars.add_vecbnd_cnst(bnd_w, (0.0, 0.0, 0.0), "")?;  // velocity

    // transient momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // add_pres_bnd - register boundary with pressure
    let mut phys = TransientFlow::new()?;
    phys.add_flow_dom(dom_f, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_i, vel_i)?;  // arguments: boundary, v
    phys.add_pres_bnd(bnd_o, pres_o)?;  // arguments: boundary, p
    phys.add_vel_bnd(bnd_w, vel_w)?;  // arguments: boundary, v

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    // num_ts_write - writes an output file every num_ts_write time steps; the initial condition is output 0
    // max_iter, tol - maximum number of iterations and tolerance in each time step
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearGmres::new(1)?;
    phys.solve(&mut vars, &mut lin, 20.0, 50, 1, 30, 1e-3, 0.5)?;

    Ok(())
}
