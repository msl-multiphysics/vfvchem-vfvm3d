use vfvchem_vfvm3d::*;
use std::fs::create_dir_all;

/// Steady-state Navier-Stokes equation with lid-driven cavity flow.
/// Run with: `cargo run --release --example flow_lid`
///
/// Geometry:
/// Cube domain with smaller point spacing near the walls.
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
/// dom - main domain (1.0 m x 1.0 m x 1.0 m)
///
/// Boundaries:
/// bnd_l - left boundary (x = 0.0 m)
/// bnd_r - right boundary (x = 1.0 m)
/// bnd_f - front boundary (y = 0.0 m)
/// bnd_k - back boundary (y = 1.0 m)
/// bnd_b - bottom boundary (z = 0.0 m)
/// bnd_t - top boundary (z = 1.0 m)
///
/// Properties:
/// dom - density (rho = 1000.0 kg m-3)
/// dom - viscosity (mu = 0.001 Pa s)
/// dom - body force (f = <0.0, 0.0, 0.0> N m-3)
///
/// Boundary conditions:
/// bnd_l - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_r - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_f - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_k - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_b - no-slip (v = <0.0, 0.0, 0.0> m s-1)
/// bnd_t - velocity (v = <4.0e-4, 0.0, 0.0> m s-1)
/// Zero reference pressure at node 0 of dom
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_lid").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points inside a bounding box
    // arguments: x_min, y_min, z_min, x_max, y_max, z_max, h_min, h_max, ref_len, max_try, seed
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    // max_try - attempts to place a point before giving up; seed - random seed for point sampling
    // points are refined near the walls, where the velocity gradients are largest
    let mut vars = Variables::new_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.05, 0.1, 0.2, 50, 42)?;

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
    // vector values are given as (x, y, z) tuples
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_0000.vtu)
    let vel = vars.add_vecdom_unkn(dom, (0.0, 0.0, 0.0), "examples/output_flow_lid/vel.vtu")?;
    let pres = vars.add_scldom_unkn(dom, 0.0, "examples/output_flow_lid/pres.vtu")?;

    // constant domain values
    // arguments: domain, value, output_file
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;  // density
    let visc = vars.add_scldom_cnst(dom, 0.001, "")?;  // viscosity
    let fce = vars.add_vecdom_cnst(dom, (0.0, 0.0, 0.0), "")?;  // body force (den * g; not acceleration g)

    // constant boundary values
    // arguments: boundary, value, output_file
    let vel_l = vars.add_vecbnd_cnst(bnd_l, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_r = vars.add_vecbnd_cnst(bnd_r, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_f = vars.add_vecbnd_cnst(bnd_f, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_k = vars.add_vecbnd_cnst(bnd_k, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_b = vars.add_vecbnd_cnst(bnd_b, (0.0, 0.0, 0.0), "")?;  // velocity
    let vel_t = vars.add_vecbnd_cnst(bnd_t, (4e-4, 0.0, 0.0), "")?;  // velocity

    // steady-state momentum transfer solver
    // add_flow_dom - register domain with momentum transfer
    // add_vel_bnd - register boundary with velocity
    // set_pres_ref - set reference pressure at a node; needed if no boundary sets the pressure
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, pres, den, visc, fce)?;  // arguments: domain, v, p, rho, mu, f
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_r, vel_r)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_f, vel_f)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_k, vel_k)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_b, vel_b)?;  // arguments: boundary, v
    phys.add_vel_bnd(bnd_t, vel_t)?;  // arguments: boundary, v
    phys.set_pres_ref(dom, 0, 0.0)?;  // arguments: domain, node_index, p

    // physics solver
    // LinearGmres::new - Jacobi-preconditioned restarted GMRES linear solver; arguments: num_thread
    // GMRES is much faster than LU (LinearLu::new) in 3D because the LU factors fill in
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearGmres::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
