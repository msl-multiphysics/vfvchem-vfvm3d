__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
3D Vertex-Centered Finite Volume Method Module - Examples

How to Use Examples
-------------------
First-time users are encouraged to run quickstart.rs

Change the directory to vfvchem_vfvm3d then run quickstart.rs:
`cargo run --release --example quickstart`

The other examples are registered in `Cargo.toml` and can be run in place:
`cargo run --release --example name_of_file`

Geometries used by the examples are stored as `.obj` files in the `examples/obj` folder.
The Python scripts in the same folder generate these files with VFVChem geom3d.
The box examples use `Variables::new_bounds` and need no `.obj` file.

3D problems have many more vertices than 2D problems at the same spacing, so the meshes of the
examples are coarse. Lower h_min and h_max for more accurate results at a higher cost.

List of Examples
----------------
Below is a list of examples and the concepts they introduce:

Heat Transfer
- quickstart.rs       : introduction to VFVChem code
- heat_transient.rs   : transient heat transfer simulation
- heat_threereg.rs    : importing `.obj` geometries; multiple domains joined by interfaces
- heat_analytic.rs    : steady-state conduction in a uniformly heated slab

Mass Transfer
- mass_cube.rs     : steady-state multicomponent diffusion with cross-diffusion
- mass_analytic.rs : steady-state diffusion with a first-order reaction using a production coefficient

Momentum Transfer
- flow_lid.rs      : lid-driven cavity flow
- flow_step.rs     : transient flow past a backward-facing step
- flow_threereg.rs : lid-driven cavity flow through three domains joined by interfaces
- flow_analytic.rs : Poiseuille flow between parallel plates

Features
- feat_prod.rs  : production coefficient versus mass source for a fast first-order reaction
- feat_itf.rs   : heat and mass transfer through regions with different properties joined by interfaces
