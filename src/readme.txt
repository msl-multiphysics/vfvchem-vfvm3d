__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
3D Vertex-Centered Finite Volume Method Module - Source Code

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.

Nodes are placed on boundary and interface faces, and points are placed inside domains via Poisson-disk sampling.
Each domain is tetrahedralized with a Delaunay tetrahedralization, and a Voronoi control volume is built around every vertex.
Currently, only backward Euler time stepping is implemented, with constant or adaptive time steps.
Advective transport is discretized via the power law scheme.
Pressure-velocity coupling in the Navier-Stokes equation is handled via Rhie-Chow interpolation.

The source code is organized into the following folders:
- base      : core data structures (e.g., geometry, vertices, scalars, vectors, writers)
- linear    : linear matrix solvers (e.g., LU, GMRES)
- operator  : applies discretized conservation equations into global matrix equation
- steady    : steady-state solvers
- transient : transient solvers

The operator, steady, and transient folders are split into two subfolders:
- core   : heat, mass, and momentum transfer in single-phase materials
- extend : additional physics built on the core (currently empty)
The base classes (operator_base, steady_base, transient_base, transient_adaptive)
and the preludes stay at the top of each folder. The preludes re-export both
core and extend types, so `use vfvchem_vfvm3d::*` covers both.

Base
----
- delaunay : Delaunay tetrahedralization of a point set (Bowyer-Watson with exact predicates)
- error    : defines errors returned by the library
- geom_bnd : stores 2D boundary geometry
- geom_dom : stores 3D domain geometry
- geom_itf : stores 2D interface geometry; interfaces join two domains
- geom_tet : stores the tetrahedralization and control volume data of a domain
- nds      : stores nodes on boundaries and interfaces; placed via Poisson-disk sampling on each face
- pts      : stores points inside domains; placed via Poisson-disk sampling
- read_obj : imports boundary, interface, and domain surfaces from an obj file
- scaling  : per-field reference scales for the Picard residual and the adaptive error norm
- scl_bnd  : stores a scalar field defined on a boundary
- scl_dom  : stores a scalar field defined on a domain
- scl_flw  : stores the Rhie-Chow face mass flux shared by the momentum and continuity equations; filled by the solvers
- scl_itf  : stores a scalar field defined on an interface
- vars     : owns all geometry, scalar, and vector data
- vec_bnd  : stores a 3D vector field defined on a boundary
- vec_dom  : stores a 3D vector field defined on a domain
- vec_itf  : stores a 3D vector field defined on an interface
- write_csv: writes scalar and vector values to CSV files
- write_vtu: writes scalar and vector values to VTU files

Vertices are either points or nodes. Within a domain, vertex ids are ordered with
points first, then nodes, i.e., vid = pid for points and vid = num_pts + nid for nodes.

Mesh generation
- read_obj groups the polygons of an obj file by name; groups starting with dom are domains,
  and groups starting with bnd or itf are boundaries or interfaces. Coplanar neighboring polygons of a group are
  merged into one face, so a face of a box is one face even if it is split into triangles.
- nds places nodes on every face with Poisson-disk sampling. Nodes are inset from the face edges
  by 0.25 h; edges and corners themselves carry no nodes, so they are slightly chamfered.
- pts places points inside each domain with Poisson-disk sampling; spacing grows from h_min at the
  surface to h_max at ref_len. Points are kept outside the circumspheres of the surface triangles
  (protecting balls) so that the surface triangles appear in the Delaunay tetrahedralization.
- geom_tet tetrahedralizes the points and nodes of each domain, keeps the tetrahedra whose centroid
  is inside the domain, and builds the Voronoi dual: control volume faces join edge midpoints,
  face circumcenters, and tet circumcenters. On the surface, the control volume of a vertex has one
  sub-face per surface triangle around it; each sub-face carries interpolation weights for the vertex
  and its two neighbors on the surface triangle.
- new_bounds builds an axis-aligned box with the six boundaries bnd_l (x_min), bnd_r (x_max),
  bnd_f (y_min), bnd_k (y_max), bnd_b (z_min), and bnd_t (z_max).

Field types are registered through Variables helpers:
- *_cnst    : constant value
- *_func    : function of time, position (x, y, z), and unknown scalar values;
              interface functions receive time and unknown values from both domains
- *_unkn    : unknown field (solved for)

Vector values are (x, y, z) tuples. Vector gradients are 9-tuples in row-major order,
(u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z). In the global matrix, a vector unknown takes
three blocks of num_tot rows, one per component.

Flow data is registered through add_sclflw and is used by the momentum and
continuity operators so that both use the same face mass flux.
SteadyBase and TransientBase call update_func, update_grad, and update_flow
before each matrix assembly so that functions, gradients, and face mass fluxes
stay consistent with the unknowns. update_flow is a required method of each
solver because the Rhie-Chow coefficient depends on the momentum equation being
solved: flow solvers call flow_helper::update_flow, and solvers without flow leave it empty.

Operator
--------
Operator names use the form op + field type + location + term. For example,
opscl_dom_diff is a scalar domain diffusion operator. Files ending in unity
are variants whose weighting coefficient is fixed at one, and files ending in
timeprev use the weighting coefficient of the previous time step.
Domain operators apply to the whole domain, boundary operators apply to the
control volumes of boundary nodes, and interface operators apply to the
combined control volumes of interface nodes.

Helpers hold code shared by several solvers. They are not re-exported by the
prelude and are called by path (e.g., flow_helper::update_flow).

- operator_base       : defines the common operator trait and global matrix assembly helpers
- prelude             : re-exports operator types for convenient internal imports

Core operators (operator/core)
- flow_helper         : computes the Rhie-Chow coefficient and face mass flux of a free fluid
- opscl_bnd_adv       : adds weighted advection to the control volumes of boundary nodes
- opscl_bnd_advunity  : adds advection with unit weight to the control volumes of boundary nodes
- opscl_bnd_diff      : adds diffusion to the control volumes of boundary nodes
- opscl_bnd_dir       : enforces a scalar Dirichlet boundary condition
- opscl_bnd_div       : adds velocity divergence to the control volumes of boundary nodes
- opscl_bnd_neu       : applies a prescribed scalar outward flux
- opscl_bnd_prod      : adds a production term (s * c) to the control volumes of boundary nodes
- opscl_bnd_pstab     : adds the Rhie-Chow correction to the control volumes of boundary nodes
- opscl_bnd_src       : adds a source to the control volumes of boundary nodes
- opscl_bnd_time      : adds a weighted time derivative to the control volumes of boundary nodes
- opscl_bnd_timeprev  : adds a time derivative with the previous weight to the control volumes of boundary nodes
- opscl_bnd_timeunity : adds a time derivative with unit weight to the control volumes of boundary nodes
- opscl_dom_adv       : adds weighted advection to a scalar transport equation
- opscl_dom_advunity  : adds advection with unit weight to a scalar transport equation
- opscl_dom_diff      : adds diffusion to a scalar transport equation
- opscl_dom_div       : adds velocity divergence to the continuity equation
- opscl_dom_pref      : sets the pressure at one reference node
- opscl_dom_prod      : adds a production term (s * c) to a scalar transport equation
- opscl_dom_pstab     : adds the Rhie-Chow correction to the continuity equation
- opscl_dom_src       : adds a source to a scalar transport equation
- opscl_dom_time      : adds a weighted time derivative to a scalar transport equation
- opscl_dom_timeprev  : adds a time derivative with the previous weight to a scalar transport equation
- opscl_dom_timeunity : adds a time derivative with unit weight to a scalar transport equation
- opscl_itf_adv       : adds weighted advection to the combined control volumes of interface nodes
- opscl_itf_advunity  : adds advection with unit weight to the combined control volumes of interface nodes
- opscl_itf_cont      : enforces scalar continuity between two domains
- opscl_itf_diff      : adds diffusion to the combined control volumes of interface nodes
- opscl_itf_div       : adds velocity divergence to the combined control volumes of interface nodes
- opscl_itf_prod      : adds a production term (s * c) to the combined control volumes of interface nodes
- opscl_itf_pstab     : adds the Rhie-Chow correction to the combined control volumes of interface nodes
- opscl_itf_src       : adds a source to the combined control volumes of interface nodes
- opscl_itf_time      : adds a weighted time derivative to the combined control volumes of interface nodes
- opscl_itf_timeprev  : adds a time derivative with the previous weight to the combined control volumes of interface nodes
- opscl_itf_timeunity : adds a time derivative with unit weight to the combined control volumes of interface nodes
- opvec_bnd_adv       : adds momentum advection to the control volumes of boundary nodes
- opvec_bnd_ddil      : adds the dilatatory viscous term to the control volumes of boundary nodes
- opvec_bnd_diff      : adds the viscous diffusion term to the control volumes of boundary nodes
- opvec_bnd_dir       : enforces a vector Dirichlet boundary condition
- opvec_bnd_dtrn      : adds the transpose viscous term to the control volumes of boundary nodes
- opvec_bnd_grad      : adds the pressure gradient to the control volumes of boundary nodes
- opvec_bnd_src       : adds a body force to the control volumes of boundary nodes
- opvec_bnd_time      : adds a momentum time derivative to the control volumes of boundary nodes
- opvec_bnd_trac      : sets the traction on a traction outflow boundary
- opvec_dom_adv       : adds advection to the momentum equation
- opvec_dom_ddil      : adds the dilatatory viscous term to the momentum equation
- opvec_dom_diff      : adds the viscous diffusion term to the momentum equation
- opvec_dom_dtrn      : adds the transpose viscous term to the momentum equation
- opvec_dom_grad      : adds the pressure gradient to the momentum equation
- opvec_dom_src       : adds a body force to the momentum equation
- opvec_dom_time      : adds a time derivative to the momentum equation
- opvec_itf_adv       : adds momentum advection to the combined control volumes of interface nodes
- opvec_itf_cont      : enforces velocity continuity between two domains
- opvec_itf_ddil      : adds the dilatatory viscous term to the combined control volumes of interface nodes
- opvec_itf_diff      : adds the viscous diffusion term to the combined control volumes of interface nodes
- opvec_itf_dtrn      : adds the transpose viscous term to the combined control volumes of interface nodes
- opvec_itf_grad      : adds the pressure gradient to the combined control volumes of interface nodes
- opvec_itf_src       : adds a body force to the combined control volumes of interface nodes
- opvec_itf_time      : adds a momentum time derivative to the combined control volumes of interface nodes

Extended operators (operator/extend)
None yet; add operators here and re-export them in operator/prelude.rs.

Linear
------
- linear_base  : defines the common linear solver trait
- linear_gmres : solves sparse systems with the Jacobi-preconditioned restarted GMRES method
- linear_lu    : solves sparse systems with an LU factorization

Steady
------
- steady_base         : implements the common steady nonlinear solution loop

Core solvers (steady/core)
- steady_flow         : solves steady momentum transfer
- steady_heat         : solves steady heat transfer
- steady_heatflow     : solves coupled steady heat and momentum transfer
- steady_heatmassflow : solves coupled steady heat, mass, and momentum transfer
- steady_mass         : solves steady mass transfer
- steady_massflow     : solves coupled steady mass and momentum transfer

Extended solvers (steady/extend)
None yet.

Transient
---------
- transient_base         : implements the common time-stepping and nonlinear solution loops
- transient_adaptive     : implements adaptive time stepping for every transient solver

Core solvers (transient/core)
- transient_flow         : solves transient momentum transfer
- transient_heat         : solves transient heat transfer
- transient_heatflow     : solves coupled transient heat and momentum transfer
- transient_heatmassflow : solves coupled transient heat, mass, and momentum transfer
- transient_mass         : solves transient mass transfer
- transient_massflow     : solves coupled transient mass and momentum transfer

Extended solvers (transient/extend)
None yet.
