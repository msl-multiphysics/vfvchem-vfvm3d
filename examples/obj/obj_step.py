from vfvchem_geom3d import Geom3D

# Channel with a backward-facing step, extruded in y
#
# Boundaries
# bnd_i - inlet
# bnd_o - outlet
# bnd_w - walls (bottom, step, top, front, and back)
#
# Domains
# dom_f - fluid
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select surfaces and domains
#

# define geometry
geom = Geom3D("obj_step")
reg1 = geom.add_extruded_polygon(  # channel; arguments: (x, z) coords, axis, c_min, c_max
    [
        (0.0, 0.2),
        (0.4, 0.2),
        (0.4, 0.0),
        (2.4, 0.0),
        (2.4, 0.4),
        (0.0, 0.4),
    ],
    1,  # extrude along y
    0.0,
    0.4,
)
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_f", [0])  # fluid
geom.group_bnd("bnd_i", [6])  # inlet
geom.group_bnd("bnd_o", [4])  # outlet
geom.group_bnd("bnd_w", [3, 1, 2, 5, 0, 7])  # walls (bottom, step, top, front, and back)
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in the solver)
# h_min is only used to check that faces are large enough for the solver; it does not modify the geometry
geom.export_geometry("obj_step.obj", 0.02)
