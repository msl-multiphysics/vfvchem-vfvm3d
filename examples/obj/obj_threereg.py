from vfvchem_geom3d import Geom3D

# Cube with lower-front-left and upper-back-right corner cuboids
#
# Boundaries and Interfaces
# bnd_lo - lower-front-left outer (left + front + bottom)
# bnd_mb - middle outer (bottom + front + right)
# bnd_up - upper-back-right outer (right + back + top)
# bnd_mt - middle outer (top + back + left)
# itf_lo - lower-front-left interface
# itf_up - upper-back-right interface
#
# Domains
# dom_m - middle
# dom_lo - lower-front-left
# dom_up - upper-back-right
#
# Group names must begin with dom_, bnd_, or itf_
# The ids shown before grouping are only used to select surfaces and domains
#

# define geometry
# box2 and box3 overlap with box1; box1 will be automatically cut into three regions
geom = Geom3D("obj_threereg")
box1 = geom.add_box(0.00, 0.00, 0.00, 1.00, 1.00, 1.00)  # arguments: x_min, y_min, z_min, x_max, y_max, z_max
box2 = geom.add_box(0.00, 0.00, 0.00, 0.75, 0.25, 0.25)  # arguments: x_min, y_min, z_min, x_max, y_max, z_max
box3 = geom.add_box(0.25, 0.75, 0.75, 1.00, 1.00, 1.00)  # arguments: x_min, y_min, z_min, x_max, y_max, z_max
geom.finalize_geometry()  # call at the end of geometry definition

# show geometry with domain, boundary, and interface ids
# this will show the automatically generated ids
geom.show_geometry()  # call to show geometry with ids

# name and/or join ids
# arguments: group name, list of ids to join
geom.group_dom("dom_m", [2])  # middle
geom.group_dom("dom_lo", [0])  # lower-front-left
geom.group_dom("dom_up", [1])  # upper-back-right
geom.group_bnd("bnd_lo", [6, 7, 8])  # lower-front-left outer (left + front + bottom)
geom.group_bnd("bnd_mb", [4, 1, 5])  # middle outer (bottom + front + right)
geom.group_bnd("bnd_up", [9, 10, 11])  # upper-back-right outer (right + back + top)
geom.group_bnd("bnd_mt", [2, 3, 0])  # middle outer (top + back + left)
geom.group_itf("itf_lo", [12, 13, 14])  # lower-front-left interface
geom.group_itf("itf_up", [15, 16, 17])  # upper-back-right interface
geom.update_group()  # call to check the groups

# show group names
geom.show_geometry()

# export geometry
# arguments: file name, h_min (optional; boundary node spacing used in the solver)
# h_min is only used to check that faces are large enough for the solver; it does not modify the geometry
geom.export_geometry("obj_threereg.obj", 0.02)
