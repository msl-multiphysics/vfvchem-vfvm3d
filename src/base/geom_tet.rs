use crate::base::prelude::*;
use std::collections::{HashMap, HashSet};

#[derive(Default)]
pub struct Tetrahedron {
    // ids
    pub tet_id: usize,
    pub dom_id: usize,  // domain this tet refers to

    // tetrahedralization data
    // vid is the index of combined point and node vectors
    // i.e., vid = if vid < dom.num_pts() {pid} else {dom.num_pts() + nid}
    pub tet_vid: Vec<[usize; 4]>,  // [vid, vid, vid, vid] positively oriented; the fourth vertex is below the counterclockwise first three

    // control volume data
    pub cvol_vol: Vec<f64>,  // vid (vertex id) -> CV (control volume) volume
    pub cvol_skew: Vec<(f64, f64, f64)>,  // vid -> CV centroid minus vertex position

    // face data
    pub face_area: Vec<Vec<f64>>,  // vid, loc (local neighbor index) -> signed face areas; negative if circumcenters cross the edge
    pub face_norm: Vec<Vec<(f64, f64, f64)>>,  // vid, loc -> outward normal on face
    pub face_skew: Vec<Vec<(f64, f64, f64)>>,  // vid, loc -> face center minus edge midpoint

    // neighbor vertex data
    pub ngh_vid: Vec<Vec<usize>>,  // vid -> neighbor vid
    pub ngh_dist: Vec<Vec<f64>>,  // vid -> edge length to neighbor

    // boundary vertex data
    // boundary vertices own the part of each surface triangle closer to them than to the other two vertices (circumcentric)
    // values on a sub-face are interpolated linearly from the triangle vertices to the sub-face centroid
    pub bnd_sum: Vec<f64>,  // vid -> CV area on boundary; zero for points
    pub bnd_vid: Vec<Vec<(usize, usize)>>,  // vid, blc (boundary loc) -> other two vids of the surface triangle; empty for points
    pub bnd_area: Vec<Vec<f64>>,  // vid, blc -> signed sub-face area
    pub bnd_norm: Vec<Vec<(f64, f64, f64)>>,  // vid, blc -> outward normal on sub-face
    pub bnd_wgt: Vec<Vec<(f64, f64, f64)>>,  // vid, blc -> interpolation weights (w_v, w_j, w_k) of the sub-face centroid

    // gradient weight data
    // delta_u = [u1 - u0, u2 - u0, ...]; u0 and uj are point and neighbor values
    // du/dx = wgt_x . delta_u; du/dy = wgt_y . delta_u; du/dz = wgt_z . delta_u
    pub grad_wgt: Vec<Vec<(f64, f64, f64)>>,  // vid, loc -> least-squares weights (wgt_x, wgt_y, wgt_z)

    // pressure gradient data
    // int_vol grad(p) dV = sum_ngh p_f * n * dA + sum_bnd p_b * n * dA = sum_m (c_x, c_y, c_z) * p_m
    // shared by the momentum and continuity operators so that both use the same pressure gradient
    pub grad_coeff: Vec<Vec<(usize, f64, f64, f64)>>,  // vid -> (vid_m, c_x, c_y, c_z) divergence theorem coefficients
}

impl Tetrahedron {
    pub fn new(tet_id: usize, dom: &Domain) -> Result<Tetrahedron, VFVChemError> {
        // default struct
        let mut tet = Tetrahedron::default();
        tet.tet_id = tet_id;
        tet.dom_id = dom.dom_id;

        Ok(tet)
    }

    pub fn update_tet(&mut self, nds: &Nodes, pts: &Points, dom: &Domain) -> Result<(), VFVChemError> {
        // clear tet data
        self.tet_vid.clear();
        self.cvol_vol = vec![0.0; dom.num_tot()];
        self.cvol_skew = vec![(0.0, 0.0, 0.0); dom.num_tot()];
        self.face_area = vec![Vec::new(); dom.num_tot()];
        self.face_norm = vec![Vec::new(); dom.num_tot()];
        self.face_skew = vec![Vec::new(); dom.num_tot()];
        self.ngh_vid = vec![Vec::new(); dom.num_tot()];
        self.ngh_dist = vec![Vec::new(); dom.num_tot()];
        self.bnd_sum = vec![0.0; dom.num_tot()];
        self.bnd_vid = vec![Vec::new(); dom.num_tot()];
        self.bnd_area = vec![Vec::new(); dom.num_tot()];
        self.bnd_norm = vec![Vec::new(); dom.num_tot()];
        self.bnd_wgt = vec![Vec::new(); dom.num_tot()];
        self.grad_wgt = vec![Vec::new(); dom.num_tot()];
        self.grad_coeff = vec![Vec::new(); dom.num_tot()];

        // collect points in domain
        let mut pts_pos = Vec::with_capacity(dom.num_tot());  // vid -> (x, y, z)
        for &gpt in &dom.pts_gpt {
            pts_pos.push(pts.pts_pos[gpt]);
        }
        for &gnd in &dom.nds_gnd {
            pts_pos.push(nds.nds_pos[gnd]);
        }

        // step 1: Delaunay tetrahedralization

        // tetrahedralize the convex hull of all vertices
        let tet_all = delaunay_tet(&pts_pos).ok_or(VFVChemError::TetraDegenerateInput { dom_id: self.dom_id })?;

        // keep tetrahedra whose centroid is inside the closed surface of the domain
        // the surface is not required to be made of tetrahedron faces; the domain boundary is taken from the kept tetrahedra
        let surf_tri: Vec<[usize; 3]> = nds.grp3d_tri[&dom.grp3d_id].iter().map(|t| t.map(|gnd| dom.nds_nid[&gnd])).collect();
        let surf = SurfaceGrid::new(&pts_pos[dom.num_pts()..], &surf_tri, pts.h_min.max(0.25 * pts.h_max));
        for t in tet_all {
            let cent = scale3(0.25, add3(add3(pts_pos[t[0]], pts_pos[t[1]]), add3(pts_pos[t[2]], pts_pos[t[3]])));
            if surf.is_inside(cent) { self.tet_vid.push(t); }
        }

        // boundary faces are faces of exactly one kept tetrahedron
        // faces opposite each vertex of a positive tetrahedron, counterclockwise seen from outside
        let mut face_cnt: HashMap<[usize; 3], usize> = HashMap::new();  // sorted face vids -> number of kept tetrahedra
        for t in &self.tet_vid {
            for k in 0..4 {
                let mut key = [t[(k + 1) % 4], t[(k + 2) % 4], t[(k + 3) % 4]];
                key.sort_unstable();
                *face_cnt.entry(key).or_insert(0) += 1;
            }
        }
        let mut surf_vid = Vec::new();  // [vid, vid, vid] counterclockwise seen from outside
        for t in &self.tet_vid {
            for (k, [a, b, c]) in [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]].into_iter().enumerate() {
                let mut key = [t[(k + 1) % 4], t[(k + 2) % 4], t[(k + 3) % 4]];
                key.sort_unstable();
                if face_cnt[&key] == 1 { surf_vid.push([t[a], t[b], t[c]]); }
            }
        }

        // boundary faces must only join nodes
        let num_bnd = surf_vid.iter().flatten().filter(|&&vid| vid < dom.num_pts()).collect::<HashSet<_>>().len();
        if num_bnd > 0 { return Err(VFVChemError::TetraPointsOnBoundary { dom_id: self.dom_id, num_bnd }) }

        // step 2: Voronoi faces

        // iterate through edges of each tetrahedron
        // the face piece of edge (i, j) on the side of face (i, j, k) is the right triangle (m, c_ijk, c_T)
        // m is the edge midpoint, c_ijk the face circumcenter, and c_T the tetrahedron circumcenter
        // its signed area is 0.5 * s1 * s2 with s1 = (c_ijk - m) . u toward k and s2 = (c_T - c_ijk) . n toward l
        // the signed pyramids from vertices i and j over each piece tile the tetrahedron
        let mut edge_face = HashMap::new();  // (vid_min, vid_max) -> (signed face area, sum of area * (piece center - m))
        let mut cent_sum = vec![(0.0, 0.0, 0.0); dom.num_tot()];  // vid -> sum of volume * (pyramid centroid - vertex)
        for t in &self.tet_vid {
            // get circumcenter
            let c_t = circumsphere3(pts_pos[t[0]], pts_pos[t[1]], pts_pos[t[2]], pts_pos[t[3]]);

            for (loc_i, loc_j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                // get edge (i, j) and the other two vertices
                let vid_i = t[loc_i];
                let vid_j = t[loc_j];
                let mut other = (0..4).filter(|&l| l != loc_i && l != loc_j).map(|l| t[l]);
                let (vid_k, vid_l) = (other.next().unwrap(), other.next().unwrap());
                let (pi, pj) = (pts_pos[vid_i], pts_pos[vid_j]);
                let dist = len3(sub3(pj, pi));
                let e = scale3(1.0 / dist, sub3(pj, pi));
                let m = scale3(0.5, add3(pi, pj));

                for (vid_a, vid_b) in [(vid_k, vid_l), (vid_l, vid_k)] {
                    // in-plane direction toward a and face normal toward b
                    let (pa, pb) = (pts_pos[vid_a], pts_pos[vid_b]);
                    let (c_f, _) = circumcircle3(pi, pj, pa);
                    let u = sub3(sub3(pa, m), scale3(dot3(sub3(pa, m), e), e));
                    let u = scale3(1.0 / len3(u), u);
                    let n = cross3(e, u);
                    let n = if dot3(sub3(pb, m), n) < 0.0 { scale3(-1.0, n) } else { n };
                    let s1 = dot3(sub3(c_f, m), u);
                    let s2 = dot3(sub3(c_t, c_f), n);
                    let area = 0.5 * s1 * s2;

                    // accumulate face pieces; the piece center is (m + c_f + c_t) / 3
                    let cent = scale3(1.0 / 3.0, add3(m, add3(c_f, c_t)));
                    let entry = edge_face.entry((vid_i.min(vid_j), vid_i.max(vid_j))).or_insert((0.0, (0.0, 0.0, 0.0)));
                    entry.0 += area;
                    entry.1 = add3(entry.1, scale3(area, sub3(cent, m)));

                    // accumulate pyramid centroids; each has signed volume area * dist / 6 and centroid 3/4 of the way to the piece center
                    let vol = area * dist / 6.0;
                    for (vid_v, pv) in [(vid_i, pi), (vid_j, pj)] {
                        cent_sum[vid_v] = add3(cent_sum[vid_v], scale3(vol * 0.75, sub3(cent, pv)));
                    }
                }
            }
        }

        // step 3: neighbor, face, and control volume data
        // edges are sorted so that neighbor order does not depend on hash order
        let mut edge_list: Vec<_> = edge_face.into_iter().collect();
        edge_list.sort_unstable_by_key(|&(key, _)| key);
        for ((vid_a, vid_b), (area, skew_sum)) in edge_list {
            // get edge geometry
            let (pa, pb) = (pts_pos[vid_a], pts_pos[vid_b]);
            let dist = len3(sub3(pb, pa));
            let norm = scale3(1.0 / dist, sub3(pb, pa));

            // get face skewness; zero-area faces carry no flux
            let skew = if area.abs() > 1e-12 * dist * dist { scale3(1.0 / area, skew_sum) } else { (0.0, 0.0, 0.0) };

            // store for both directions; normals point out of each control volume
            self.ngh_vid[vid_a].push(vid_b);
            self.ngh_dist[vid_a].push(dist);
            self.face_area[vid_a].push(area);
            self.face_norm[vid_a].push(norm);
            self.face_skew[vid_a].push(skew);
            self.ngh_vid[vid_b].push(vid_a);
            self.ngh_dist[vid_b].push(dist);
            self.face_area[vid_b].push(area);
            self.face_norm[vid_b].push(scale3(-1.0, norm));
            self.face_skew[vid_b].push(skew);

            // control volume volume: each side of the face covers a pyramid of volume area * (0.5 * dist) / 3
            self.cvol_vol[vid_a] += area * dist / 6.0;
            self.cvol_vol[vid_b] += area * dist / 6.0;
        }

        // control volume skewness
        for vid in 0..dom.num_tot() {
            self.cvol_skew[vid] = scale3(1.0 / self.cvol_vol[vid], cent_sum[vid]);
        }

        // step 4: boundary sub-faces

        // the sub-face of vertex i in boundary face (i, j, k) is made of the signed triangles (i, m_ij, c) and (i, m_ik, c)
        // each has area 0.25 * |ij| * s_ij, where s_ij = (c - m_ij) . u_ij and u_ij points from edge ij toward k in the plane
        // centroids are kept in barycentric coordinates (w_i, w_j, w_k) for interpolation
        for t in &surf_vid {
            // get triangle geometry
            let (p0, p1, p2) = (pts_pos[t[0]], pts_pos[t[1]], pts_pos[t[2]]);
            let norm = cross3(sub3(p1, p0), sub3(p2, p0));
            let norm = scale3(1.0 / len3(norm), norm);
            let (c, _) = circumcircle3(p0, p1, p2);

            // barycentric coordinates of the circumcenter; lam_i is proportional to a^2 (b^2 + c^2 - a^2) with a opposite i
            let sq = [dot3(sub3(p2, p1), sub3(p2, p1)), dot3(sub3(p0, p2), sub3(p0, p2)), dot3(sub3(p1, p0), sub3(p1, p0))];
            let lam_raw = [sq[0] * (sq[1] + sq[2] - sq[0]), sq[1] * (sq[2] + sq[0] - sq[1]), sq[2] * (sq[0] + sq[1] - sq[2])];
            let lam_sum = lam_raw[0] + lam_raw[1] + lam_raw[2];
            let lam = lam_raw.map(|v| v / lam_sum);

            for loc in 0..3 {
                // get vertex i and the other two vertices in counterclockwise order
                let (vid_i, vid_j, vid_k) = (t[loc], t[(loc + 1) % 3], t[(loc + 2) % 3]);
                let (pi, pj, pk) = (pts_pos[vid_i], pts_pos[vid_j], pts_pos[vid_k]);

                // signed areas of the two triangles
                let mut area = [0.0; 2];
                for (side, (pa, pb)) in [(pj, pk), (pk, pj)].into_iter().enumerate() {
                    let m = scale3(0.5, add3(pi, pa));
                    let e = sub3(pa, pi);
                    let len = len3(e);
                    let u = sub3(sub3(pb, m), scale3(dot3(sub3(pb, m), e) / (len * len), e));
                    area[side] = 0.25 * len * dot3(sub3(c, m), scale3(1.0 / len3(u), u));
                }
                let area_sum = area[0] + area[1];

                // centroid weights; (i + m_ij + c) / 3 and (i + m_ik + c) / 3 in barycentric coordinates (i, j, k)
                let lam_loc = (lam[loc], lam[(loc + 1) % 3], lam[(loc + 2) % 3]);
                let cent_j = ((1.5 + lam_loc.0) / 3.0, (0.5 + lam_loc.1) / 3.0, lam_loc.2 / 3.0);
                let cent_k = ((1.5 + lam_loc.0) / 3.0, lam_loc.1 / 3.0, (0.5 + lam_loc.2) / 3.0);
                let wgt = if area_sum.abs() > 1e-12 * len3(sub3(pj, pi)) * len3(sub3(pk, pi)) {
                    scale3(1.0 / area_sum, add3(scale3(area[0], cent_j), scale3(area[1], cent_k)))
                } else {
                    (1.0, 0.0, 0.0)
                };

                // store sub-face
                self.bnd_sum[vid_i] += area_sum;
                self.bnd_vid[vid_i].push((vid_j, vid_k));
                self.bnd_area[vid_i].push(area_sum);
                self.bnd_norm[vid_i].push(norm);
                self.bnd_wgt[vid_i].push(wgt);
            }
        }

        // step 5: least-squares gradient weights
        // minimize sum(omega_n * (u_c + grad(u) . dx_n - u_n)^2) with omega_n = 1 / d_n^2
        // grad(u) = M^-1 sum(omega_n * dx_n * (u_n - u_c)), where M = sum(omega_n * dx_n dx_n^T)
        for vid_c in 0..dom.num_tot() {
            // get center point
            let pc = pts_pos[vid_c];

            // construct moment matrix M
            let mut m = [[0.0; 3]; 3];
            for &vid_n in &self.ngh_vid[vid_c] {
                let d = sub3(pts_pos[vid_n], pc);
                let d = [d.0, d.1, d.2];
                let omega = 1.0 / (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                for r in 0..3 { for s in 0..3 { m[r][s] += omega * d[r] * d[s]; } }
            }

            // invert M by the adjugate; omega * dx^2 <= 1, so M is dimensionless
            let adj = [
                [m[1][1] * m[2][2] - m[1][2] * m[2][1], m[0][2] * m[2][1] - m[0][1] * m[2][2], m[0][1] * m[1][2] - m[0][2] * m[1][1]],
                [m[1][2] * m[2][0] - m[1][0] * m[2][2], m[0][0] * m[2][2] - m[0][2] * m[2][0], m[0][2] * m[1][0] - m[0][0] * m[1][2]],
                [m[1][0] * m[2][1] - m[1][1] * m[2][0], m[0][1] * m[2][0] - m[0][0] * m[2][1], m[0][0] * m[1][1] - m[0][1] * m[1][0]],
            ];
            let det = m[0][0] * adj[0][0] + m[0][1] * adj[1][0] + m[0][2] * adj[2][0];
            if !(det > 1e-12 * (m[0][0] + m[1][1] + m[2][2]).powi(3)) { return Err(VFVChemError::TetraSingularGradient { dom_id: dom.dom_id, vid: vid_c }) }

            // get weights w_n = omega_n * M^-1 dx_n
            for &vid_n in &self.ngh_vid[vid_c] {
                let d = sub3(pts_pos[vid_n], pc);
                let omega = 1.0 / dot3(d, d);
                let w = [0, 1, 2].map(|r| omega * (adj[r][0] * d.0 + adj[r][1] * d.1 + adj[r][2] * d.2) / det);
                self.grad_wgt[vid_c].push((w[0], w[1], w[2]));
            }
        }

        // step 6: pressure gradient coefficients
        // p_f = (p_v + p_n) / 2 + (grad(p_v) + grad(p_n)) / 2 . skew
        // p_b = w_v * p_v + w_j * p_j + w_k * p_k
        // grad(p_v) = sum_ngh w * (p_n - p_v), so the skewness correction is fully implicit
        for vid_v in 0..dom.num_tot() {
            // iterate through faces
            let mut coeff_m = Vec::new();
            for (loc, &vid_n) in self.ngh_vid[vid_v].iter().enumerate() {
                // get geometry
                let area = self.face_area[vid_v][loc];
                let (nx, ny, nz) = self.face_norm[vid_v][loc];
                let (sx, sy, sz) = self.face_skew[vid_v][loc];

                // arithmetic mean
                coeff_m.push((vid_v, 0.5 * area * nx, 0.5 * area * ny, 0.5 * area * nz));
                coeff_m.push((vid_n, 0.5 * area * nx, 0.5 * area * ny, 0.5 * area * nz));

                // skewness correction; (p_m - p_s) terms for s = v, n
                for vid_s in [vid_v, vid_n] {
                    for (loc_m, &vid_m) in self.ngh_vid[vid_s].iter().enumerate() {
                        // project gradient weight onto skewness vector
                        let (w_x, w_y, w_z) = self.grad_wgt[vid_s][loc_m];
                        let w_skew = 0.5 * (w_x * sx + w_y * sy + w_z * sz) * area;
                        coeff_m.push((vid_m,  w_skew * nx,  w_skew * ny,  w_skew * nz));
                        coeff_m.push((vid_s, -w_skew * nx, -w_skew * ny, -w_skew * nz));
                    }
                }
            }

            // iterate through boundary sub-faces
            for (blc, &(vid_j, vid_k)) in self.bnd_vid[vid_v].iter().enumerate() {
                // get geometry
                let area = self.bnd_area[vid_v][blc];
                let (nx, ny, nz) = self.bnd_norm[vid_v][blc];
                let (w_v, w_j, w_k) = self.bnd_wgt[vid_v][blc];

                // linear interpolation to the sub-face centroid
                coeff_m.push((vid_v, w_v * area * nx, w_v * area * ny, w_v * area * nz));
                coeff_m.push((vid_j, w_j * area * nx, w_j * area * ny, w_j * area * nz));
                coeff_m.push((vid_k, w_k * area * nx, w_k * area * ny, w_k * area * nz));
            }

            // merge repeated vertices
            coeff_m.sort_by_key(|&(vid_m, _, _, _)| vid_m);
            coeff_m.dedup_by(|cur, prev| {
                if cur.0 != prev.0 { return false; }
                prev.1 += cur.1;
                prev.2 += cur.2;
                prev.3 += cur.3;
                true
            });
            self.grad_coeff[vid_v] = coeff_m;
        }

        Ok(())
    }

}
