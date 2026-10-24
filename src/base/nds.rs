use crate::base::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

const OFFSET_RATIO: f64 = 0.25;
const FACE_TRY: usize = 30;

#[derive(Default)]
pub struct Nodes {
    // spacing
    pub space: f64,

    // nodes
    // nodes are on boundaries and interfaces; every node lies inside one face
    pub nds_pos: Vec<(f64, f64, f64)>,  // gnd (global node id) -> (x, y, z)

    // 2d groups
    // defines boundaries and interfaces
    pub grp2d_gnd: HashMap<String, Vec<usize>>,  // grp2d_id -> [gnd, ...]

    // 3d groups
    // defines domains, outward normals, and the closed surface around each domain
    pub grp3d_gnd: HashMap<String, Vec<usize>>,  // grp3d_id -> [gnd, ...]
    pub grp3d_norm: HashMap<String, Vec<(f64, f64, f64)>>,  // grp3d_id -> [(nx, ny, nz), ...]; same order as grp3d_gnd
    pub grp3d_tri: HashMap<String, Vec<[usize; 3]>>,  // grp3d_id -> [[gnd, gnd, gnd], ...]; counterclockwise seen from outside
}

struct Face {
    // a face is a set of coplanar polygons used by the same groups, merged across their shared edges
    loops: Vec<Vec<usize>>,  // lid (loop id) -> [vertex id, ...]; the face is on the left of each loop seen against norm
    norm: (f64, f64, f64),  // unit normal by the right-hand rule
    corner_gnd: Vec<Vec<usize>>,  // lid, loc (loop vertex) -> inset corner node
    edge_gnd: Vec<Vec<Vec<usize>>>,  // lid, loc -> rim nodes between the corners of loc and loc + 1
    gnd: Vec<usize>,  // all nodes of the face
    tri: Vec<[usize; 3]>,  // face triangles in gnd; counterclockwise about norm
}

impl Nodes {
    pub fn new_from_obj(file_path: &str, h: f64) -> Result<Nodes, VFVChemError> {
        if !(h > 0.0) { return Err(VFVChemError::NodesNonPositiveSpacing { h }) }
        let nds = read_nds_obj(file_path, h, OFFSET_RATIO)?;
        Ok(nds)
    }

    pub fn new_from_bounds(x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64, h: f64) -> Result<Nodes, VFVChemError> {
        // creates boundary nodes on the faces of the bounding box
        // the box is passed through the same face pipeline as an obj file
        // boundaries: bnd_l and bnd_r (x), bnd_f and bnd_k (y; k for back), bnd_b and bnd_t (z)

        // input checking
        if !(x_max > x_min) { return Err(VFVChemError::NodesInvalidBoundsX { x_min, x_max }) }
        if !(y_max > y_min) { return Err(VFVChemError::NodesInvalidBoundsY { y_min, y_max }) }
        if !(z_max > z_min) { return Err(VFVChemError::NodesInvalidBoundsZ { z_min, z_max }) }
        if !(h > 0.0) { return Err(VFVChemError::NodesNonPositiveSpacing { h }) }

        // box vertices; vertex id = i + 2 * j + 4 * k for (x, y, z) = (x_i, y_j, z_k)
        let mut verts = Vec::with_capacity(8);
        for k in 0..2 {
        for j in 0..2 {
        for i in 0..2 {
            let x = if i == 0 { x_min } else { x_max };
            let y = if j == 0 { y_min } else { y_max };
            let z = if k == 0 { z_min } else { z_max };
            verts.push((x, y, z));
        }}}

        // faces with outward normals
        let face_l = vec![0, 4, 6, 2];  // left (x_min)
        let face_r = vec![1, 3, 7, 5];  // right (x_max)
        let face_f = vec![0, 1, 5, 4];  // front (y_min)
        let face_k = vec![2, 6, 7, 3];  // back (y_max)
        let face_b = vec![0, 2, 3, 1];  // bottom (z_min)
        let face_t = vec![4, 5, 7, 6];  // top (z_max)
        let groups = vec![
            ObjGroup { id: "bnd_l".to_string(), is_dom: false, faces: vec![face_l.clone()] },
            ObjGroup { id: "bnd_r".to_string(), is_dom: false, faces: vec![face_r.clone()] },
            ObjGroup { id: "bnd_f".to_string(), is_dom: false, faces: vec![face_f.clone()] },
            ObjGroup { id: "bnd_k".to_string(), is_dom: false, faces: vec![face_k.clone()] },
            ObjGroup { id: "bnd_b".to_string(), is_dom: false, faces: vec![face_b.clone()] },
            ObjGroup { id: "bnd_t".to_string(), is_dom: false, faces: vec![face_t.clone()] },
            ObjGroup { id: "dom".to_string(), is_dom: true, faces: vec![face_l, face_r, face_f, face_k, face_b, face_t] },
        ];

        Nodes::new_from_faces(&verts, &groups, h, OFFSET_RATIO)
    }

    pub fn new_from_faces(verts: &[(f64, f64, f64)], groups: &[ObjGroup], h: f64, offset_ratio: f64) -> Result<Nodes, VFVChemError> {
        // places nodes on every face and closes the surface around every domain group
        // coplanar polygons used by the same groups are merged into one face, which may have holes
        // nodes are inset from the face edges by dh, so every node lies inside exactly one face
        // the gaps along the face edges are closed by seam ribbons and the gaps at the polygon vertices by corner fans
        // faces sharing an edge stagger their rim nodes along it so that no four seam nodes are cocircular

        // default struct
        let mut nds = Nodes::default();
        nds.space = h;
        let dh = offset_ratio * h;

        // step 1: unique polygons
        // a polygon is identified by its sorted vertex ids; groups may use either winding
        let mut poly_key: HashMap<Vec<usize>, usize> = HashMap::new();  // sorted vertex ids -> pid (polygon id)
        let mut poly_vid: Vec<Vec<usize>> = Vec::new();  // pid -> vertex ids
        let mut grp_poly: Vec<Vec<(usize, bool)>> = Vec::new();  // group -> [(pid, true if same winding as poly_vid), ...]
        for grp in groups {
            let mut list = Vec::new();
            let mut seen = HashSet::new();
            for (loc, vids) in grp.faces.iter().enumerate() {
                let mut key = vids.clone();
                key.sort_unstable();
                if key.windows(2).any(|w| w[0] == w[1]) { return Err(_geom(format!("face {loc} of group '{}' repeats a vertex", grp.id))) }
                let pid = *poly_key.entry(key).or_insert_with(|| { poly_vid.push(vids.clone()); poly_vid.len() - 1 });
                if !seen.insert(pid) { return Err(_geom(format!("group '{}' uses face {loc} twice", grp.id))) }
                match _same_winding(&poly_vid[pid], vids) {
                    Some(same) => list.push((pid, same)),
                    None => return Err(_geom(format!("face {loc} of group '{}' uses the vertices of another face in a different order", grp.id))),
                }
            }
            grp_poly.push(list);
        }

        // step 2: merge coplanar polygons into faces

        // polygon normals and signatures (groups and windings that use the polygon)
        let mut poly_sig = vec![Vec::new(); poly_vid.len()];
        for (idx, list) in grp_poly.iter().enumerate() {
            for &(pid, same) in list { poly_sig[pid].push((idx, same)); }
        }
        let mut poly_norm = Vec::with_capacity(poly_vid.len());
        for (pid, vids) in poly_vid.iter().enumerate() {
            let norm = _newell(&vids.iter().map(|&v| verts[v]).collect::<Vec<_>>());
            if !(len3(norm) > 0.0) { return Err(_geom(format!("polygon {pid} has zero area"))) }
            poly_norm.push(scale3(1.0 / len3(norm), norm));
        }

        // join polygons that share an edge in opposite directions, have the same signature, and lie in the same plane
        // edges shared by more than two polygons are junctions and are never joined across
        let mut edge_poly: HashMap<(usize, usize), Vec<usize>> = HashMap::new();  // (vid_min, vid_max) -> [pid, ...]
        for (pid, vids) in poly_vid.iter().enumerate() {
            for loc in 0..vids.len() {
                let (a, b) = (vids[loc], vids[(loc + 1) % vids.len()]);
                edge_poly.entry((a.min(b), a.max(b))).or_default().push(pid);
            }
        }
        let mut root: Vec<usize> = (0..poly_vid.len()).collect();  // pid -> parent pid (union-find)
        for (&(a, b), list) in &edge_poly {
            if list.len() != 2 { continue; }
            let (p, q) = (list[0], list[1]);
            if poly_sig[p] != poly_sig[q] || dot3(poly_norm[p], poly_norm[q]) < 1.0 - 1e-12 { continue; }
            if _has_edge(&poly_vid[p], a, b) == _has_edge(&poly_vid[q], a, b) { continue; }
            let size = len3(sub3(verts[b], verts[a]));
            if poly_vid[q].iter().any(|&v| dot3(sub3(verts[v], verts[a]), poly_norm[p]).abs() > 1e-9 * size) { continue; }
            let (rp, rq) = (_root(&mut root, p), _root(&mut root, q));
            root[rp.max(rq)] = rp.min(rq);
        }

        // faces and their boundary loops
        // edges used in both directions inside a face are removed; the rest are chained into loops
        let mut face_poly: Vec<Vec<usize>> = Vec::new();  // fid (face id) -> [pid, ...]
        let mut poly_face = vec![0; poly_vid.len()];  // pid -> fid
        let mut root_face = HashMap::new();  // root pid -> fid
        for pid in 0..poly_vid.len() {
            let r = _root(&mut root, pid);
            let fid = *root_face.entry(r).or_insert_with(|| { face_poly.push(Vec::new()); face_poly.len() - 1 });
            face_poly[fid].push(pid);
            poly_face[pid] = fid;
        }
        let mut face_loops = Vec::with_capacity(face_poly.len());
        for (fid, pids) in face_poly.iter().enumerate() {
            let mut edges = HashSet::new();
            for &pid in pids {
                let vids = &poly_vid[pid];
                for loc in 0..vids.len() { edges.insert((vids[loc], vids[(loc + 1) % vids.len()])); }
            }
            let mut next = HashMap::new();  // vertex -> next vertex along the boundary
            for &(a, b) in &edges {
                if edges.contains(&(b, a)) { continue; }
                if next.insert(a, b).is_some() { return Err(_geom(format!("face {fid} touches itself at a vertex"))) }
            }
            let mut starts: Vec<usize> = next.keys().copied().collect();
            starts.sort_unstable();
            let mut used = HashSet::new();
            let mut loops = Vec::new();
            for start in starts {
                if used.contains(&start) { continue; }
                let mut lp = Vec::new();
                let mut v = start;
                loop {
                    used.insert(v);
                    lp.push(v);
                    v = next[&v];
                    if v == start { break; }
                }
                loops.push(lp);
            }
            face_loops.push(loops);
        }

        // faces used by each group
        let grp_face: Vec<Vec<(usize, bool)>> = grp_poly.iter().map(|list| {
            let mut seen = HashSet::new();
            list.iter().filter(|&&(pid, _)| seen.insert(poly_face[pid])).map(|&(pid, same)| (poly_face[pid], same)).collect()
        }).collect();

        // step 3: faces around each loop edge; used to stagger the rim nodes
        let mut edge_face: HashMap<(usize, usize), Vec<usize>> = HashMap::new();  // (vid_min, vid_max) -> [fid, ...]
        for (fid, loops) in face_loops.iter().enumerate() {
            for lp in loops {
                for loc in 0..lp.len() {
                    let (a, b) = (lp[loc], lp[(loc + 1) % lp.len()]);
                    edge_face.entry((a.min(b), a.max(b))).or_default().push(fid);
                }
            }
        }

        // step 4: nodes and triangles of each face
        let mut faces = Vec::with_capacity(face_loops.len());
        for (fid, loops) in face_loops.into_iter().enumerate() {
            faces.push(_face_nodes(&mut nds, verts, loops, fid, &edge_face, h, dh)?);
        }

        // step 5: groups
        for (grp, list) in groups.iter().zip(grp_face.iter()) {
            if list.is_empty() { continue; }

            // boundary and interface groups hold the nodes of their faces
            if !grp.is_dom {
                let mut gnds = Vec::new();
                for &(fid, _) in list { gnds.extend(faces[fid].gnd.iter().copied()); }
                nds.grp2d_gnd.insert(grp.id.clone(), gnds);
                continue;
            }

            // domain nodes, outward normals, and face triangles
            let mut gnds = Vec::new();
            let mut norms = Vec::new();
            let mut tris = Vec::new();
            for &(fid, same) in list {
                let face = &faces[fid];
                let sgn = if same { 1.0 } else { -1.0 };
                for &gnd in &face.gnd {
                    gnds.push(gnd);
                    norms.push(scale3(sgn, face.norm));
                }
                for &[a, b, c] in &face.tri { tris.push(if same { [a, b, c] } else { [a, c, b] }); }
            }

            // directed loop edges in the winding of the group
            // a closed, consistently oriented surface uses every edge once in each direction
            let mut edge_own: HashMap<(usize, usize), usize> = HashMap::new();  // (vid_a, vid_b) -> fid
            for &(fid, same) in list {
                for lp in _group_winding(&faces[fid], same) {
                    for loc in 0..lp.len() {
                        let (a, b) = (lp[loc], lp[(loc + 1) % lp.len()]);
                        if edge_own.insert((a, b), fid).is_some() { return Err(_geom(format!("faces of group '{}' are not consistently oriented or an edge is shared by more than two faces", grp.id))) }
                    }
                }
            }
            for &(a, b) in edge_own.keys() {
                if !edge_own.contains_key(&(b, a)) { return Err(_geom(format!("surface of group '{}' is not closed", grp.id))) }
            }

            // seam ribbons
            // walking the chain of the face that owns a -> b from a to b, the ribbon is on the right seen from outside
            for &(fid_p, same) in list {
                for lp in _group_winding(&faces[fid_p], same) {
                    for loc in 0..lp.len() {
                        let (a, b) = (lp[loc], lp[(loc + 1) % lp.len()]);
                        if a > b { continue; }
                        let fid_q = edge_own[&(b, a)];
                        let chain_p = _chain(&faces[fid_p], a, b);
                        let chain_q = _chain(&faces[fid_q], a, b);
                        _ribbon(&nds, &chain_p, &chain_q, &mut tris);
                    }
                }
            }

            // corner fans
            // the faces around a loop vertex are visited through the edges leaving the vertex
            let face_same: HashMap<usize, bool> = list.iter().copied().collect();  // fid -> same winding
            let mut visited = HashSet::new();  // (vid, fid)
            for &(fid_0, same_0) in list {
                for lp_0 in _group_winding(&faces[fid_0], same_0) {
                for &v in &lp_0 {
                    if visited.contains(&(v, fid_0)) { continue; }
                    let mut corners = Vec::new();
                    let mut norm_sum = (0.0, 0.0, 0.0);
                    let mut fid = fid_0;
                    loop {
                        visited.insert((v, fid));
                        let face = &faces[fid];
                        let same = face_same[&fid];
                        let (lid, loc) = _find_vertex(face, v);
                        let lp = &face.loops[lid];
                        let w = if same { lp[(loc + 1) % lp.len()] } else { lp[(loc + lp.len() - 1) % lp.len()] };
                        corners.push(face.corner_gnd[lid][loc]);
                        norm_sum = add3(norm_sum, scale3(if same { 1.0 } else { -1.0 }, face.norm));
                        fid = edge_own[&(w, v)];
                        if fid == fid_0 || corners.len() > list.len() { break; }
                    }
                    _fan(&nds, &corners, norm_sum, &mut tris);
                }}
            }

            // outward normals enclose a positive volume
            let mut vol = 0.0;
            for &[a, b, c] in &tris { vol += dot3(nds.nds_pos[a], cross3(nds.nds_pos[b], nds.nds_pos[c])) / 6.0; }
            if !(vol > 0.0) { return Err(_geom(format!("faces of group '{}' must point out of the domain", grp.id))) }

            nds.grp3d_gnd.insert(grp.id.clone(), gnds);
            nds.grp3d_norm.insert(grp.id.clone(), norms);
            nds.grp3d_tri.insert(grp.id.clone(), tris);
        }

        Ok(nds)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_pos.len()
    }
}

fn _face_nodes(nds: &mut Nodes, verts: &[(f64, f64, f64)], loops: Vec<Vec<usize>>, fid: usize, edge_face: &HashMap<(usize, usize), Vec<usize>>, h: f64, dh: f64) -> Result<Face, VFVChemError> {
    // places rim, corner, and interior nodes on one planar face and triangulates them
    // rim nodes run parallel to each loop edge at a distance dh inside the face; corner nodes join the rims of adjacent edges
    // interior nodes are placed by uniform poisson-disk sampling at spacing h

    // step 1: plane of the face
    let pos: Vec<Vec<(f64, f64, f64)>> = loops.iter().map(|lp| lp.iter().map(|&v| verts[v]).collect()).collect();
    let mut norm = (0.0, 0.0, 0.0);
    for lp_pos in &pos { norm = add3(norm, _newell(lp_pos)); }
    if !(len3(norm) > 0.0) { return Err(_geom(format!("face {fid} has zero area"))) }
    let area = 0.5 * len3(norm);
    let norm = scale3(1.0 / len3(norm), norm);

    // limit the inset of small faces to half the inradius (area / perimeter for a triangle), so the inset polygon stays valid
    let perim: f64 = pos.iter().map(|lp_pos| (0..lp_pos.len()).map(|loc| len3(sub3(lp_pos[(loc + 1) % lp_pos.len()], lp_pos[loc]))).sum::<f64>()).sum();
    let dh = dh.min(area / perim);

    // check planarity against the size of the face
    let p0 = pos[0][0];
    let size = pos.iter().flatten().map(|&p| len3(sub3(p, p0))).fold(0.0, f64::max);
    if pos.iter().flatten().any(|&p| dot3(sub3(p, p0), norm).abs() > 1e-6 * size) { return Err(_geom(format!("face {fid} is not planar"))) }

    // local frame; e1 along the first edge, e2 = n x e1, so the face is on the left of each loop in (e1, e2)
    let d01 = sub3(pos[0][1], p0);
    let e1 = sub3(d01, scale3(dot3(d01, norm), norm));
    let e1 = scale3(1.0 / len3(e1), e1);
    let e2 = cross3(norm, e1);
    let to_2d = |p: (f64, f64, f64)| (dot3(sub3(p, p0), e1), dot3(sub3(p, p0), e2));
    let to_3d = |(u, v): (f64, f64)| add3(p0, add3(scale3(u, e1), scale3(v, e2)));

    // step 2: corner and rim nodes of each loop
    // the inset lines of the edges before and after a vertex meet at q + dh * (n_a + n_b) / (1 + n_a . n_b)
    // edges are split into intervals no longer than h; faces sharing an edge shift their nodes by rank / count of an interval
    let mut corner = Vec::with_capacity(loops.len());  // lid, loc -> (u, v)
    let mut rim = Vec::with_capacity(loops.len());  // lid, loc -> [(u, v), ...]
    for (lid, lp) in loops.iter().enumerate() {
        let num = lp.len();
        let q: Vec<(f64, f64)> = pos[lid].iter().map(|&p| to_2d(p)).collect();
        let mut tang = Vec::with_capacity(num);  // loc -> unit tangent of edge (loc, loc + 1)
        let mut lens = Vec::with_capacity(num);  // loc -> edge length
        for loc in 0..num {
            let (x0, y0) = q[loc];
            let (x1, y1) = q[(loc + 1) % num];
            let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
            if len <= 2.0 * dh { return Err(_geom(format!("edge of face {fid} must be longer than two corner offsets. Got length = {len} and 2*dh = {}", 2.0 * dh))) }
            tang.push(((x1 - x0) / len, (y1 - y0) / len));
            lens.push(len);
        }
        let mut corner_lp = Vec::with_capacity(num);
        for loc in 0..num {
            let (ta, tb) = (tang[(loc + num - 1) % num], tang[loc]);
            let (na, nb) = ((-ta.1, ta.0), (-tb.1, tb.0));  // inward normals; the face is on the left of each edge
            let den = 1.0 + na.0 * nb.0 + na.1 * nb.1;
            if den < 1e-6 { return Err(_geom(format!("face {fid} folds back on itself at a vertex"))) }
            corner_lp.push((q[loc].0 + dh * (na.0 + nb.0) / den, q[loc].1 + dh * (na.1 + nb.1) / den));
        }
        let mut rim_lp = Vec::with_capacity(num);
        for loc in 0..num {
            let (a, b) = (lp[loc], lp[(loc + 1) % num]);
            let list = &edge_face[&(a.min(b), a.max(b))];
            let rank = list.iter().position(|&f| f == fid).unwrap();
            let shift = rank as f64 / list.len() as f64;
            let (tx, ty) = tang[loc];
            let (nx, ny) = (-ty, tx);
            let cnt = ((lens[loc] / h - 1e-6).ceil() as usize).max(1);
            let ds = lens[loc] / (cnt as f64);

            // keep rim nodes between the corner nodes of the edge
            let (ca, cb) = (corner_lp[loc], corner_lp[(loc + 1) % num]);
            let s_a = (ca.0 - q[loc].0) * tx + (ca.1 - q[loc].1) * ty;
            let s_b = (cb.0 - q[loc].0) * tx + (cb.1 - q[loc].1) * ty;
            let mut rim_loc = Vec::new();
            for k in 0..cnt {
                let s = (k as f64 + shift) * ds;
                if s < s_a + 2.0 * dh || s > s_b - 2.0 * dh { continue; }
                rim_loc.push((q[loc].0 + s * tx + dh * nx, q[loc].1 + s * ty + dh * ny));
            }
            rim_lp.push(rim_loc);
        }
        corner.push(corner_lp);
        rim.push(rim_lp);
    }

    // step 3: interior nodes
    let mut loop_pos = Vec::with_capacity(loops.len());  // lid -> [corner, rim, corner, rim, ...]
    for lid in 0..loops.len() {
        let mut lp_pos = Vec::new();
        for loc in 0..loops[lid].len() {
            lp_pos.push(corner[lid][loc]);
            lp_pos.extend(rim[lid][loc].iter().copied());
        }
        loop_pos.push(lp_pos);
    }
    let mut rng = StdRng::seed_from_u64(fid as u64);
    let inner = _face_sample(&loop_pos, h, &mut rng);

    // step 4: register nodes
    let mut face = Face { loops, norm, corner_gnd: Vec::new(), edge_gnd: Vec::new(), gnd: Vec::new(), tri: Vec::new() };
    let mut pos_2d = Vec::new();  // face node index -> (u, v)
    let mut loop_idx = Vec::new();  // lid -> [face node index, ...]
    for lid in 0..face.loops.len() {
        let mut corner_gnd = Vec::new();
        let mut edge_gnd = Vec::new();
        let mut idx_lp = Vec::new();
        for loc in 0..face.loops[lid].len() {
            idx_lp.push(pos_2d.len());
            corner_gnd.push(nds.nds_pos.len());
            face.gnd.push(nds.nds_pos.len());
            nds.nds_pos.push(to_3d(corner[lid][loc]));
            pos_2d.push(corner[lid][loc]);
            let mut edge_loc = Vec::new();
            for &uv in &rim[lid][loc] {
                idx_lp.push(pos_2d.len());
                edge_loc.push(nds.nds_pos.len());
                face.gnd.push(nds.nds_pos.len());
                nds.nds_pos.push(to_3d(uv));
                pos_2d.push(uv);
            }
            edge_gnd.push(edge_loc);
        }
        face.corner_gnd.push(corner_gnd);
        face.edge_gnd.push(edge_gnd);
        loop_idx.push(idx_lp);
    }
    for &uv in &inner {
        face.gnd.push(nds.nds_pos.len());
        nds.nds_pos.push(to_3d(uv));
        pos_2d.push(uv);
    }

    // step 5: constrained Delaunay triangulation of the face
    // same labeling as the 2D domains: loops keep the face on their left, so a counterclockwise triangle containing a loop edge is inside
    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
    let mut hdl_idx = vec![usize::MAX; pos_2d.len()];  // vertex handle index -> face node index
    let mut idx_hdl = Vec::with_capacity(pos_2d.len());  // face node index -> vertex handle
    for (idx, &(u, v)) in pos_2d.iter().enumerate() {
        let hdl = cdt.insert(Point2::new(u, v)).map_err(|source| VFVChemError::FaceInsertionError { face_id: fid, source })?;
        if hdl.index() >= hdl_idx.len() { hdl_idx.resize(hdl.index() + 1, usize::MAX); }
        hdl_idx[hdl.index()] = idx;
        idx_hdl.push(hdl);
    }
    let mut edge_loop = HashSet::new();  // directed loop edges (idx_0, idx_1)
    for idx_lp in &loop_idx {
        for loc in 0..idx_lp.len() {
            let (idx_0, idx_1) = (idx_lp[loc], idx_lp[(loc + 1) % idx_lp.len()]);
            if !cdt.can_add_constraint(idx_hdl[idx_0], idx_hdl[idx_1]) { return Err(VFVChemError::FaceConstraintIntersection { face_id: fid }) }
            cdt.add_constraint(idx_hdl[idx_0], idx_hdl[idx_1]);
            edge_loop.insert((idx_0, idx_1));
        }
    }
    let tri_all: Vec<[usize; 3]> = cdt.inner_faces().map(|f| {
        let [h0, h1, h2] = f.vertices();
        [hdl_idx[h0.fix().index()], hdl_idx[h1.fix().index()], hdl_idx[h2.fix().index()]]
    }).collect();
    let mut edge_tri = HashMap::new();  // directed edge -> triangle index
    for (idx, tri) in tri_all.iter().enumerate() {
        for loc in 0..3 { edge_tri.insert((tri[loc], tri[(loc + 1) % 3]), idx); }
    }
    let mut is_inside = vec![None; tri_all.len()];
    let mut stack = Vec::new();
    for (idx, tri) in tri_all.iter().enumerate() {
        for loc in 0..3 {
            let (i, j) = (tri[loc], tri[(loc + 1) % 3]);
            if edge_loop.contains(&(i, j)) { is_inside[idx] = Some(true); }
            else if edge_loop.contains(&(j, i)) || !edge_tri.contains_key(&(j, i)) { is_inside[idx] = Some(false); }
        }
        if is_inside[idx].is_some() { stack.push(idx); }
    }
    while let Some(idx) = stack.pop() {
        let tri = tri_all[idx];
        for loc in 0..3 {
            let (i, j) = (tri[loc], tri[(loc + 1) % 3]);
            if edge_loop.contains(&(i, j)) || edge_loop.contains(&(j, i)) { continue; }
            if let Some(&idx_ngh) = edge_tri.get(&(j, i)) {
                if is_inside[idx_ngh].is_none() {
                    is_inside[idx_ngh] = is_inside[idx];
                    stack.push(idx_ngh);
                }
            }
        }
    }
    face.tri = tri_all.into_iter().zip(is_inside).filter(|(_, is_in)| *is_in == Some(true)).map(|(t, _)| [face.gnd[t[0]], face.gnd[t[1]], face.gnd[t[2]]]).collect();

    Ok(face)
}

fn _face_sample(loop_pos: &[Vec<(f64, f64)>], h: f64, rng: &mut StdRng) -> Vec<(f64, f64)> {
    // uniform bridson poisson-disk sampling inside the loops of one face
    // the loop nodes are loaded first, so new nodes keep a distance h from them

    // bounding box and grid
    let all = || loop_pos.iter().flatten();
    let x_min = all().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let y_min = all().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let x_max = all().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let y_max = all().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let cell = h / f64::sqrt(2.0);
    let nx = (((x_max - x_min) / cell).ceil() as usize).max(1);
    let ny = (((y_max - y_min) / cell).ceil() as usize).max(1);
    let mut grid = vec![Vec::new(); nx * ny];  // (i + j*nx) -> [idx, ...]
    let cell_of = |x: f64, y: f64| ((((x - x_min) / cell).max(0.0) as usize).min(nx - 1), (((y - y_min) / cell).max(0.0) as usize).min(ny - 1));

    // load loop nodes
    let mut pts_pos = Vec::new();
    for &(x, y) in all() {
        let (i, j) = cell_of(x, y);
        grid[i + j * nx].push(pts_pos.len());
        pts_pos.push((x, y));
    }
    let num_loop = pts_pos.len();

    // candidate test: inside the loops and at least h from every node
    let is_valid = |x: f64, y: f64, pts_pos: &Vec<(f64, f64)>, grid: &Vec<Vec<usize>>| -> bool {
        if !_point_in_loops(x, y, loop_pos) { return false; }
        let (ic, jc) = cell_of(x, y);
        for j in jc.saturating_sub(2)..=(jc + 2).min(ny - 1) {
        for i in ic.saturating_sub(2)..=(ic + 2).min(nx - 1) {
            for &idx in &grid[i + j * nx] {
                let (xn, yn) = pts_pos[idx];
                if (x - xn).powi(2) + (y - yn).powi(2) < h * h { return false; }
            }
        }}
        true
    };

    // starting point: the centroid of the nodes, then a distance h left of each edge midpoint
    // the face is connected, so one starting point fills it
    let mut cand = vec![(all().map(|p| p.0).sum::<f64>() / (num_loop as f64), all().map(|p| p.1).sum::<f64>() / (num_loop as f64))];
    for lp_pos in loop_pos {
        for idx in 0..lp_pos.len() {
            let (x0, y0) = lp_pos[idx];
            let (x1, y1) = lp_pos[(idx + 1) % lp_pos.len()];
            let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
            if len == 0.0 { continue; }
            cand.push((0.5 * (x0 + x1) - h * (y1 - y0) / len, 0.5 * (y0 + y1) + h * (x1 - x0) / len));
        }
    }
    let mut pts_act = Vec::new();
    for (x, y) in cand {
        if !is_valid(x, y, &pts_pos, &grid) { continue; }
        let (i, j) = cell_of(x, y);
        grid[i + j * nx].push(pts_pos.len());
        pts_act.push(pts_pos.len());
        pts_pos.push((x, y));
        break;
    }

    // insert remaining points
    while !pts_act.is_empty() {
        let a = rng.random_range(0..pts_act.len());
        let (xa, ya) = pts_pos[pts_act[a]];
        let mut is_insert = false;
        for _ in 0..FACE_TRY {
            let ang = rng.random_range(0.0..(2.0 * PI));
            let rad = rng.random_range(h..(2.0 * h));
            let (x, y) = (xa + rad * ang.cos(), ya + rad * ang.sin());
            if !is_valid(x, y, &pts_pos, &grid) { continue; }
            let (i, j) = cell_of(x, y);
            grid[i + j * nx].push(pts_pos.len());
            pts_act.push(pts_pos.len());
            pts_pos.push((x, y));
            is_insert = true;
            break;
        }
        if !is_insert { pts_act.swap_remove(a); }
    }

    pts_pos.split_off(num_loop)
}

fn _point_in_loops(xc: f64, yc: f64, loop_pos: &[Vec<(f64, f64)>]) -> bool {
    // winding number over each closed loop; positive if inside
    let mut wind = 0;
    for lp_pos in loop_pos {
        for idx in 0..lp_pos.len() {
            let (x0, y0) = lp_pos[idx];
            let (x1, y1) = lp_pos[(idx + 1) % lp_pos.len()];
            let side = (x1 - x0) * (yc - y0) - (xc - x0) * (y1 - y0);
            if y0 <= yc && y1 > yc && side > 0.0 { wind += 1; }
            if y0 > yc && y1 <= yc && side < 0.0 { wind -= 1; }
        }
    }
    wind > 0
}

fn _ribbon(nds: &Nodes, chain_p: &[usize], chain_q: &[usize], tris: &mut Vec<[usize; 3]>) {
    // zips two parallel rim chains running from vertex a to vertex b into a strip of triangles
    // the next triangle advances the chain whose new node keeps the circumcircle of the other choice empty (Delaunay)
    // triangles (p_i, q_j, p_i+1) and (p_i, q_j, q_j+1) are counterclockwise seen from outside
    let (mut i, mut j) = (0, 0);
    let (m, n) = (chain_p.len() - 1, chain_q.len() - 1);
    while i < m || j < n {
        let adv_p = if i == m { false } else if j == n { true } else {
            let (p0, q0, p1, q1) = (nds.nds_pos[chain_p[i]], nds.nds_pos[chain_q[j]], nds.nds_pos[chain_p[i + 1]], nds.nds_pos[chain_q[j + 1]]);
            let (cen, rad) = circumcircle3(p0, q0, p1);
            len3(sub3(q1, cen)) > rad
        };
        if adv_p {
            tris.push([chain_p[i], chain_q[j], chain_p[i + 1]]);
            i += 1;
        } else {
            tris.push([chain_p[i], chain_q[j], chain_q[j + 1]]);
            j += 1;
        }
    }
}

fn _fan(nds: &Nodes, corners: &[usize], norm_sum: (f64, f64, f64), tris: &mut Vec<[usize; 3]>) {
    // closes the polygon of corner nodes around one vertex by ear clipping
    // the polygon is projected onto the plane normal to the summed outward normals of the faces around the vertex
    // and made counterclockwise about that normal; the ear with the largest smallest angle is clipped first
    if corners.len() < 3 { return; }
    let pos: Vec<(f64, f64, f64)> = corners.iter().map(|&gnd| nds.nds_pos[gnd]).collect();
    let norm = if len3(norm_sum) > 0.0 { norm_sum } else { _newell(&pos) };
    let norm = scale3(1.0 / len3(norm), norm);
    let axis = if norm.0.abs() < 0.5 { (1.0, 0.0, 0.0) } else { (0.0, 1.0, 0.0) };
    let e1 = cross3(norm, axis);
    let e1 = scale3(1.0 / len3(e1), e1);
    let e2 = cross3(norm, e1);
    let uv: Vec<(f64, f64)> = pos.iter().map(|&p| (dot3(sub3(p, pos[0]), e1), dot3(sub3(p, pos[0]), e2))).collect();
    let cross = |a: usize, b: usize, c: usize| (uv[b].0 - uv[a].0) * (uv[c].1 - uv[a].1) - (uv[b].1 - uv[a].1) * (uv[c].0 - uv[a].0);

    // counterclockwise polygon
    let mut poly: Vec<usize> = (0..corners.len()).collect();
    let area: f64 = (1..poly.len() - 1).map(|k| cross(0, k, k + 1)).sum();
    if area < 0.0 { poly.reverse(); }

    // clip ears; an ear is a convex vertex whose triangle holds no other polygon vertex
    while poly.len() > 3 {
        let num = poly.len();
        let mut best = (f64::NEG_INFINITY, 1);  // (smallest angle, position)
        for loc in 0..num {
            let (a, b, c) = (poly[(loc + num - 1) % num], poly[loc], poly[(loc + 1) % num]);
            if cross(a, b, c) <= 0.0 { continue; }
            if poly.iter().any(|&d| d != a && d != b && d != c && cross(a, b, d) >= 0.0 && cross(b, c, d) >= 0.0 && cross(c, a, d) >= 0.0) { continue; }
            let ang = |p: usize, q: usize, r: usize| {
                let (u, v) = ((uv[q].0 - uv[p].0, uv[q].1 - uv[p].1), (uv[r].0 - uv[p].0, uv[r].1 - uv[p].1));
                f64::atan2((u.0 * v.1 - u.1 * v.0).abs(), u.0 * v.0 + u.1 * v.1)
            };
            let ang_min = ang(a, b, c).min(ang(b, c, a)).min(ang(c, a, b));
            if ang_min > best.0 { best = (ang_min, loc); }
        }
        let loc = best.1;
        tris.push([corners[poly[(loc + num - 1) % num]], corners[poly[loc]], corners[poly[(loc + 1) % num]]]);
        poly.remove(loc);
    }
    tris.push([corners[poly[0]], corners[poly[1]], corners[poly[2]]]);
}

fn _newell(pos: &[(f64, f64, f64)]) -> (f64, f64, f64) {
    // area-weighted normal of a closed polygon (newell's method); twice the area in length
    let mut norm = (0.0, 0.0, 0.0);
    for loc in 0..pos.len() {
        let (x0, y0, z0) = pos[loc];
        let (x1, y1, z1) = pos[(loc + 1) % pos.len()];
        norm.0 += (y0 - y1) * (z0 + z1);
        norm.1 += (z0 - z1) * (x0 + x1);
        norm.2 += (x0 - x1) * (y0 + y1);
    }
    norm
}

fn _same_winding(a: &[usize], b: &[usize]) -> Option<bool> {
    // true if b is a rotation of a, false if b is a rotation of a reversed, None otherwise
    let n = a.len();
    if b.len() != n { return None; }
    let p = a.iter().position(|&v| v == b[0])?;
    if (0..n).all(|i| b[i] == a[(p + i) % n]) { return Some(true); }
    if (0..n).all(|i| b[i] == a[(p + n - i) % n]) { return Some(false); }
    None
}

fn _has_edge(vids: &[usize], a: usize, b: usize) -> bool {
    // true if the polygon runs from a to b
    let loc = vids.iter().position(|&v| v == a).unwrap();
    vids[(loc + 1) % vids.len()] == b
}

fn _root(root: &mut Vec<usize>, mut pid: usize) -> usize {
    // union-find root with path halving
    while root[pid] != pid {
        root[pid] = root[root[pid]];
        pid = root[pid];
    }
    pid
}

fn _find_vertex(face: &Face, v: usize) -> (usize, usize) {
    // loop and position of a vertex on the face boundary
    for (lid, lp) in face.loops.iter().enumerate() {
        if let Some(loc) = lp.iter().position(|&u| u == v) { return (lid, loc); }
    }
    unreachable!()
}

fn _group_winding(face: &Face, same: bool) -> Vec<Vec<usize>> {
    // loop vertices in the winding used by a group
    if same { face.loops.clone() } else { face.loops.iter().map(|lp| lp.iter().rev().copied().collect()).collect() }
}

fn _chain(face: &Face, a: usize, b: usize) -> Vec<usize> {
    // rim nodes of the face along edge (a, b), from the corner at a to the corner at b
    let (lid, loc_a) = _find_vertex(face, a);
    let lp = &face.loops[lid];
    let n = lp.len();
    let loc_b = lp.iter().position(|&v| v == b).unwrap();
    let mut chain = vec![face.corner_gnd[lid][loc_a]];
    if lp[(loc_a + 1) % n] == b { chain.extend(face.edge_gnd[lid][loc_a].iter().copied()); }
    else { chain.extend(face.edge_gnd[lid][loc_b].iter().rev().copied()); }
    chain.push(face.corner_gnd[lid][loc_b]);
    chain
}

fn _geom(reason: String) -> VFVChemError {
    VFVChemError::NodesInvalidGeometry { reason }
}
