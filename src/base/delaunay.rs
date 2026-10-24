use robust::{insphere, orient3d, Coord3D};
use std::collections::HashMap;

const NONE: usize = usize::MAX;

pub fn delaunay_tet(pos: &[(f64, f64, f64)]) -> Option<Vec<[usize; 4]>> {
    // Delaunay tetrahedralization of a point set by Bowyer-Watson insertion
    // returns tetrahedra as vertex indices with orient3d(p0, p1, p2, p3) > 0 (p3 below the counterclockwise face p0 p1 p2)
    // returns None if the points coincide or do not span a volume
    // all predicates are exact (robust crate), so the result does not depend on round-off

    // pseudocode:
    // enclose the points in a large super-tetrahedron
    // insert points in morton order; locate each point by a visibility walk from the last new tetrahedron
    // the cavity holds the tetrahedra whose circumsphere contains the point
    // the cavity is grown until every boundary face sees the point with a positive orientation, which only matters for cospherical ties
    // the cavity is replaced by tetrahedra joining the point to the cavity boundary
    // tetrahedra touching the super-tetrahedron are removed at the end

    // step 1: super-tetrahedron
    let num = pos.len();
    if num < 4 { return None; }
    let mut lo = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y, z) in pos {
        lo = (lo.0.min(x), lo.1.min(y), lo.2.min(z));
        hi = (hi.0.max(x), hi.1.max(y), hi.2.max(z));
    }
    let cen = scale3(0.5, add3(lo, hi));
    let size = (hi.0 - lo.0).max(hi.1 - lo.1).max(hi.2 - lo.2).max(f64::MIN_POSITIVE);
    let mut pts: Vec<Coord3D<f64>> = pos.iter().map(|&(x, y, z)| Coord3D { x, y, z }).collect();
    for dir in [(1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0)] {
        let (x, y, z) = add3(cen, scale3(50.0 * size, dir));
        pts.push(Coord3D { x, y, z });
    }
    let mut tets = vec![[num, num + 1, num + 2, num + 3]];  // tet -> [vid, ...]
    if orient3d(pts[num], pts[num + 1], pts[num + 2], pts[num + 3]) < 0.0 { tets[0].swap(0, 1); }
    let mut adj = vec![[NONE; 4]];  // tet -> [tet opposite vertex k, ...]
    let mut dead = vec![false];

    // step 2: insertion order
    // morton order of the quantized coordinates keeps consecutive points close, so walks are short
    let morton = |p: (f64, f64, f64)| -> u64 {
        let q = |v: f64, l: f64| ((v - l) / size * 1023.0).clamp(0.0, 1023.0) as u64;
        let (qx, qy, qz) = (q(p.0, lo.0), q(p.1, lo.1), q(p.2, lo.2));
        let mut code = 0;
        for b in 0..10 { code |= ((qx >> b) & 1) << (3 * b) | ((qy >> b) & 1) << (3 * b + 1) | ((qz >> b) & 1) << (3 * b + 2); }
        code
    };
    let mut order: Vec<usize> = (0..num).collect();
    order.sort_by_key(|&vid| (morton(pos[vid]), vid));

    // step 3: insert points
    let mut tet_last = 0;
    let mut turn = 0;  // rotates the first face tested by the walk to avoid cycles
    let mut in_cav = vec![false];
    for &vid in &order {
        let p = pts[vid];

        // locate the tetrahedron containing p
        let mut tet = tet_last;
        let mut num_step = 0;
        'walk: loop {
            num_step += 1;
            if num_step > 4 * tets.len() {
                // fall back to a linear search
                tet = (0..tets.len()).find(|&t| !dead[t] && (0..4).all(|k| { let mut v = tets[t].map(|u| pts[u]); v[k] = p; orient3d(v[0], v[1], v[2], v[3]) >= 0.0 }))?;
                break;
            }
            turn = (turn + 1) % 4;
            for loc in 0..4 {
                let k = (loc + turn) % 4;
                let mut v = tets[tet].map(|u| pts[u]);
                v[k] = p;
                if orient3d(v[0], v[1], v[2], v[3]) < 0.0 {
                    if adj[tet][k] == NONE { return None; }
                    tet = adj[tet][k];
                    continue 'walk;
                }
            }
            break;
        }
        if tets[tet].iter().any(|&u| pts[u].x == p.x && pts[u].y == p.y && pts[u].z == p.z) { return None; }

        // cavity of tetrahedra whose circumsphere contains p
        let mut cav = vec![tet];
        in_cav[tet] = true;
        let mut idx = 0;
        while idx < cav.len() {
            let t = cav[idx];
            idx += 1;
            for k in 0..4 {
                let t_n = adj[t][k];
                if t_n == NONE || in_cav[t_n] { continue; }
                let v = tets[t_n].map(|u| pts[u]);
                if insphere(v[0], v[1], v[2], v[3], p) > 0.0 {
                    in_cav[t_n] = true;
                    cav.push(t_n);
                }
            }
        }

        // grow the cavity until p sees every boundary face from the inside
        loop {
            let mut is_star = true;
            for idx in 0..cav.len() {
                let t = cav[idx];
                for k in 0..4 {
                    let t_n = adj[t][k];
                    if t_n != NONE && in_cav[t_n] { continue; }
                    let mut v = tets[t].map(|u| pts[u]);
                    v[k] = p;
                    if orient3d(v[0], v[1], v[2], v[3]) > 0.0 { continue; }
                    if t_n == NONE { return None; }
                    in_cav[t_n] = true;
                    cav.push(t_n);
                    is_star = false;
                }
            }
            if is_star { break; }
        }

        // replace the cavity by tetrahedra joining p to each boundary face
        // the face opposite p keeps the outside neighbor; faces containing p are matched by the edge they share with the boundary face
        let mut edge_new: HashMap<(usize, usize), (usize, usize)> = HashMap::new();  // edge -> (new tet, face index)
        for idx in 0..cav.len() {
            let t = cav[idx];
            for k in 0..4 {
                let t_n = adj[t][k];
                if t_n != NONE && in_cav[t_n] { continue; }
                let mut tet_new = tets[t];
                tet_new[k] = vid;
                let t_new = tets.len();
                tets.push(tet_new);
                adj.push([NONE; 4]);
                dead.push(false);
                in_cav.push(false);
                adj[t_new][k] = t_n;
                if t_n != NONE {
                    let m = (0..4).find(|&m| adj[t_n][m] == t).unwrap();
                    adj[t_n][m] = t_new;
                }
                for m in 0..4 {
                    if m == k { continue; }
                    let (a, b) = {
                        let mut e = (0..4).filter(|&l| l != k && l != m).map(|l| tet_new[l]);
                        let (a, b) = (e.next().unwrap(), e.next().unwrap());
                        (a.min(b), a.max(b))
                    };
                    if let Some((t_o, m_o)) = edge_new.remove(&(a, b)) {
                        adj[t_new][m] = t_o;
                        adj[t_o][m_o] = t_new;
                    } else {
                        edge_new.insert((a, b), (t_new, m));
                    }
                }
                tet_last = t_new;
            }
        }
        for &t in &cav {
            dead[t] = true;
            in_cav[t] = false;
        }
    }

    // step 4: keep tetrahedra without super vertices
    let tet_vid: Vec<[usize; 4]> = (0..tets.len()).filter(|&t| !dead[t] && tets[t].iter().all(|&u| u < num)).map(|t| tets[t]).collect();
    if tet_vid.is_empty() { return None; }
    Some(tet_vid)
}

// vector helpers shared by the geometry modules

#[inline]
pub fn add3(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 + b.0, a.1 + b.1, a.2 + b.2)
}

#[inline]
pub fn sub3(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.0 - b.0, a.1 - b.1, a.2 - b.2)
}

#[inline]
pub fn scale3(s: f64, a: (f64, f64, f64)) -> (f64, f64, f64) {
    (s * a.0, s * a.1, s * a.2)
}

#[inline]
pub fn dot3(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

#[inline]
pub fn cross3(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.1 * b.2 - a.2 * b.1, a.2 * b.0 - a.0 * b.2, a.0 * b.1 - a.1 * b.0)
}

#[inline]
pub fn len3(a: (f64, f64, f64)) -> f64 {
    dot3(a, a).sqrt()
}

pub fn circumcircle3(a: (f64, f64, f64), b: (f64, f64, f64), c: (f64, f64, f64)) -> ((f64, f64, f64), f64) {
    // center and radius of the circle through three points
    let (u, v) = (sub3(b, a), sub3(c, a));
    let w = cross3(u, v);
    let num = add3(scale3(dot3(u, u), cross3(v, w)), scale3(dot3(v, v), cross3(w, u)));
    let cen = add3(a, scale3(0.5 / dot3(w, w), num));
    (cen, len3(sub3(cen, a)))
}

pub fn circumsphere3(a: (f64, f64, f64), b: (f64, f64, f64), c: (f64, f64, f64), d: (f64, f64, f64)) -> (f64, f64, f64) {
    // center of the sphere through four points
    let (u, v, w) = (sub3(b, a), sub3(c, a), sub3(d, a));
    let num = add3(add3(scale3(dot3(u, u), cross3(v, w)), scale3(dot3(v, v), cross3(w, u))), scale3(dot3(w, w), cross3(u, v)));
    add3(a, scale3(0.5 / dot3(u, cross3(v, w)), num))
}
