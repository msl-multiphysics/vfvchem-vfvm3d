use crate::base::prelude::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use robust::{orient3d, Coord3D};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::f64::consts::PI;

pub const INVALID: usize = usize::MAX;
const MAX_CELL: f64 = 4e6;

#[derive(Default)]
pub struct Points {
    // sampling parameters
    // spacing grows linearly from h_min at the boundary to h_max at a distance ref_len
    pub h_min: f64,
    pub h_max: f64,
    pub ref_len: f64,
    pub max_try: usize,
    pub seed: u64,

    // points
    // points are located inside boundaries
    pub pts_pos: Vec<(f64, f64, f64)>,  // gpt (global point id) -> (x, y, z)

    // 3d groups
    // uses the same grp3d_id in nodes
    pub grp3d_gpt: HashMap<String, Vec<usize>>,  // grp3d_id -> [gpt, ...]
}

pub struct SurfaceGrid {
    // closed triangulated surface sorted into a background grid
    // surface triangles are listed in every cell their bounding box touches
    pub lo: (f64, f64, f64),  // bounding box
    pub hi: (f64, f64, f64),
    pub cell: f64,  // cell size
    pub num: (usize, usize, usize),  // number of cells along x, y, and z
    pos: Vec<Coord3D<f64>>,  // surface vertex -> position
    tri: Vec<[usize; 3]>,  // triangle -> [surface vertex, ...]; counterclockwise seen from outside
    cell_tri: Vec<Vec<usize>>,  // cell -> [triangle, ...]
    tri_lo: Vec<(usize, usize, usize)>,  // triangle -> first cell of its bounding box
}

impl Points {
    pub fn new(nds: &Nodes, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Result<Points, VFVChemError> {
        // input checking
        if max_try == 0 { return Err(VFVChemError::PointsNonPositiveAttempts { max_try }) }
        if !(h_min > 0.0 && h_min.is_finite()) { return Err(VFVChemError::PointsInvalidMinSpacing { h_min }) }
        if !(h_max >= h_min && h_max.is_finite()) { return Err(VFVChemError::PointsInvalidMaxSpacing { h_min, h_max }) }
        if !(ref_len > 0.0 && ref_len.is_finite()) { return Err(VFVChemError::PointsInvalidRefinementLength { ref_len }) }

        // default struct
        let mut pts = Points::default();
        pts.h_min = h_min;
        pts.h_max = h_max;
        pts.ref_len = ref_len;
        pts.max_try = max_try;
        pts.seed = seed;

        // iterate through 3d groups
        // sorted so that point ids do not depend on hash order
        let mut grp3d_ids: Vec<&String> = nds.grp3d_gnd.keys().collect();
        grp3d_ids.sort();
        for grp3d_id in grp3d_ids {
            // create internal points
            let pts_pos = _point_sample(nds, grp3d_id, h_min, h_max, ref_len, max_try, seed);
            let mut pts_id = Vec::new();

            // register internal points
            for pos in pts_pos.iter() {
                let gpt = pts.pts_pos.len();
                pts.pts_pos.push(*pos);
                pts_id.push(gpt);
            }
            pts.grp3d_gpt.insert(grp3d_id.clone(), pts_id);
        }

        Ok(pts)
    }

    #[inline]
    pub fn num_pts(&self) -> usize {
        self.pts_pos.len()
    }
}

impl SurfaceGrid {
    pub fn new(pos: &[(f64, f64, f64)], tri: &[[usize; 3]], cell: f64) -> SurfaceGrid {
        // bounding box
        let mut lo = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &(x, y, z) in pos {
            lo = (lo.0.min(x), lo.1.min(y), lo.2.min(z));
            hi = (hi.0.max(x), hi.1.max(y), hi.2.max(z));
        }

        // cells; the cell size is raised if the grid would be too large
        let (sx, sy, sz) = (hi.0 - lo.0, hi.1 - lo.1, hi.2 - lo.2);
        let cell = if sx * sy * sz / cell.powi(3) > MAX_CELL { (sx * sy * sz / MAX_CELL).cbrt() } else { cell };
        let num = (((sx / cell).ceil() as usize).max(1), ((sy / cell).ceil() as usize).max(1), ((sz / cell).ceil() as usize).max(1));
        let mut surf = SurfaceGrid {
            lo, hi, cell, num,
            pos: pos.iter().map(|&(x, y, z)| Coord3D { x, y, z }).collect(),
            tri: tri.to_vec(),
            cell_tri: vec![Vec::new(); num.0 * num.1 * num.2],
            tri_lo: Vec::with_capacity(tri.len()),
        };

        // load surface triangles into grid
        for (idx, t) in tri.iter().enumerate() {
            let (pa, pb, pc) = (pos[t[0]], pos[t[1]], pos[t[2]]);
            let (i0, i1) = surf.span(pa.0.min(pb.0).min(pc.0), pa.0.max(pb.0).max(pc.0), 0);
            let (j0, j1) = surf.span(pa.1.min(pb.1).min(pc.1), pa.1.max(pb.1).max(pc.1), 1);
            let (k0, k1) = surf.span(pa.2.min(pb.2).min(pc.2), pa.2.max(pb.2).max(pc.2), 2);
            for k in k0..=k1 {
            for j in j0..=j1 {
            for i in i0..=i1 {
                surf.cell_tri[i + j * num.0 + k * num.0 * num.1].push(idx);
            }}}
            surf.tri_lo.push((i0, j0, k0));
        }

        surf
    }

    pub fn span(&self, a: f64, b: f64, axis: usize) -> (usize, usize) {
        // cells from coordinate a to b along one axis
        let (o, n) = [(self.lo.0, self.num.0), (self.lo.1, self.num.1), (self.lo.2, self.num.2)][axis];
        ((((a - o) / self.cell).max(0.0) as usize).min(n - 1), (((b - o) / self.cell).max(0.0) as usize).min(n - 1))
    }

    pub fn cells(&self, a: (f64, f64, f64), b: (f64, f64, f64)) -> impl Iterator<Item = usize> + use<> {
        // cells touching the box from a to b
        let (i0, i1) = self.span(a.0, b.0, 0);
        let (j0, j1) = self.span(a.1, b.1, 1);
        let (k0, k1) = self.span(a.2, b.2, 2);
        let (nx, ny) = (self.num.0, self.num.1);
        (k0..=k1).flat_map(move |k| (j0..=j1).flat_map(move |j| (i0..=i1).map(move |i| i + j * nx + k * nx * ny)))
    }

    pub fn is_inside(&self, p: (f64, f64, f64)) -> bool {
        // true if the segment from p to just past the nearest side of the bounding box crosses the surface an odd number of times
        // the segment is tilted slightly, and tilted again if it touches a triangle edge or vertex
        let dist = [p.0 - self.lo.0, self.hi.0 - p.0, p.1 - self.lo.1, self.hi.1 - p.1, p.2 - self.lo.2, self.hi.2 - p.2];
        if dist.iter().any(|&d| d < 0.0) { return false; }
        let side = (0..6).min_by(|&a, &b| dist[a].total_cmp(&dist[b])).unwrap();
        let len = dist[side] + self.cell;
        let s = if side % 2 == 0 { -1.0 } else { 1.0 };
        for (t0, t1) in [(0.0123, 0.0071), (-0.0217, 0.0137), (0.0093, -0.0251), (-0.0311, -0.0189)] {
            let dir = match side / 2 { 0 => (s, t0, t1), 1 => (t0, s, t1), _ => (t0, t1, s) };
            if let Some(num_cross) = self.cross(p, add3(p, scale3(len, dir))) { return num_cross % 2 == 1; }
        }
        false
    }

    fn cross(&self, p: (f64, f64, f64), q: (f64, f64, f64)) -> Option<usize> {
        // number of surface triangles crossed by the segment pq; None if pq touches a triangle edge or vertex
        // a triangle listed in several cells is only counted in the first cell shared by the triangle and the segment
        let (i0, i1) = self.span(p.0.min(q.0), p.0.max(q.0), 0);
        let (j0, j1) = self.span(p.1.min(q.1), p.1.max(q.1), 1);
        let (k0, k1) = self.span(p.2.min(q.2), p.2.max(q.2), 2);
        let (cp, cq) = (Coord3D { x: p.0, y: p.1, z: p.2 }, Coord3D { x: q.0, y: q.1, z: q.2 });
        let mut num_cross = 0;
        for k in k0..=k1 {
        for j in j0..=j1 {
        for i in i0..=i1 {
            for &idx in &self.cell_tri[i + j * self.num.0 + k * self.num.0 * self.num.1] {
                let (ti, tj, tk) = self.tri_lo[idx];
                if (i, j, k) != (ti.max(i0), tj.max(j0), tk.max(k0)) { continue; }
                let [a, b, c] = self.tri[idx].map(|v| self.pos[v]);

                // endpoints on the same side of the plane do not cross
                let o_p = orient3d(a, b, c, cp);
                let o_q = orient3d(a, b, c, cq);
                if o_p * o_q > 0.0 { continue; }

                // the line through pq passes inside the triangle if all three edges are seen with the same sign
                let s = [orient3d(cp, cq, a, b), orient3d(cp, cq, b, c), orient3d(cp, cq, c, a)];
                let num_pos = s.iter().filter(|&&v| v > 0.0).count();
                let num_neg = s.iter().filter(|&&v| v < 0.0).count();
                if num_pos > 0 && num_neg > 0 { continue; }
                if (num_pos == 3 || num_neg == 3) && o_p != 0.0 && o_q != 0.0 { num_cross += 1; }
                else { return None; }
            }
        }}}
        Some(num_cross)
    }
}

fn _point_sample(nds: &Nodes, grp3d_id: &str, h_min: f64, h_max: f64, ref_len: f64, max_try: usize, seed: u64) -> Vec<(f64, f64, f64)> {
    // pseudocode:
    // for one 3d group, sort the surface triangles into a background grid
    // use the distance to the nearest surface node to get the local spacing of each candidate
    // a candidate is inside if a segment from it to outside the bounding box crosses the surface an odd number of times
    // seed from each surface triangle, then use variable-radius bridson poisson-disk sampling to fill the region of the seed
    // any number of regions and cavities are allowed

    // step 1: get surface and grid

    // surface nodes and triangles
    // vid (vertex id) refers to the combined node and point ids
    let mut gnd_vid = HashMap::new();  // gnd -> vid
    let mut surf_pos = Vec::new();  // vid -> (x, y, z); nodes only
    for &gnd in &nds.grp3d_gnd[grp3d_id] {
        gnd_vid.insert(gnd, surf_pos.len());
        surf_pos.push(nds.nds_pos[gnd]);
    }
    let tri_vid: Vec<[usize; 3]> = nds.grp3d_tri[grp3d_id].iter().map(|t| t.map(|gnd| gnd_vid[&gnd])).collect();

    // cells are sized near the smallest spacing; each cell keeps a list of vertices
    let surf = SurfaceGrid::new(&surf_pos, &tri_vid, h_min.max(0.25 * h_max));
    let lo = surf.lo;
    let mut grid = vec![Vec::new(); surf.cell_tri.len()];  // cell -> [vid, ...]
    let cell_of = |p: (f64, f64, f64)| surf.cells(p, p).next().unwrap();

    // load node points into grid
    // nodes are on the surface, so their spacing is h_min
    let mut pts_pos = Vec::new();  // vid -> (x, y, z)
    let mut pts_h = Vec::new();  // vid -> local spacing
    let mut pts_act = Vec::new();  // point indices where neighbors are to be added
    for &pos in &surf_pos {
        grid[cell_of(pos)].push(pts_pos.len());
        pts_pos.push(pos);
        pts_h.push(h_min);
    }

    // step 2: distance field

    // each cell keeps the surface node nearest to its center, found by a dijkstra sweep from the cells holding nodes
    // the sweep stops beyond ref_len, where the spacing is h_max anyway
    let cell = surf.cell;
    let (nx, ny) = (surf.num.0, surf.num.1);
    let center = |c: usize| (lo.0 + ((c % nx) as f64 + 0.5) * cell, lo.1 + (((c / nx) % ny) as f64 + 0.5) * cell, lo.2 + ((c / (nx * ny)) as f64 + 0.5) * cell);
    let mut cell_seed = vec![INVALID; grid.len()];  // cell -> vid of nearest node
    if h_max > h_min {
        let mut cell_dist = vec![f64::INFINITY; grid.len()];
        let mut heap = BinaryHeap::new();
        for (vid, &pos) in surf_pos.iter().enumerate() {
            let c = cell_of(pos);
            let d = len3(sub3(center(c), pos));
            if d < cell_dist[c] {
                cell_dist[c] = d;
                cell_seed[c] = vid;
            }
        }
        for c in 0..grid.len() {
            if cell_seed[c] != INVALID { heap.push(Reverse((cell_dist[c].to_bits(), c))); }
        }
        while let Some(Reverse((bits, c))) = heap.pop() {
            // skip stale entries and cells beyond the refinement layer
            if f64::from_bits(bits) > cell_dist[c] || cell_dist[c] > ref_len + 2.0 * cell { continue; }
            let cen = center(c);
            for c_n in surf.cells(sub3(cen, (cell, cell, cell)), add3(cen, (cell, cell, cell))) {
                let d = len3(sub3(center(c_n), surf_pos[cell_seed[c]]));
                if d < cell_dist[c_n] {
                    cell_dist[c_n] = d;
                    cell_seed[c_n] = cell_seed[c];
                    heap.push(Reverse((d.to_bits(), c_n)));
                }
            }
        }
    }

    // local spacing h(d) = h_min + (h_max - h_min) * min(d / ref_len, 1)
    // d is the distance to the nearest node among the seeds of the surrounding cells
    let space = |pos: (f64, f64, f64)| -> f64 {
        if h_max == h_min { return h_min; }
        let mut d = f64::INFINITY;
        for c in surf.cells(sub3(pos, (cell, cell, cell)), add3(pos, (cell, cell, cell))) {
            if cell_seed[c] != INVALID { d = d.min(len3(sub3(pos, surf_pos[cell_seed[c]]))); }
        }
        h_min + (h_max - h_min) * (d / ref_len).min(1.0)
    };

    // two points are far enough apart if their distance is at least the mean of their local spacings
    // existing spacings are at most h_max, so only cells within (hc + h_max) / 2 of the candidate are checked
    let is_far = |pos: (f64, f64, f64), hc: f64, pts_pos: &Vec<(f64, f64, f64)>, pts_h: &Vec<f64>, grid: &Vec<Vec<usize>>| -> bool {
        let r = 0.5 * (hc + h_max);
        for c in surf.cells(sub3(pos, (r, r, r)), add3(pos, (r, r, r))) {
            for &vid_n in &grid[c] {
                let d = sub3(pos, pts_pos[vid_n]);
                let hs = 0.5 * (hc + pts_h[vid_n]);
                if dot3(d, d) < hs * hs { return false; }
            }
        }
        true
    };

    // protecting balls
    // a point inside the smallest sphere through a surface triangle would move the circumcenters of the tetrahedra on that triangle outside the surface
    // each ball is centered on the circumcenter of its triangle; balls of very obtuse triangles are skipped to keep the region filled
    let mut ball = Vec::with_capacity(tri_vid.len());  // ball -> (center, radius)
    let mut ball_grid = vec![Vec::new(); grid.len()];  // cell -> [ball, ...]
    for tri in &tri_vid {
        let (cen, rad) = circumcircle3(surf_pos[tri[0]], surf_pos[tri[1]], surf_pos[tri[2]]);
        if !(rad <= 2.0 * h_max) { continue; }
        for c in surf.cells(sub3(cen, (rad, rad, rad)), add3(cen, (rad, rad, rad))) {
            ball_grid[c].push(ball.len());
        }
        ball.push((cen, rad));
    }
    let is_free = |pos: (f64, f64, f64)| -> bool {
        ball_grid[cell_of(pos)].iter().all(|&b| {
            let (cen, rad) = ball[b];
            let d = sub3(pos, cen);
            dot3(d, d) >= rad * rad
        })
    };

    // step 3: insert points

    // seed rng
    let mut rng = StdRng::seed_from_u64(seed);

    // growth only reaches twice the local spacing from existing points, so separate regions need their own starting point
    // candidates are a distance h_min or h_max inside the centroid of each surface triangle
    // the h_max offset always clears the surface nodes, while the h_min offset fits in thin regions
    for tri in &tri_vid {
        // get inward offsets of the centroid
        let (pa, pb, pc) = (surf_pos[tri[0]], surf_pos[tri[1]], surf_pos[tri[2]]);
        let norm = cross3(sub3(pb, pa), sub3(pc, pa));
        let len = len3(norm);
        if len == 0.0 { continue; }
        let cent = scale3(1.0 / 3.0, add3(pa, add3(pb, pc)));

        // register the first valid candidate
        let mut is_seed = false;
        for h in [h_min, h_max] {
            // check if point is (1) far enough from other points at its local spacing, (2) outside the protecting balls, and (3) in the domain
            let pos = sub3(cent, scale3(h / len, norm));
            let hc = space(pos);
            if !is_far(pos, hc, &pts_pos, &pts_h, &grid) || !is_free(pos) || !surf.is_inside(pos) { continue; }

            // register initial point
            grid[cell_of(pos)].push(pts_pos.len());
            pts_act.push(pts_pos.len());
            pts_pos.push(pos);
            pts_h.push(hc);
            is_seed = true;
            break;
        }
        if !is_seed { continue; }

        // insert points in the region of the seed
        while !pts_act.is_empty() {
            // choose random active point
            let a = rng.random_range(0..pts_act.len());
            let vid_a = pts_act[a];
            let pos_a = pts_pos[vid_a];
            let ha = pts_h[vid_a];

            // try to insert point next to it
            let mut is_insert = false;  // true if successfully inserted
            for _ in 0..max_try {
                // candidate point a distance [ha, 2ha] from active point in a uniformly random direction
                let cos_pol = rng.random_range(-1.0..1.0);
                let sin_pol = f64::sqrt(1.0 - cos_pol * cos_pol);
                let ang = rng.random_range(0.0..(2.0 * PI));
                let rad = rng.random_range(ha..(2.0 * ha));
                let pos = add3(pos_a, scale3(rad, (sin_pol * ang.cos(), sin_pol * ang.sin(), cos_pol)));

                // check if point is (1) far enough from other points at its local spacing, (2) outside the protecting balls, and (3) in the domain
                let hc = space(pos);
                if !is_far(pos, hc, &pts_pos, &pts_h, &grid) || !is_free(pos) || !surf.is_inside(pos) { continue; }

                // register candidate point
                grid[cell_of(pos)].push(pts_pos.len());
                pts_act.push(pts_pos.len());
                pts_pos.push(pos);
                pts_h.push(hc);
                is_insert = true;
                break;
            }

            // flag active point as inactive if candidates cannot be inserted
            if !is_insert {
                pts_act.swap_remove(a);
            }
        }
    }

    // return internal points
    // skip first num_nds points
    pts_pos.split_off(surf_pos.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    const H_MIN: f64 = 0.05;
    const H_MAX: f64 = 0.15;
    const REF_LEN: f64 = 0.3;

    fn _cube(seed: u64) -> (Nodes, Points) {
        let nds = Nodes::new_from_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, H_MIN).unwrap();
        let pts = Points::new(&nds, H_MIN, H_MAX, REF_LEN, 30, seed).unwrap();
        (nds, pts)
    }

    #[test]
    fn invalid_parameters_are_rejected() {
        let nds = Nodes::new_from_bounds(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, H_MIN).unwrap();
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, REF_LEN, 0, 42), Err(VFVChemError::PointsNonPositiveAttempts { .. })));
        assert!(matches!(Points::new(&nds, 0.0, H_MAX, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMinSpacing { .. })));
        assert!(matches!(Points::new(&nds, f64::NAN, H_MAX, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMinSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, 0.5 * H_MIN, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMaxSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, f64::INFINITY, REF_LEN, 30, 42), Err(VFVChemError::PointsInvalidMaxSpacing { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, 0.0, 30, 42), Err(VFVChemError::PointsInvalidRefinementLength { .. })));
        assert!(matches!(Points::new(&nds, H_MIN, H_MAX, f64::NAN, 30, 42), Err(VFVChemError::PointsInvalidRefinementLength { .. })));
    }

    #[test]
    fn sampling_is_deterministic() {
        let (_, pts_a) = _cube(7);
        let (_, pts_b) = _cube(7);
        assert!(pts_a.num_pts() > 0);
        assert_eq!(pts_a.pts_pos, pts_b.pts_pos);
    }

    #[test]
    fn points_respect_min_spacing() {
        let (nds, pts) = _cube(42);
        for (p, &pos_p) in pts.pts_pos.iter().enumerate() {
            assert!(pos_p.0 > 0.0 && pos_p.1 > 0.0 && pos_p.2 > 0.0 && pos_p.0 < 1.0 && pos_p.1 < 1.0 && pos_p.2 < 1.0);
            for (q, &pos_q) in pts.pts_pos.iter().enumerate() {
                if q == p { continue; }
                assert!(len3(sub3(pos_p, pos_q)) >= H_MIN - 1e-12);
            }
            for &pos_n in &nds.nds_pos {
                assert!(len3(sub3(pos_p, pos_n)) >= H_MIN - 1e-12);
            }
        }
    }

    #[test]
    fn points_are_denser_near_boundary() {
        // band within 0.1 of the boundary has volume 0.488; center cube farther than 0.3 has volume 0.064
        let (_, pts) = _cube(42);
        let dist = |&(x, y, z): &(f64, f64, f64)| x.min(1.0 - x).min(y).min(1.0 - y).min(z).min(1.0 - z);
        let num_near = pts.pts_pos.iter().filter(|p| dist(p) < 0.1).count() as f64;
        let num_far = pts.pts_pos.iter().filter(|p| dist(p) > 0.3).count() as f64;
        assert!(num_near / 0.488 > 3.0 * num_far / 0.064);
    }
}
