use crate::base::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub struct ObjGroup {
    pub id: String,
    pub is_dom: bool,
    pub faces: Vec<Vec<usize>>,  // face -> [vertex index, ...]; 0-based, in winding order
}

pub fn read_nds_obj(file_path: &str, h: f64, offset_ratio: f64) -> Result<Nodes, VFVChemError> {
    // read boundary, interface, and domain polygons from an obj file
    // every f record must belong to a group named bnd*, itf*, or dom*
    // each f record is one planar polygon; the winding sets its normal by the right-hand rule
    // the faces of a domain group form closed surfaces with every normal pointing out of the domain
    // a face shared by groups (same vertex indices) reuses the same node ids
    // nodes are placed by Nodes::new_from_faces

    // file reader
    let file = File::open(file_path).map_err(|e| _obj_file(file_path, e))?;
    let reader = BufReader::new(file);

    // obj contents
    let mut verts = Vec::new();
    let mut groups: Vec<ObjGroup> = Vec::new();
    let mut cur: Option<usize> = None;

    // iterate over lines
    for (line_idx, line) in reader.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line.map_err(|e| _obj_file(file_path, e))?;
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() { continue; }

        let mut tok = line.split_whitespace();
        let tag = tok.next().unwrap();
        match tag {
            "v" => {
                let coords: Vec<&str> = tok.collect();
                if !(3..=4).contains(&coords.len()) {
                    return Err(_obj_parse(file_path, format!("line {line_no}: vertex must have x, y, and z coordinates")));
                }
                let x = _obj_f64(file_path, line_no, coords[0])?;
                let y = _obj_f64(file_path, line_no, coords[1])?;
                let z = _obj_f64(file_path, line_no, coords[2])?;
                verts.push((x, y, z));
            }
            "g" => {
                let Some(name) = tok.next() else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group must have one name")));
                };
                if tok.next().is_some() {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group must have one name")));
                }
                let is_dom = if name.starts_with("dom") {
                    true
                } else if name.starts_with("bnd") || name.starts_with("itf") {
                    false
                } else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: group name '{name}' must begin with 'bnd', 'itf', or 'dom'")));
                };

                // repeated declarations continue the same group
                if let Some(idx) = groups.iter().position(|grp| grp.id == name) {
                    cur = Some(idx);
                } else {
                    groups.push(ObjGroup { id: name.to_string(), is_dom, faces: Vec::new() });
                    cur = Some(groups.len() - 1);
                }
            }
            "f" => {
                let Some(idx) = cur else {
                    return Err(_obj_parse(file_path, format!("line {line_no}: face is not in a group")));
                };
                let mut vids = Vec::new();
                for tok in tok {
                    // keep the vertex index of v/vt/vn records
                    let tok = tok.split('/').next().unwrap_or("");
                    let vid: usize = tok.parse().map_err(|_| {
                        _obj_parse(file_path, format!("line {line_no}: invalid vertex index '{tok}'"))
                    })?;
                    if vid == 0 || vid > verts.len() {
                        return Err(_obj_parse(file_path, format!("line {line_no}: vertex index {vid} is out of range")));
                    }
                    vids.push(vid - 1);
                }
                if vids.len() < 3 {
                    return Err(_obj_parse(file_path, format!("line {line_no}: face must contain at least three vertices")));
                }
                groups[idx].faces.push(vids);
            }
            _ => {
                return Err(_obj_parse(file_path, format!("line {line_no}: unsupported record '{tag}'")));
            }
        }
    }

    // place nodes; geometry errors are reported against the file
    Nodes::new_from_faces(&verts, &groups, h, offset_ratio).map_err(|e| match e {
        VFVChemError::NodesInvalidGeometry { reason } => _obj_parse(file_path, reason),
        e => e,
    })
}

fn _obj_f64(file_path: &str, line_no: usize, tok: &str) -> Result<f64, VFVChemError> {
    let val: f64 = tok.parse().map_err(|_| {
        _obj_parse(file_path, format!("line {line_no}: invalid number '{tok}'"))
    })?;
    if !val.is_finite() {
        return Err(_obj_parse(file_path, format!("line {line_no}: non-finite number '{tok}'")));
    }
    Ok(val)
}

fn _obj_file(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::ReadObjFile { file_path: file_path.to_string(), source }
}

fn _obj_parse(file_path: &str, reason: String) -> VFVChemError {
    VFVChemError::ReadObjParse { file_path: file_path.to_string(), reason }
}
