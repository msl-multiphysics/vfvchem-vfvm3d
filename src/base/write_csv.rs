use crate::base::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_scldom_csv(scldom: &ScalarDomain, pts: &Points, nds: &Nodes, dom: &Domain, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", scldom.file_name, ts, scldom.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (pid, &gpt) in dom.pts_gpt.iter().enumerate() {
        let (x, y, z) = pts.pts_pos[gpt];
        let value = scldom.value[pid];  // vid = pid for points
        writeln!(w, "{:.9},{:.9},{:.9},{:.9}", x, y, z, value).map_err(|e| _csv_err(&file_path, e))?;
    }
    for (nid, &gnd) in dom.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let vid= dom.num_pts() + nid;  // vid = num_pts + nid for nodes
        let value = scldom.value[vid];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9}", x, y, z, value).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_sclbnd_csv(sclbnd: &ScalarBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclbnd.file_name, ts, sclbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in bnd.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let value = sclbnd.value[idx];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9}", x, y, z, value).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_sclitf_csv(sclitf: &ScalarInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclitf.file_name, ts, sclitf.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in itf.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let value = sclitf.value[idx];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9}", x, y, z, value).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_vecdom_csv(vecdom: &VectorDomain, pts: &Points, nds: &Nodes, dom: &Domain, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecdom.file_name, ts, vecdom.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value_x,value_y,value_z").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (pid, &gpt) in dom.pts_gpt.iter().enumerate() {
        let (x, y, z) = pts.pts_pos[gpt];
        let (value_x, value_y, value_z) = vecdom.value[pid];  // vid = pid for points
        writeln!(w, "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}", x, y, z, value_x, value_y, value_z).map_err(|e| _csv_err(&file_path, e))?;
    }
    for (nid, &gnd) in dom.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let vid= dom.num_pts() + nid;  // vid = num_pts + nid for nodes
        let (value_x, value_y, value_z) = vecdom.value[vid];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}", x, y, z, value_x, value_y, value_z).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_vecbnd_csv(vecbnd: &VectorBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecbnd.file_name, ts, vecbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value_x,value_y,value_z").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in bnd.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let (value_x, value_y, value_z) = vecbnd.value[idx];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}", x, y, z, value_x, value_y, value_z).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_vecitf_csv(vecitf: &VectorInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecitf.file_name, ts, vecitf.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,y,z,value_x,value_y,value_z").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in itf.nds_gnd.iter().enumerate() {
        let (x, y, z) = nds.nds_pos[gnd];
        let (value_x, value_y, value_z) = vecitf.value[idx];  
        writeln!(w, "{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}", x, y, z, value_x, value_y, value_z).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

fn _csv_err(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteCsv { file_path: file_path.to_string(), source }
}
