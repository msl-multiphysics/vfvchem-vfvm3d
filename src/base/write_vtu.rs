use crate::base::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_scldom_vtu(scldom: &ScalarDomain, pts: &Points, nds: &Nodes, dom: &Domain, tet: &Tetrahedron, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", scldom.file_name, ts, scldom.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", dom.num_tot(), tet.tet_vid.len()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gpt in &dom.pts_gpt {
        let (x, y, z) = pts.pts_pos[gpt];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    for &gnd in &dom.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (Delaunay tetrahedra)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &[vid_0, vid_1, vid_2, vid_3] in &tet.tet_vid {
        writeln!(w, "{} {} {} {}", vid_0, vid_1, vid_2, vid_3).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for tet_idx in 0..tet.tet_vid.len() {
        writeln!(w, "{}", 4 * (tet_idx + 1)).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..tet.tet_vid.len() {
        writeln!(w, "10").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for vid in 0..dom.num_tot() {
        writeln!(w, "{:.9}", scldom.value[vid]).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_sclbnd_vtu(sclbnd: &ScalarBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclbnd.file_name, ts, sclbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", bnd.num_nds(), bnd.num_nds()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gnd in &bnd.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (degenerate to 1 point per cell)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        writeln!(w, "{}", nid).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        writeln!(w, "{}", nid + 1).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..bnd.num_nds() {
        writeln!(w, "1").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        writeln!(w, "{:.9}", sclbnd.value[nid]).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_sclitf_vtu(sclitf: &ScalarInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclitf.file_name, ts, sclitf.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", itf.num_nds(), itf.num_nds()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gnd in &itf.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (degenerate to 1 point per cell)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..itf.num_nds() {
        writeln!(w, "{}", nid).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..itf.num_nds() {
        writeln!(w, "{}", nid + 1).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..itf.num_nds() {
        writeln!(w, "1").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_vecdom_vtu(vecdom: &VectorDomain, pts: &Points, nds: &Nodes, dom: &Domain, tet: &Tetrahedron, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecdom.file_name, ts, vecdom.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", dom.num_tot(), tet.tet_vid.len()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gpt in &dom.pts_gpt {
        let (x, y, z) = pts.pts_pos[gpt];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    for &gnd in &dom.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (Delaunay tetrahedra)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &[vid_0, vid_1, vid_2, vid_3] in &tet.tet_vid {
        writeln!(w, "{} {} {} {}", vid_0, vid_1, vid_2, vid_3).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for tet_idx in 0..tet.tet_vid.len() {
        writeln!(w, "{}", 4 * (tet_idx + 1)).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..tet.tet_vid.len() {
        writeln!(w, "10").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for vid in 0..dom.num_tot() {
        let (vx, vy, vz) = vecdom.value[vid];
        writeln!(w, "{:.9} {:.9} {:.9}", vx, vy, vz).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_vecbnd_vtu(vecbnd: &VectorBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecbnd.file_name, ts, vecbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", bnd.num_nds(), bnd.num_nds()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gnd in &bnd.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (degenerate to 1 point per cell)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        writeln!(w, "{}", nid).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        writeln!(w, "{}", nid + 1).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..bnd.num_nds() {
        writeln!(w, "1").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..bnd.num_nds() {
        let (vx, vy, vz) = vecbnd.value[nid];
        writeln!(w, "{:.9} {:.9} {:.9}", vx, vy, vz).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_vecitf_vtu(vecitf: &VectorInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", vecitf.file_name, ts, vecitf.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", itf.num_nds(), itf.num_nds()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gnd in &itf.nds_gnd {
        let (x, y, z) = nds.nds_pos[gnd];
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (degenerate to 1 point per cell)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..itf.num_nds() {
        writeln!(w, "{}", nid).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..itf.num_nds() {
        writeln!(w, "{}", nid + 1).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..itf.num_nds() {
        writeln!(w, "1").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    
    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for nid in 0..itf.num_nds() {
        let (vx, vy, vz) = vecitf.value[nid];
        writeln!(w, "{:.9} {:.9} {:.9}", vx, vy, vz).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;    
    
    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

fn _vtu_err(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteVtu { file_path: file_path.to_string(), source }
}
