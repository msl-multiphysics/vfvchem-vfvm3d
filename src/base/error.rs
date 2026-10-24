use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum VFVChemError {
    // placeholder
    #[error("Placeholder")]
    Placeholder {},

    // vertices
    #[error("Points: Maximum sampling attempts `max_try` must be greater than zero. Got max_try = {max_try}.")]
    PointsNonPositiveAttempts {
        max_try: usize,
    },
    #[error("Points: Minimum spacing `h_min` must be positive and finite. Got h_min = {h_min}.")]
    PointsInvalidMinSpacing {
        h_min: f64,
    },
    #[error("Points: Maximum spacing `h_max` must be finite and at least h_min. Got h_min = {h_min} and h_max = {h_max}.")]
    PointsInvalidMaxSpacing {
        h_min: f64,
        h_max: f64,
    },
    #[error("Points: Refinement layer thickness `ref_len` must be positive and finite. Got ref_len = {ref_len}.")]
    PointsInvalidRefinementLength {
        ref_len: f64,
    },
    #[error("Nodes: x_max must be greater than x_min. Got x_min = {x_min} and x_max = {x_max}.")]
    NodesInvalidBoundsX {
        x_min: f64,
        x_max: f64,
    },
    #[error("Nodes: y_max must be greater than y_min. Got y_min = {y_min} and y_max = {y_max}.")]
    NodesInvalidBoundsY {
        y_min: f64,
        y_max: f64,
    },
    #[error("Nodes: z_max must be greater than z_min. Got z_min = {z_min} and z_max = {z_max}.")]
    NodesInvalidBoundsZ {
        z_min: f64,
        z_max: f64,
    },
    #[error("Nodes: Spacing `h` must be greater than zero. Got h = {h}.")]
    NodesNonPositiveSpacing {
        h: f64,
    },
    #[error("Nodes: Invalid geometry: {reason}")]
    NodesInvalidGeometry {
        reason: String,
    },

    // geometry
    #[error("Domain: Group ID {grp3d_id} not found.")]
    DomainGroupIDNotFound {
        grp3d_id: String,
    },
    #[error("Boundary: Group ID {grp2d_id} not found in Nodes.")]
    BoundaryGroupIDNotFound {
        grp2d_id: String,
    },
    #[error("Interface: Group ID {grp2d_id} not found in Nodes.")]
    InterfaceGroupIDNotFound {
        grp2d_id: String,
    },
    #[error("Interface: Domains {dom_a} and {dom_b} do not share an interface at group ID {grp2d_id}.")]
    InterfaceNoCommonInterface {
        dom_a: usize,
        dom_b: usize,
        grp2d_id: String,
    },

    // scalars and vectors
    #[error("Scalar{geom}: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    ScalarFunctionDependsOnNonUnknown {
        geom: &'static str,
        scldom_id: usize,
    },
    #[error("Scalar{geom}: Output file must be CSV or VTU. Found {file_type}.")]
    ScalarUnsupportedOutputType {
        geom: &'static str,
        file_type: String,
    },
    #[error("Scalar{geom}: Cannot deduce output file type from path: {file_path}.")]
    ScalarUnknownOutputType {
        geom: &'static str,
        file_path: String,
    },
    #[error("Vector{geom}: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    VectorFunctionDependsOnNonUnknown {
        geom: &'static str,
        scldom_id: usize,
    },
    #[error("Vector{geom}: Output file must be CSV or VTU. Found {file_type}.")]
    VectorUnsupportedOutputType {
        geom: &'static str,
        file_type: String,
    },
    #[error("Vector{geom}: Cannot deduce output file type from path: {file_path}.")]
    VectorUnknownOutputType {
        geom: &'static str,
        file_path: String,
    },

    // read and write
    #[error("ReadObj: Cannot open or read file '{file_path}': {source}")]
    ReadObjFile {
        file_path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("ReadObj: Parse error in '{file_path}': {reason}")]
    ReadObjParse {
        file_path: String,
        reason: String,
    },
    #[error("WriteCsv: Cannot write file '{file_path}': {source}")]
    WriteCsv {
        file_path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("WriteVtu: Cannot write file '{file_path}': {source}")]
    WriteVtu {
        file_path: String,
        #[source]
        source: std::io::Error,
    },

    // triangulation and tetrahedralization
    #[error("Face: Could not insert node into the triangulation of face {face_id}.")]
    FaceInsertionError {
        face_id: usize,
        #[source]
        source: spade::InsertionError,
    },
    #[error("Face: Rim edges of face {face_id} intersect each other.")]
    FaceConstraintIntersection {
        face_id: usize,
    },
    #[error("Tetrahedron: Vertices of domain {dom_id} do not span a volume.")]
    TetraDegenerateInput {
        dom_id: usize,
    },
    #[error("Tetrahedron: {num_bnd} interior points of domain {dom_id} lie on the boundary of the tetrahedralization. Try a smaller h_min or another seed.")]
    TetraPointsOnBoundary {
        dom_id: usize,
        num_bnd: usize,
    },
    #[error("Tetrahedron: Could not compute gradient weights at domain {dom_id} vertex {vid}; neighbors are coplanar.")]
    TetraSingularGradient {
        dom_id: usize,
        vid: usize,
    },

    // solvers
    #[error("Solver: Time step must be greater than zero. Got dt = {dt}.")]
    SolverNonPositiveTimeStep {
        dt: f64,
    },
    #[error("Solver: Number of time steps must be greater than zero. Got num_ts = {num_ts}.")]
    SolverNonPositiveNumTimeSteps {
        num_ts: usize,
    },
    #[error("Solver: Write frequency must be greater than zero. Got num_ts_write = {num_ts_write}.")]
    SolverNonPositiveWriteFrequency {
        num_ts_write: usize,
    },
    #[error("Solver: Maximum iterations must be greater than zero. Got max_iter = {max_iter}.")]
    SolverNonPositiveMaxIter {
        max_iter: usize,
    },
    #[error("Solver: Tolerance must be greater than zero. Got tol = {tol}.")]
    SolverNonPositiveTolerance {
        tol: f64,
    },
    #[error("Variables: Reference scale must be greater than zero. Got scale = {scale}.")]
    VariablesNonPositiveScale {
        scale: f64,
    },
    #[error("Solver: Damping must be in (0, 1]. Got damp = {damp}.")]
    SolverInvalidDamping {
        damp: f64,
    },
    #[error("SteadyBase: Nonlinear iteration did not converge within {max_iter} iterations. Residual = {res:.4e}.")]
    SteadyDidNotConverge {
        max_iter: usize,
        res: f64,
    },
    #[error("TransientBase: Nonlinear iteration did not converge at time step {ts} within {max_iter} iterations. Residual = {res:.4e}.")]
    TransientDidNotConverge {
        ts: usize,
        max_iter: usize,
        res: f64,
    },
    #[error("SteadyBase: Nonlinear iteration diverged; the solution is not finite (NaN or infinity) at iteration {iter}.")]
    SteadyNonFiniteResidual {
        iter: usize,
    },
    #[error("TransientBase: Nonlinear iteration diverged; the solution is not finite (NaN or infinity) at time step {ts}, iteration {iter}.")]
    TransientNonFiniteResidual {
        ts: usize,
        iter: usize,
    },
    #[error("Mass: Number of components must be greater than zero. Got num_comp = {num_comp}.")]
    MassNonPositiveNumComponents {
        num_comp: usize,
    },
    #[error("Mass: Component index must be less than the number of components. Got comp = {comp} and num_comp = {num_comp}.")]
    MassInvalidComponent {
        comp: usize,
        num_comp: usize,
    },
    #[error("Mass: Interface ID {itf_id} joins a domain with flow and a domain without flow. Concentrations are per unit mass on domains with flow and per unit volume otherwise.")]
    MassInterfaceMixedFlow {
        itf_id: usize,
    },

    // linear solvers
    #[error("LinearLu: Could not construct symbolic sparse matrix.")]
    LinearLuSymbolicError {
        #[source]
        source: faer::sparse::CreationError,
    },
    #[error("LinearLu: Could not assemble sparse matrix.")]
    LinearLuMatrixError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearLu: Could not compute symbolic LU factorization.")]
    LinearLuSymbolicLuError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearLu: Could not compute LU factorization.")]
    LinearLuLuError {
        #[source]
        source: faer::sparse::linalg::LuError,
    },
    #[error("LinearGmres: Number of threads must be greater than zero. Got num_thread = {num_thread}.")]
    LinearGmresNonPositiveThreads {
        num_thread: usize,
    },
    #[error("LinearGmres: Tolerance must be greater than zero. Got tol = {tol}.")]
    LinearGmresNonPositiveTolerance {
        tol: f64,
    },
    #[error("LinearGmres: Restart size must be greater than zero. Got restart = {restart}.")]
    LinearGmresNonPositiveRestart {
        restart: usize,
    },
    #[error("LinearGmres: Maximum outer iterations must be greater than zero. Got max_outer = {max_outer}.")]
    LinearGmresNonPositiveMaxOuter {
        max_outer: usize,
    },
    #[error("LinearGmres: Could not construct symbolic sparse matrix.")]
    LinearGmresSymbolicError {
        #[source]
        source: faer::sparse::CreationError,
    },
    #[error("LinearGmres: Could not assemble sparse matrix.")]
    LinearGmresMatrixError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearGmres: Solver did not converge. {message}")]
    LinearGmresDidNotConverge {
        message: String,
    },

    // operators
    #[error("Operator: Use apply_time instead of apply.")]
    OperatorUseApplyTime {},
    #[error("Operator: Only unknown-type scalars and vectors can be solved for in the global linear system.")]
    OperatorNonUnknownInLinearSystem {},
    #[error("Operator: Node ID {nid} for the reference pressure is not on Domain ID {dom_id}.")]
    OperatorNodeNotFound {
        dom_id: usize,
        nid: usize,
    },
    #[error("Operator: Node ID {nid} of Domain ID {dom_id} for the reference pressure is on Interface ID {itf_id}; choose a node off the interface.")]
    OperatorNodeOnInterface {
        dom_id: usize,
        nid: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarDomainNotOnDomain {
        scldom_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarDomainNotOnBoundary {
        scldom_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarBoundary ID {sclbnd_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarBoundaryNotOnBoundary {
        sclbnd_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarDomain IDs {scldom_a_id} and {scldom_b_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarDomainNotOnInterface {
        scldom_a_id: usize,
        scldom_b_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarInterface ID {sclitf_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarInterfaceNotOnInterface {
        sclitf_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarFlowNotOnDomain {
        sclflw_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarFlowNotOnBoundary {
        sclflw_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorDomain ID {vecdom_id} is not defined on Domain ID {dom_id}.")]
    OperatorVectorDomainNotOnDomain {
        vecdom_id: usize,
        dom_id: usize,
    },
    #[error("Operator: VectorDomain ID {vecdom_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorVectorDomainNotOnBoundary {
        vecdom_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorBoundary ID {vecbnd_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorVectorBoundaryNotOnBoundary {
        vecbnd_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: VectorDomain IDs {vecdom_a_id} and {vecdom_b_id} is not defined on Interface ID {itf_id}.")]
    OperatorVectorDomainNotOnInterface {
        vecdom_a_id: usize,
        vecdom_b_id: usize,
        itf_id: usize,
    },
    #[error("Operator: VectorInterface ID {vecitf_id} is not defined on Interface ID {itf_id}.")]
    OperatorVectorInterfaceNotOnInterface {
        vecitf_id: usize,
        itf_id: usize,
    },
}
