template<typename T, typename M>
inline bool cit::strictlyEqual(const Element<T, M>& e1, const Element<T, M>& e2) {
  return e1.getMesh() == e2.getMesh() && e1.getIndex() == e2.getIndex();
}

inline bool cit::isEdgeCompatible(const CoInTri& cointri, Edge e) {
  return cointri.correspondingEdge[e] != Edge();
}

inline bool cit::isFaceCompatible(const CoInTri& cointri, Face f) {
  return cointri.correspondingFace[f] != Face();
}

inline bool cit::isFaceReversed(const CoInTri& cointri, Face f) {
  return cointri.reversedFace[f] != Face();
}

template <typename T>
inline Vector3 cit::getRandomColor(const T& obj) {
  std::hash<T> h;
  std::mt19937 gen(h(obj) + (size_t)cit::tweaks.i["seed"]);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  return { dist(gen), dist(gen), dist(gen) };
}

template <typename Vector2T>
typename Vector2T::Scalar cit::computeTriangleEnergy(const std::array<Vector2T, 3>& A_triangle, const std::array<Vector2T, 3>& B_triangle, typename Vector2T::Scalar* A_area_ptr, typename Vector2T::Scalar* B_area_ptr) {
  using Scalar = typename Vector2T::Scalar;
  using Matrix2T = Eigen::Matrix<Scalar, 2, 2>;

  Vector2T A_d0 = A_triangle[1] - A_triangle[0];
  Vector2T A_d1 = A_triangle[2] - A_triangle[0];
  Vector2T B_d0 = B_triangle[1] - B_triangle[0];
  Vector2T B_d1 = B_triangle[2] - B_triangle[0];

  Scalar A_area = cross(A_d0, A_d1) / 2;
  Scalar B_area = cross(B_d0, B_d1) / 2;
  if (A_area_ptr) *A_area_ptr = A_area;
  if (B_area_ptr) *B_area_ptr = B_area;

  Matrix2T A_M, A_Minv;
  Matrix2T B_M, B_Minv;
  A_M << A_d0, A_d1;
  B_M << B_d0, B_d1;

  A_Minv <<
    A_d1[1], -A_d1[0],
    -A_d0[1], A_d0[0];
  B_Minv <<
    B_d1[1], -B_d1[0],
    -B_d0[1], B_d0[0];
  A_Minv /= 2 * A_area;
  B_Minv /= 2 * B_area;

  Matrix2T AB_J = B_M * A_Minv;
  Matrix2T BA_J = A_M * B_Minv;

  Scalar res{0};
  res += B_area * AB_J.squaredNorm();
  res += A_area * BA_J.squaredNorm();
  return res;
}

template <typename Vector2T>
inline typename Vector2T::Scalar cit::cross(const Vector2T& u, const Vector2T& v) { return u[0] * v[1] - u[1] * v[0]; }

inline double cit::angleAbsDifference(double angle1, double angle2) {
  while (angle1 < 0) angle1 += 2. * M_PI;
  while (angle2 < 0) angle2 += 2. * M_PI;
  angle1 = std::fmod(angle1, 2. * M_PI);
  angle2 = std::fmod(angle2, 2. * M_PI);
  if (angle1 > angle2) std::swap(angle1, angle2);
  return std::min<double>(angle2 - angle1, angle1 + 2. * M_PI - angle2);
}

inline double cit::standardizeAngle(double angle, double modulo) {
  while (angle < 0)
    angle += modulo;
  return std::fmod(angle, modulo);
}

inline Matrix6d cit::getHessian(const ADScalar_6v_2nd& y) {
  Matrix6d res;
  for (int i = 0; i < 6; ++i)
    res.col(i) = y.derivatives()[i].derivatives();
  return res;
}

template <typename MatrixT>
inline MatrixT cit::projectHessian(const MatrixT& H) {
  Eigen::SelfAdjointEigenSolver<MatrixT> eigensolver(H);
  CIT_ASSERT(eigensolver.info() == Eigen::Success);
  auto eigenvalues = eigensolver.eigenvalues();
  auto eigenvectors = eigensolver.eigenvectors();
  for (int i = 0; i < MatrixT::RowsAtCompileTime; ++i) {
    if (eigenvalues[i] < 0)
      eigenvalues[i] = sysParam.hessianProjectionValue;
  }
  return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

inline std::array<Vertex, 2> cit::getOppositeVertices(Edge e) {
  return {
    e.halfedge().next().next().vertex(),
    e.halfedge().twin().next().next().vertex(),
  };
}

inline SurfacePoint cit::convertEdgePointToFacePoint(SurfacePoint edgePoint) {
  CIT_ASSERT(edgePoint.type == SurfacePointType::Edge);
  return edgePoint.inFace(edgePoint.edge.halfedge().face());
}

inline SurfacePoint cit::getCorrespondingEdgePoint(const CoInTri& cointriP, SurfacePoint P_intrinsicSP) {
  CIT_ASSERT(P_intrinsicSP.type == SurfacePointType::Edge);
  Edge P_intrinsicE = P_intrinsicSP.edge;
  Edge Q_intrinsicE = cointriP.correspondingEdge[P_intrinsicE];
  CIT_ASSERT(cointriP.correspondingVertex[P_intrinsicE.halfedge().vertex()] == Q_intrinsicE.halfedge().vertex());
  CIT_ASSERT(cointriP.correspondingVertex[P_intrinsicE.halfedge().tipVertex()] == Q_intrinsicE.halfedge().tipVertex());
  return {Q_intrinsicE, P_intrinsicSP.tEdge};
}

template <typename T>
inline T* cit::staticValuePtr(T value) {
  static std::map<T, T> m;
  m[value] = value;
  return &m[value];
}

inline bool cit::isVertexMerged(int uniqueVertexID) {
  const size_t A_nV = mdataA.nV;
  const size_t B_nV = mdataB.nV;
  return uniqueVertexID >= (int)(A_nV < B_nV ? A_nV : B_nV) + 2;
}

inline bool cit::isVertexMerged(const CoInTri& cointri, Vertex v) {
  return isVertexMerged(cointri.uniqueID_per_Vertex[v]);
}

inline bool cit::isVertexMovable(const CoInTri& cointri, Vertex v) {
  return v.getIndex() >= cointri.mdata->nV && !isVertexMerged(cointri, v);
}

inline double cit::getRadiusRatio(const SignpostIntrinsicTriangulation& geometry, Face f) {
  // Philippe Pébay and Timothy Baker. 2003. Analysis of triangle quality measures. Math. Comp. 72, 244 (2003), 1817--1839
  double sineAlpha = std::sin(geometry.cornerAngle(f.halfedge().corner()));
  double sineBeta  = std::sin(geometry.cornerAngle(f.halfedge().next().corner()));
  double sineGamma = std::sin(geometry.cornerAngle(f.halfedge().next().next().corner()));
  return (sineAlpha + sineBeta + sineGamma) / (2. * sineAlpha * sineBeta * sineGamma);
}

inline std::string cit::getResultFileNameBase() {
  std::string res = "results/" + caseName;
  if (sysParam.fixAnchors)
    res += "-fixAnchors";
  res += "/" + std::to_string(sessionID);
  return res;
}

inline std::string cit::getResultFileName(ResultContext context, uint64_t id, uint64_t previous_id, int l) {
  if (context == ResultContext::UpdateConfiguration)
    return getResultFileNameBase() + "/" + std::to_string(id);

  return getResultFileNameBase() + "/" + std::to_string(previous_id) +
    "/" + (context == ResultContext::FindMaxStep ? "findMaxStep" : "findOptimalStep") +
    "/l=" + std::to_string(l) + "/" + std::to_string(id);
}
