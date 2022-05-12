#pragma once

#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/combining_hash_functions.h"

#include "kt84/graphics/DisplayList.hh"
#include "kt84/graphics/TextureObjectT.hh"
#include "kt84/geometry/CameraFree.hh"
#include "kt84/MaxMinSelector.hh"
#include "kt84/Timer.hh"
#include "kt84/ScopeExit.hh"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <filesystem>
namespace fs = std::filesystem;

namespace spdlog {
  using logger_ptr = std::shared_ptr<spdlog::logger>;   // so commonly used, why not already defined?
}

// Logging with given logger
namespace cit { namespace detail {
template<typename FormatString, typename... Args> void logTrace   (spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::trace   , std::string(2 * depth, ' ') + fmt, args...); }
template<typename FormatString, typename... Args> void logDebug   (spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::debug   , std::string(2 * depth, ' ') + fmt, args...); }
template<typename FormatString, typename... Args> void logInfo    (spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::info    , std::string(2 * depth, ' ') + fmt, args...); }
template<typename FormatString, typename... Args> void logWarn    (spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::warn    , std::string(2 * depth, ' ') + fmt, args...); }
template<typename FormatString, typename... Args> void logError   (spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::err     , std::string(2 * depth, ' ') + fmt, args...); }
template<typename FormatString, typename... Args> void logCritical(spdlog::logger_ptr logger, size_t depth, spdlog::source_loc loc, const FormatString &fmt, const Args &...args) { logger->log(loc, spdlog::level::critical, std::string(2 * depth, ' ') + fmt, args...); }
} }
#define LOG_TRACE(   logger, depth, ...) cit::detail::logTrace   (logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
#define LOG_DEBUG(   logger, depth, ...) cit::detail::logDebug   (logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
#define LOG_INFO(    logger, depth, ...) cit::detail::logInfo    (logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
#define LOG_WARN(    logger, depth, ...) cit::detail::logWarn    (logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
#define LOG_ERROR(   logger, depth, ...) cit::detail::logError   (logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
#define LOG_CRITICAL(logger, depth, ...) cit::detail::logCritical(logger, depth, spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, __VA_ARGS__)
// Logging with default logger
#define DLOG_TRACE(   depth, ...) LOG_TRACE   (spdlog::default_logger(), depth, __VA_ARGS__)
#define DLOG_DEBUG(   depth, ...) LOG_DEBUG   (spdlog::default_logger(), depth, __VA_ARGS__)
#define DLOG_INFO(    depth, ...) LOG_INFO    (spdlog::default_logger(), depth, __VA_ARGS__)
#define DLOG_WARN(    depth, ...) LOG_WARN    (spdlog::default_logger(), depth, __VA_ARGS__)
#define DLOG_ERROR(   depth, ...) LOG_ERROR   (spdlog::default_logger(), depth, __VA_ARGS__)
#define DLOG_CRITICAL(depth, ...) LOG_CRITICAL(spdlog::default_logger(), depth, __VA_ARGS__)


#define CIT_ASSERT(CONDITION) do { GC_SAFETY_ASSERT((CONDITION), ""); } while (false)

#define MAX_NTHREADS 8

#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <deque>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using kt84::Camera;
using kt84::MaxMinAverage;
using Eigen::Matrix2d;
using Eigen::Matrix4d;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector4f;
using Eigen::VectorXd;
using Eigen::AlignedBox3d;
using SparseMatrixd = Eigen::SparseMatrix<double>;

#include <unsupported/Eigen/AutoDiff>
// typenames for auto differentiating twice with 6 variables
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using ADScalar_6v_1st = Eigen::AutoDiffScalar<Vector6d>;
using ADScalar_6v_2nd = Eigen::AutoDiffScalar<Eigen::Matrix<ADScalar_6v_1st, 6, 1>>;
using ADVector2d_6v_2nd = Eigen::Matrix<ADScalar_6v_2nd, 2, 1>;
using ADVector3d_6v_2nd = Eigen::Matrix<ADScalar_6v_2nd, 3, 1>;
using ADMatrix2d_6v_2nd = Eigen::Matrix<ADScalar_6v_2nd, 2, 2>;

using VectorXsp = Eigen::Matrix<SurfacePoint, -1, 1>;

struct GLFWwindow;
struct ImGuiIO;

#define MAKE_FORMATTABLE(T) \
  template <> \
  struct fmt::formatter<T> : public formatter<std::string> { \
    template <typename FormatContext> \
    auto format(const T &value, FormatContext &ctx) -> decltype(ctx.out()) { \
      std::ostringstream oss; \
      oss << value; \
      return formatter<std::string>::format(oss.str(), ctx); \
    } \
  }

MAKE_FORMATTABLE(Vector2);
MAKE_FORMATTABLE(Vector3);
MAKE_FORMATTABLE(Vertex);
MAKE_FORMATTABLE(Edge);
MAKE_FORMATTABLE(Halfedge);
MAKE_FORMATTABLE(Face);
MAKE_FORMATTABLE(SurfacePoint);
MAKE_FORMATTABLE(MaxMinAverage);

template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct fmt::formatter<Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>> : public formatter<std::string> {
  template <typename FormatContext>
  auto format(const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols> &value, FormatContext &ctx) -> decltype(ctx.out()) {
    std::ostringstream oss;
    oss << value;
    return formatter<std::string>::format(oss.str(), ctx);
  }
};


namespace cit {

//---------+
// structs |
//---------+

struct Tweaks {
  std::map<std::string, int> i;
  std::map<std::string, bool> b;
  std::map<std::string, float> f;
  std::map<std::string, double> d;
  std::map<std::string, Vector3f> vec3;
  std::map<std::string, std::string> s;
};
extern Tweaks tweaks;

struct LightParam {
  bool enabled = false;
  float ambient = 0.f;
  float diffuse = 0.f;
  float specular = 0.f;
  Vector4f position = {0.f, 0.f, 1.f, 0.f};
};
extern std::array<LightParam, 8> lightParam;

struct MaterialParam {
  float ambient = 0.2f;
  float specular = 0.5f;
  float shininess = 20.f;
};
extern MaterialParam materialParam;

struct SystemParameter {
  // double angleThreshold = 0.05;
  double maxSignpostError = 1.e-5;

  // Face point whose smallest (resp. largest) barycentric coordinate is below (resp. above)
  // snapThreshold (resp. 1-snapThreshold) gets snapped to the closest edge (vertex)
  double snapThreshold = 0.005;

  size_t nThreads = 1;

  bool fixAnchors = false;

  double hessianProjectionValue = 1.0e-10;

  double smax_scaleFactor = 1.2;
  double s_scaleFactor = 0.8;
  double armijo = 1.e-6;
  size_t nMaxIter = 1000;

  double wL_initial = 1.e+9;
  double alpha = 1.2;

  size_t historyMaxSize = 5;

  bool firstOrderMode = false;

  double energyDeltaThreshold = 1.e-5;      // Optimization stops when energy improvement gets lower than this
};
extern SystemParameter sysParam;

// Connected set of incompatbile faces (used by convexifyConcaveIncompatiblePatch & mergeInconsistentVertexPairs)
struct Patch {
  std::set<Face>   faces;
  std::set<Edge>   edges;
  std::set<Vertex> vertices;
  std::set<Edge>   interiorEdges;
  std::set<Vertex> interiorVertices;
  std::set<Edge>        boundaryEdges;            // In no particular order
  std::set<Vertex>      boundaryVertices;         // Can be smaller than boundaryEdges when some vertices have non-half-disk topology
  std::vector<std::vector<Halfedge>> boundaryHalfedges;        // Ordered counter-clockwise; Can be more than one loops in case of holes
  std::map<Vertex, double> vertexAngleSums;       // Normalized (i.e. divided by intrinsicVertexAngleSums[v])
};
using PatchPair = std::pair<Patch, Patch>;

struct ModelData {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::unique_ptr<CornerData<Vector2>> texCoords;

  size_t nV;
  std::set<size_t> anchorVertexIDs;                // Input vertex indices which are "fixed" (vertex-vertex correspondence in the original HOT data to be kept intact)
  EdgeData<Matrix2d> transportOperator;

  AlignedBox3d bbox;
  std::string shortName;    // "A" or "B"
  std::string fullName;
  kt84::TextureObject texture;
  kt84::CameraFree camera;
  kt84::DisplayList dispList;
  Vector3f inputEdgeColor;
};
extern ModelData mdataA, mdataB;

struct OverlayVertex {
  /*
    Classification of overlay vertex and its corresponding non-null fields
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | Category                              | vA | vB | eA | eB | eIntrinsic | tEdgeA | tEdgeB |
    |=======================================|====|====|====|====|============|========|========|
    | A's vertex                            | *  | -  | -  | -  |     -      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | B's vertex                            | -  | *  | -  | -  |     -      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | Anchor vertex                         | *  | *  | -  | -  |     -      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | Crossing of A's edge & B's edge       | -  | -  | *  | *  |     -      |    *   |    *   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | Crossing of A's edge & intrinsic edge | -  | -  | *  | -  |     *      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | Crossing of B's edge & intrinsic edge | -  | -  | -  | *  |     *      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | A's vertex mapped onto B's edge       | *  | -  | -  | *  |     -      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
    | B's vertex mapped onto A's edge       | -  | *  | *  | -  |     -      |    -   |    -   |
    +---------------------------------------+----+----+----+----+------------+--------+--------+
  */
  Vertex vA;
  Vertex vB;
  Edge eA;
  Edge eB;

  // Because intrinsic meshes for both A and B are compatible, this is always pointed to A's mesh
  // (same idea for OverlayWedge)
  Edge eIntrinsic;

  // These are used to make distinctions between potentially more than one intersections between A's input edge and B's input edge
  double tEdgeA = -1.;
  double tEdgeB = -1.;

  bool operator==(const OverlayVertex& other) const {
    return
      vA == other.vA &&
      vB == other.vB &&
      eA == other.eA &&
      eB == other.eB &&
      eIntrinsic == other.eIntrinsic &&
      tEdgeA == other.tEdgeA &&
      tEdgeB == other.tEdgeB;
  }
  bool operator!=(const OverlayVertex& other) const { return !(*this == other); }
};

// Notes on original intrinsic edges (input edges preserved in the intrinsic mesh):
// Such input edges are represented by their corresponding intrinsic edges. This reduces the type of edge
// intersection into the following three:
//    - intrinsic-inputA
//    - intrinsic-inputB
//    - inputA-inputB

struct OverlayPolygon;
struct OverlayWedge {
  OverlayVertex overlayVertex;
  Halfedge halfedgeIn;
  Halfedge halfedgeOut;
  double tHalfedgeIn = 0.;
  double tHalfedgeOut = 0.;
  SurfacePoint facePointA;
  SurfacePoint facePointB;

  OverlayPolygon* overlayPolygon_ptr = nullptr;     // Non-owning

  bool operator==(const OverlayWedge& other) const {
    return
      overlayVertex == other.overlayVertex &&
      halfedgeIn == other.halfedgeIn &&
      halfedgeOut == other.halfedgeOut &&
      tHalfedgeIn == other.tHalfedgeIn &&
      tHalfedgeOut == other.tHalfedgeOut &&
      facePointA == other.facePointA &&
      facePointB == other.facePointB;
  }
  bool operator!=(const OverlayWedge& other) const { return !(*this == other); }
};

// >>>>> To make OverlayVertex/OverlayWedge hashable >>>>>
} // namespace cit
namespace std {
  template <> struct hash<cit::OverlayVertex> {
    size_t operator()(const cit::OverlayVertex& ov) const {
      size_t seed = 0;
      hash_combine(seed, ov.vA);
      hash_combine(seed, ov.vB);
      hash_combine(seed, ov.eA);
      hash_combine(seed, ov.eB);
      hash_combine(seed, ov.eIntrinsic);
      hash_combine(seed, ov.tEdgeA);
      hash_combine(seed, ov.tEdgeB);
      return seed;
    }
  };
  template <> struct hash<cit::OverlayWedge> {
    size_t operator()(const cit::OverlayWedge& ow) const {
      size_t seed = 0;
      hash_combine(seed, ow.overlayVertex);
      hash_combine(seed, ow.halfedgeIn);
      hash_combine(seed, ow.halfedgeOut);
      hash_combine(seed, ow.tHalfedgeIn);
      hash_combine(seed, ow.tHalfedgeOut);
      hash_combine(seed, ow.facePointA);
      hash_combine(seed, ow.facePointB);
      return seed;
    }
  };
}
namespace cit{
// <<<<< To make OverlayVertex/OverlayWedge hashable <<<<<

struct OverlayPolygon {
  std::vector<OverlayWedge> wedges;
  Face fIntrinsic;
};

struct CoInTri {
  ModelData* mdata = nullptr;

  std::shared_ptr<SignpostIntrinsicTriangulation> signpostTri;

  VertexData<int> uniqueID_per_Vertex;

  // >>>>> Fields computed by updateCorrespondence() >>>>>
  std::map<int, Vertex> uniqueID_to_Vertex;

  // Stores null element if corresponding element doesn't exist
  EdgeData<Edge> correspondingEdge;
  VertexData<Vertex> correspondingVertex;
  HalfedgeData<Halfedge> correspondingHalfedge;
  FaceData<Face> correspondingFace;

  // Reversed faces are in correspondence in a weak sense that the two faces in A & B span the same set of vertices
  FaceData<Face> reversedFace;
  // <<<<< Fields computed by updateCorrespondence() <<<<<

  // Obtained by simply performing cointri.signpostTri->traceEdges()
  EdgeData<std::vector<SurfacePoint>> intrinsicEdgePath;

  // Obtained by rearranging intrinsicEdgePath, while making correspondence between edge points on input/intrinsic meshes
  EdgeData<std::vector<SurfacePoint>> inputEdgePathOnIntrinsic;
  std::unordered_map<SurfacePoint, SurfacePoint> edgePointMap_input2intrinsic;
  std::unordered_map<SurfacePoint, SurfacePoint> edgePointMap_intrinsic2input;

  // Obtained by tracing geodesics on the other input mesh, going through points corresponding to those stored in inputEdgePathOnIntrinsic
  EdgeData<std::vector<SurfacePoint>> inputEdgePathOnOtherInput;

  EdgeData<Matrix2d> laplacianTransportOperator;

  VertexData<std::vector<SurfacePoint>> vertexPath;
  VertexData<Matrix2d> temporalTransportOperator;
};

struct Configuration {
  CoInTri cointriA, cointriB;
  bool topologyValid = false;
  // bool angleValid = false;

  VectorXsp x;
  double energy = std::numeric_limits<double>::infinity();
  FaceData<double> energyDensity;
  VectorXd gradient;
  SparseMatrixd hessian;

  SparseMatrixd S, L, P;              // Laplacian preconditioner: P = (S*L)^T * M * S * L

  SparseMatrixd T, Tinv;              // Temporal transport matrix: g' = T * g,    H' = T^-T * H * T^-1
  VectorXd smoothedGradient;
  SparseMatrixd smoothedHessian;

  // Wedges for every overlay vertex corresponding to either:
  //    intrinsicMesh vertex, or
  //    edge intersections (intrinsicMesh-inputMesh{A,B} or inputMeshA-inputMeshB)
  std::unordered_map<OverlayVertex, std::vector<OverlayWedge>> overlayWedgesPerOverlayVertex;

  std::list<OverlayPolygon> overlayPolygons;

  bool edgePathComputed = false;
  bool derivativeComputed = false;
  bool laplacianTransportOperatorComputed = false;
  bool preconditionerComputed = false;

  size_t nCompatibleEdges = 0;
  size_t nCompatibleFaces = 0;

  uint64_t id = 0;
  uint64_t previous_id = 0;
};
extern Configuration currentConfig, initialConfig;

struct ResultInfo {
  size_t nIncompatibleEdges = 0;
  size_t nIncompatibleFaces = 0;

  double energy = 0;
  double energyDelta = 0;

  double A_minAngle = 0;
  double B_minAngle = 0;

  double A_maxSignpostError = 0;
  double B_maxSignpostError = 0;

  double previous_wL = 0;
  int previous_l = 0;
  double previous_gradientNorm = 0;
  double previous_descentDirectionNorm = 0;
  double previous_slope = 0;
  double stepSize = 0;
};
extern std::map<uint64_t, ResultInfo> resultInfoTable;


//-------------------+
// other global data |
//-------------------+

extern double current_wL;             // weight for the Laplacian preconditioner term
extern double current_smax;
extern std::list<std::pair<VectorXd, SparseMatrixd>> currentHistory;
extern size_t current_nSteps;

extern SparseMatrixd massMatrix;      // combined mass matrix (2*(mdataA.nV+mdataB.nV) rows/cols)

extern Camera* camera_active;
extern int window_width;
extern int window_height;
extern ImGuiIO* io;
extern uint64_t sessionID;
extern kt84::TextureObject colormapTexture;
extern std::string caseName;

extern int highlightedVertexID;
extern std::pair<int, int> highlightedEdgeID;
extern int highlightedFaceID;     // note: this ID is bare index of intrinsic mesh face

extern const char* const VERSIONTAG;


//-------------------+
// utility functions |
//-------------------+

void readInputMesh(ModelData& mdata, const std::string& inputFilename);

VectorXsp readBarycentricCoordinates(ModelData& mdataFrom, const ModelData& mdataTo, const std::string& filename);

std::pair<int, int> getUniqueEdgeID(const CoInTri& cointri, const std::array<Vertex, 2>& vp);

std::pair<int, int> getUniqueEdgeID(const CoInTri& cointri, Edge e);

Edge getEdgeByUniqueID(const CoInTri& cointri, const std::pair<int, int>& ue);

void setCameraMatrix(const ModelData& mdata);

// To be used when comparing elements belonging to different meshes (operator==() doesn't look at the referenced mesh)
template<typename T, typename M>
inline bool strictlyEqual(const Element<T, M>& e1, const Element<T, M>& e2);

inline bool isEdgeCompatible(const CoInTri& cointri, Edge e);

inline bool isFaceCompatible(const CoInTri& cointri, Face f);

inline bool isFaceReversed(const CoInTri& cointri, Face f);

template <typename T>
inline Vector3 getRandomColor(const T& obj);

std::array<Vector2, 4> getLayoutDiamond(const IntrinsicGeometryInterface& geometry, Halfedge iHe);

template <typename Vector2T>
typename Vector2T::Scalar computeTriangleEnergy(const std::array<Vector2T, 3>& A_triangle, const std::array<Vector2T, 3>& B_triangle, typename Vector2T::Scalar* A_area = nullptr, typename Vector2T::Scalar* B_area = nullptr);

template <typename Vector2T>
inline typename Vector2T::Scalar cross(const Vector2T& u, const Vector2T& v);

inline double angleAbsDifference(double angle1, double angle2);

inline double standardizeAngle(double angle, double modulo = 2. * M_PI);

// Collect all underlying input mesh faces that overlap with intrinsic_f, and lay them out consistently in 2D
std::unordered_map<Vertex, Vector2> layoutBaseFaces(const CoInTri& cointri, Face intrinsic_f);

ADScalar_6v_2nd computeTriangleEnergyAD(Configuration& config, Face A_f);

inline Matrix6d getHessian(const ADScalar_6v_2nd& y);

template <typename MatrixT>
inline MatrixT projectHessian(const MatrixT& H);

void setSolutionVector(Configuration& config);

bool isEdgeFlippable(const CoInTri& cointri, Edge e);

void checkVertexDegree(const CoInTri& cointri);

inline std::array<Vertex, 2> getOppositeVertices(Edge e);

inline SurfacePoint convertEdgePointToFacePoint(SurfacePoint edgePoint);

inline SurfacePoint getCorrespondingEdgePoint(const CoInTri& cointriP, SurfacePoint P_intrinsicSP);

template <typename T>
inline T* staticValuePtr(T value);  // for passing static values to ImGui using SliderScalar etc

void writeBlobToFile(const std::string& filename, const std::string& blob);

std::string readBlobFromFile(const std::string& filename);

std::vector<PatchPair> extractInconsistentPatches(const Configuration& config, bool useVertexAdjacency);

void analyzeInconsistentPatch(const CoInTri& cointri, Patch& patch);

inline bool isVertexMerged(int uniqueVertexID);

inline bool isVertexMerged(const CoInTri& cointri, Vertex v);

inline bool isVertexMovable(const CoInTri& cointri, Vertex v);

inline double getRadiusRatio(const SignpostIntrinsicTriangulation& geometry, Face f);

double getMinAngle(const CoInTri& cointri);

enum struct ResultContext {
  UpdateConfiguration = 0,
  FindMaxStep,
  FindOptimalStep,
};

inline std::string getResultFileNameBase();

inline std::string getResultFileName(ResultContext context, uint64_t id, uint64_t previous_id = 0, int l = 0);

void fillResultInfo(Configuration& newConfig, double currentConfigEnergy, const VectorXd& previousGradient, const VectorXd& previousDescentDirection, int l, double stepSize);

void writeResult(const Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, double smax, size_t nSteps, ResultContext context);

void writeOverlayMesh(const Configuration& config);

uint64_t makeConfigurationID();

spdlog::logger_ptr makeLogger(const std::string& filename, bool alongWithDefault = false, std::ostringstream* additional_oss = nullptr);

// Finds a strip of faces connecting fFrom and fTo; halfedges are oritented such that
// fFrom & fTo equal to res.front().face() & res.back().twin().face(), respectively
std::vector<Halfedge> findFaceStrip(Face fFrom, Face fTo);

bool isPatchSizeConsistent(const Patch& p1, const Patch& p2);

// Ensure that:
//    - If intrinsic edge is original on either P or Q, make it oriented the same as its corresponding input edge on P or Q
//    - For non-original intrinsic edge e, cointri.uniqueID_per_Vertex[e.halfedge().vertex()] < cointri.uniqueID_per_Vertex[e.halfedge().tipVertex()]
void makeIntrinsicEdgesConsistentlyOriented(CoInTri& cointriP, const CoInTri& cointriQ);

void copyConfiguration(const Configuration& origConfig, Configuration& copiedConfig, size_t logDepth = 0);

// If given face point is very close to an edge (with the corresponding barycentric coordinate component below
// sysParam.snapThreshold), snap it to the edge and return true
bool snapFacePointToEdge(SurfacePoint& sp);

// Estimate reasonable maximum step size based on the magnitude of descent direction (to avoid useless long geodesic tracing)
double estimateMaxStep(Configuration& config, double wL, size_t logDepth = 0);

// Useful for debugging
std::set<Edge> getCompatibleEdges(const CoInTri& cointri);
std::set<Edge> getIncompatibleEdges(const CoInTri& cointri);
std::set<Face> getCompatibleFaces(const CoInTri& cointri);
std::set<Face> getIncompatibleFaces(const CoInTri& cointri);
std::set<Face> getReversedFaces(const CoInTri& cointri);
std::vector<Vertex> adjacentVertices(Vertex v);
std::vector<Halfedge> outgoingHalfedges(Vertex v);
std::vector<Face> adjacentFaces(Vertex v);
std::array<Face, 2> adjacentFaces(Edge e);
std::array<Vertex, 3> adjacentVertices(Face f);
std::array<Halfedge, 3> adjacentHalfedges(Face f);
std::array<Face, 3> adjacentFaces(Face f);

std::string guessNiceNameFromPath(std::string fullname);  // copied from https://github.com/nmwsharp/polyscope/blob/master/src/utilities.cpp


//---------------------+
// core mesh functions |
//---------------------+

// P_x is of size cointriP.mdata->nV, each representing a point on cointriQ
void insertVertices(const VectorXsp& P_x, const CoInTri& cointriP, CoInTri& cointriQ, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

void updateCorrespondence(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

// For every compatible edge between original & inserted vertices, if they're very close (based on snapThreshold)
// both on A & B, merge them
size_t mergeNearbyVertexPairs(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t simpleFlip(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t simpleCoFlip(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t flipCompatible(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t doubleFlipCompatible(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t flipFlatPolygon(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

// flip edges to make A and B as compatible as possible. return the number of flipped edges
size_t flipToCompatible(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

// Flip compatible edges within incompatible patches as much as possible to improve compatibility
size_t flipPatchInteriorCompatibleEdgesToExploreVariations(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t convexifyConcaveIncompatiblePatch(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t mergeInconsistentVertexPairs(Configuration& config, bool useVertexAdjacency, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

bool splitMergedVertices(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

OverlayVertex getOverlayVertexAtEdgeIntersection(const Configuration& config, std::array<Edge, 2> e);

void storeOverlayWedgesAtEdgeIntersection(Configuration& config, std::array<SurfacePoint, 2> edgePoint, size_t logDepth = 0);

void computeEdgePath(Configuration& config, size_t logDepth = 0);

void computeOverlayPolygons(Configuration& config, size_t logDepth = 0);

double computeEnergy(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

void computeDerivative(Configuration& config, size_t logDepth = 0);

void computeLaplacianTransportOperator(Configuration& config, size_t logDepth = 0);

void computeTemporalTransportOperator(Configuration& config, size_t logDepth = 0);

void computePreconditioner(Configuration& config, size_t logDepth = 0);

VectorXd computeDescentDirection(Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, int l, size_t logDepth = 0);

void flipToMinimizeEnergy(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

void flipToMaximizeMinAngle(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

size_t relocateVerticesToRemoveDegenerateFaces(Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

// moves solution x by delta, which may involve crossing different triangles,
// thus the barycentric coordinates are transformed appropriately using parallel transport (see [Schmidt,SIG20] supp mat)
VectorXsp displaceSolutionVector(const VectorXsp& x, const VectorXd& delta, VertexData<std::vector<SurfacePoint>>& A_vertexPath, VertexData<std::vector<SurfacePoint>>& B_vertexPath, size_t logDepth = 0);

// Creates a configuration from given solution vector by inserting vertices.
// Returns false if x is infeasible (i.e. failed to arrive at cointri).
bool createConfiguration(uint64_t id, VectorXsp x, Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);

// Creates a configuration from given previous configuration and delta vector by relocating vertices.
// If it fails, falls back to the vertex insertion mode.
bool createConfiguration(uint64_t id, const Configuration& prevConfig, const VectorXd& delta, Configuration& config, spdlog::logger_ptr logger = spdlog::default_logger(), size_t logDepth = 0);


//-----------+
// optimizer |
//-----------+

// Find maximum step size by starting from current_smax, store result into current_smax
bool findMaxStep(const VectorXd& descentDirection, int l);

// Given current_smax, perform backtrack line search along descentDirection
bool findOptimalStep(const VectorXd& descentDirection, int l, double& s, Configuration& newConfig);

void updateConfiguration(const Configuration& newConfig, int l_chosen);

bool optimizeOneStep(bool doTemporalSmoothing);

// Optimize N steps (when -1, indefinitely until convergence)
void optimize(int N);

void optimizerMenu();


//----------------+
// event handling |
//----------------+

void callback_framebuffersize(GLFWwindow* window, int width, int height);

void callback_mousebutton(GLFWwindow* window, int button, int action, int mods);

void callback_cursorpos(GLFWwindow* window, double xpos, double ypos);

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods);


//-----------+
// rendering |
//-----------+

void draw();


//---------------+
// serialization |
//---------------+

std::string toSerializedBlob(const Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, double smax, bool firstOrderMode, size_t nSteps);
void fromSerializedBlob(const std::string& serializedBlob, Configuration& config, std::list<std::pair<VectorXd, SparseMatrixd>>& history, double& wL, double& smax, bool& firstOrderMode, size_t& nSteps);


}

#include "utility.ipp"
