#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_point.h"
#include "args/args.hxx"

#include "kt84/MaxMinSelector.hh"
#include "kt84/eigen_util.hh"

#include <spdlog/spdlog.h>

#define CIT_ASSERT(CONDITION) GC_SAFETY_ASSERT((CONDITION), "");

#define LOG_TRACE(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::trace   , __VA_ARGS__)
#define LOG_DEBUG(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::debug   , __VA_ARGS__)
#define LOG_INFO(...)     spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::info    , __VA_ARGS__)
#define LOG_WARN(...)     spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::warn    , __VA_ARGS__)
#define LOG_ERROR(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::err     , __VA_ARGS__)
#define LOG_CRITICAL(...) spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::critical, __VA_ARGS__)

#include <fstream>
#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::AlignedBox3d;
using kt84::MinSelector;
using kt84::MaxMinAverage;

using VectorXsp = Eigen::Matrix<SurfacePoint, -1, 1>;

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

MAKE_FORMATTABLE(Vertex);
MAKE_FORMATTABLE(SurfacePoint);
MAKE_FORMATTABLE(Vector3);
MAKE_FORMATTABLE(MaxMinAverage);




struct ModelData {
  std::string shortName;

  std::unique_ptr<ManifoldSurfaceMesh> inputMesh;
  std::unique_ptr<VertexPositionGeometry> imGeometry;

  std::unique_ptr<ManifoldSurfaceMesh> overlayMesh;
  std::unique_ptr<VertexPositionGeometry> omGeometry;

  VertexData<Vertex> omVertex_per_imVertex;
  VertexData<Vertex> omVertex_to_imVertex;        // Vertex() if not exist

  AlignedBox3d bbox;
} mdataA, mdataB;




void initializeModelData(ModelData& mdata, const std::string& filename_inputMesh, const std::string& filename_overlayMesh) {
  LOG_INFO("Reading meshes and initializing correspondence data for {}", mdata.shortName);

  LOG_INFO("  Reading input mesh ...");
  std::tie(mdata.inputMesh, mdata.imGeometry, std::ignore) = readManifoldSurfaceMesh(filename_inputMesh);
  LOG_INFO("  Mesh stats:");
  LOG_INFO("    nV = {}", mdata.inputMesh->nVertices());
  LOG_INFO("    nF = {}", mdata.inputMesh->nFaces());
  LOG_INFO("  Reading overlay mesh ...");
  std::tie(mdata.overlayMesh, mdata.omGeometry, std::ignore) = readManifoldSurfaceMesh(filename_overlayMesh);
  LOG_INFO("  Mesh stats:");
  LOG_INFO("    nV = {}", mdata.overlayMesh->nVertices());
  LOG_INFO("    nF = {}", mdata.overlayMesh->nFaces());

  mdata.omVertex_per_imVertex = VertexData<Vertex>(*mdata.inputMesh);
  mdata.omVertex_to_imVertex = VertexData<Vertex>(*mdata.overlayMesh);

  for (Vertex im_v : mdata.inputMesh->vertices())
    mdata.bbox.extend(mdata.imGeometry->inputVertexPositions[im_v].castTo<Vector3d>());
  LOG_INFO("  Bounding box diagonal (input): {}", Vector3::castFrom(mdata.bbox.diagonal()));

  mdata.bbox.setEmpty();
  for (Vertex om_v : mdata.overlayMesh->vertices())
    mdata.bbox.extend(mdata.omGeometry->inputVertexPositions[om_v].castTo<Vector3d>());
  LOG_INFO("  Bounding box diagonal (overlay): {}", Vector3::castFrom(mdata.bbox.diagonal()));

  mdata.imGeometry->requireEdgeLengths(); // For debugging

  mdata.imGeometry->requireFaceAreas();
  double surfaceArea = 0.;
  for (Face im_f : mdata.inputMesh->faces())
    surfaceArea += mdata.imGeometry->faceAreas[im_f];
  LOG_INFO("  Surface area : {}", surfaceArea);

  mdata.omGeometry->requireFaceAreas();
  surfaceArea = 0.;
  for (Face om_f : mdata.overlayMesh->faces())
    surfaceArea += mdata.omGeometry->faceAreas[om_f];
  LOG_INFO("  Surface area : {}", surfaceArea);

  MaxMinAverage mma;
  LOG_INFO("Finding correspondence among input & overlay mesh vertices");
  for (Vertex im_v : mdata.inputMesh->vertices()) {
    if ((im_v.getIndex() + 1) % (mdata.inputMesh->nVertices() / 10) == 0)
      LOG_DEBUG("  {}", im_v);

    MinSelector<Vertex> om_vClosest;
    for (Vertex om_v : mdata.overlayMesh->vertices())
      om_vClosest.update(norm(mdata.imGeometry->inputVertexPositions[im_v] - mdata.omGeometry->inputVertexPositions[om_v]), om_v);
    CIT_ASSERT(om_vClosest.score < 1.e-16);
    mma.update(om_vClosest.score);

    mdata.omVertex_per_imVertex[im_v] = om_vClosest.value;
    mdata.omVertex_to_imVertex[om_vClosest.value] = im_v;
  }
  LOG_INFO("  Error stats: {}", mma);
}

VectorXsp extractBarycentricCoordinates(const ModelData& mdataP, const ModelData& mdataQ) {
  LOG_INFO("Extracting barycentric coordinates for vertices in {} onto faces in {}", mdataP.shortName, mdataQ.shortName);
  VectorXsp res;
  res.resize(mdataP.inputMesh->nVertices());

  MaxMinAverage mma;
  for (Vertex P_im_v : mdataP.inputMesh->vertices()) {
    if ((P_im_v.getIndex() + 1) % (mdataP.inputMesh->nVertices() / 10) == 0)
      LOG_DEBUG("  {}", P_im_v);

    Vertex P_om_v = mdataP.omVertex_per_imVertex[P_im_v];
    Vertex Q_om_v = mdataQ.overlayMesh->vertex(P_om_v.getIndex());
    Vertex Q_im_v = mdataQ.omVertex_to_imVertex[Q_om_v];

    if (Q_im_v != Vertex()) {
      // P_im_v corresponds to Q_im_v
      res[P_im_v.getIndex()] = SurfacePoint(Q_im_v);
      LOG_DEBUG("  Vertex {} in {} corresponds to vertex {} in {}", P_im_v, mdataP.shortName, Q_im_v, mdataQ.shortName);

    } else {
      // P_im_v corresponds to a face point on Q
      Vector3d p = mdataQ.omGeometry->inputVertexPositions[Q_om_v].castTo<Vector3d>();

      MinSelector<SurfacePoint> Q_spClosest;
      for (Face Q_im_f : mdataQ.inputMesh->faces()) {
        // Get 3D points of Q's inputMesh face triangle
        std::array<Vector3d, 3> pTri;
        int i = 0;
        for (Vertex Q_im_v2 : Q_im_f.adjacentVertices())
          pTri[i++] = mdataQ.imGeometry->inputVertexPositions[Q_im_v2].castTo<Vector3d>();

        Vector3d t;
        const double safeMargin = 1.e-4;
        if (kt84::eigen_util::project_to_triangle(pTri[0], pTri[1], pTri[2], p, t) && (t.array() > -safeMargin).all()) {
          // >>>>> Project t to valid region >>>>>
          // First, project to a corner if greater than 1
          bool found = false;
          for (int j = 0; j < 3; ++j) {
            if (t[j] >= 1.) {
              CIT_ASSERT(!found);     // There can be just one such component
              found = true;
              t[j] = 1.;
              t[(j + 1) % 3] = 0.;
              t[(j + 2) % 3] = 0.;
            }
          }
          // If above isn't the case, project to an edge if negative
          if (!found) {
            for (int j = 0; j < 3; ++j) {
              if (t[j] <= 0.) {
                CIT_ASSERT(!found);   // There can be just one such component
                found = true;
                double delta = -t[j];
                t[j] = 0.;
                t[(j + 1) % 3] -= 0.5 * delta;
                t[(j + 2) % 3] -= 0.5 * delta;
              }
            }
          }
          // <<<<< Project t to valid region <<<<<

          Vector3d pOnTri = Vector3d::Zero();
          for (int j = 0; j < 3; ++j)
            pOnTri += t[j] * pTri[j];

          Q_spClosest.update((pOnTri - p).norm(), SurfacePoint(Q_im_f, Vector3::castFrom(t)).reduced());
        }
      }

      CIT_ASSERT(Q_spClosest.count);
      CIT_ASSERT(Q_spClosest.score < 1.e-5);
      mma.update(Q_spClosest.score);
      res[P_im_v.getIndex()] = Q_spClosest.value;
      LOG_TRACE("  Vertex {} in {} corresponds to a surface point {} in {}", P_im_v, mdataP.shortName, Q_spClosest.value, mdataQ.shortName);
    }
    CIT_ASSERT(res[P_im_v.getIndex()] != SurfacePoint());
  }
  LOG_INFO("  Error stats: {}", mma);

  return res;
}

void checkVertexCorrespondenceConsistency(const ModelData& mdataP, const ModelData& mdataQ, const VectorXsp& baryCoords_PtoQ, VectorXsp& baryCoords_QtoP) {
  LOG_INFO("Checking vertex-vertex map consistency: from {} to {}", mdataP.shortName, mdataQ.shortName);
  for (Vertex P_im_v : mdataP.inputMesh->vertices()) {
    const SurfacePoint& Q_sp = baryCoords_PtoQ[P_im_v.getIndex()];
    if (Q_sp.type == SurfacePointType::Vertex) {
      // Get the corresponding vertex on Q, then its mapped point on P
      Vertex Q_im_v = Q_sp.vertex;
      SurfacePoint& P_sp = baryCoords_QtoP[Q_im_v.getIndex()];
      LOG_DEBUG("{}'s vertex {} is mapped to {}'s vertex {}", mdataP.shortName, P_im_v, mdataQ.shortName, Q_im_v);

      // If the other map is also onto the same vertex, then good
      if (P_sp.type == SurfacePointType::Vertex) {
        CIT_ASSERT(P_sp.vertex == P_im_v);
        LOG_DEBUG("The other mapping is also onto the same vertex");
        continue;
      }

      // Otherwise, check adjacency & proxmity, then snap
      CIT_ASSERT(checkAdjacent(P_sp, SurfacePoint(P_im_v)));
      double d = norm(mdataP.imGeometry->inputVertexPositions[P_im_v] - P_sp.interpolate(mdataP.imGeometry->inputVertexPositions));
      CIT_ASSERT(d < 1.e-5);
      LOG_DEBUG("The other mapping is onto a nearby surface point {} (distance = {}), snapping", P_sp, d);
      P_sp = SurfacePoint(P_sp.nearestVertex());
    }
  }
}

void exportResult(const VectorXsp& result, const std::string& filename) {
  std::ofstream fout(filename.c_str());
  if (!fout.is_open())
    throw std::runtime_error("Couldn't open file " + filename);

  for (int i = 0; i < result.size(); ++i) {
    const SurfacePoint sp = result[i];

    if (sp.type == SurfacePointType::Vertex) {
      fout << (i + 1) << " " << (sp.vertex.getIndex() + 1) << " 1" << std::endl;

    } else if (sp.type == SurfacePointType::Edge) {
      fout << (i + 1) << " " << (sp.edge.halfedge().tailVertex().getIndex() + 1) << " " << (1. - sp.tEdge) << std::endl;
      fout << (i + 1) << " " << (sp.edge.halfedge().tipVertex ().getIndex() + 1) << " " << sp.tEdge        << std::endl;

    } else {
      int j = 0;
      for (Vertex v : sp.face.adjacentVertices())
        fout << (i + 1) << " " << (v.getIndex() + 1) << " " << sp.faceCoords[j++] << std::endl;
    }
  }
  LOG_INFO("Result written to {}", filename);
}

int main(int argc, char** argv) {
  // Configure the argument parser
  // clang-format off
  args::ArgumentParser parser("Barycentric Coordinate Extractor");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  std::string logLevelHelp;
  {
    std::ostringstream oss;
    oss << "Log level {";
    bool first = true;
    for (auto s : SPDLOG_LEVEL_NAMES) {
      if (!first)
        oss << ", ";
      first = false;
      oss << s.data();
    }
    oss << "}";
    logLevelHelp = oss.str();
  }
  args::ValueFlag<std::string> logLevel(parser, "level", logLevelHelp, {"log-level"});

  args::Positional<std::string> argMeshA(parser, "meshA", "The mesh of model A");
  args::Positional<std::string> argMeshB(parser, "meshB", "The mesh of model B");
  args::Positional<std::string> argMeshMonA(parser, "meshMonA", "The overlay mesh on model A");
  args::Positional<std::string> argMeshMonB(parser, "meshMonB", "The overlay mesh on model B");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!argMeshA || !argMeshB || !argMeshMonA || !argMeshMonB) {
    std::cout << parser;
    return EXIT_FAILURE;
  }

  // Configure logging
  spdlog::set_pattern("%D %T %L %s:%# %v");
  spdlog::flush_on(spdlog::level::debug);
  if (logLevel) {
    int i = 0;
    for (auto s : SPDLOG_LEVEL_NAMES) {
      if (args::get(logLevel) == s) {
        spdlog::set_level((spdlog::level::level_enum)i);
        break;
      }
      ++i;
    }
    if (i == spdlog::level::n_levels) {
      LOG_ERROR("Wrong log level was specified");
      std::cout << parser;
      return EXIT_FAILURE;
    }
  }

  mdataA.shortName = "A";
  mdataB.shortName = "B";

  initializeModelData(mdataA, args::get(argMeshA), args::get(argMeshMonA));
  initializeModelData(mdataB, args::get(argMeshB), args::get(argMeshMonB));

  VectorXsp baryCoords_AtoB = extractBarycentricCoordinates(mdataA, mdataB);
  VectorXsp baryCoords_BtoA = extractBarycentricCoordinates(mdataB, mdataA);

  checkVertexCorrespondenceConsistency(mdataA, mdataB, baryCoords_AtoB, baryCoords_BtoA);
  checkVertexCorrespondenceConsistency(mdataB, mdataA, baryCoords_BtoA, baryCoords_AtoB);

  exportResult(baryCoords_AtoB, "AtoB.txt");
  exportResult(baryCoords_BtoA, "BtoA.txt");

  return EXIT_SUCCESS;
}
