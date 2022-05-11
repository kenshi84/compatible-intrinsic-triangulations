#include <iostream>

#include "kt84/MaxMinSelector.hh"
#include "kt84/graphics/DisplayList.hh"
#include "kt84/graphics/TextureObjectT.hh"
#include "kt84/graphics/graphics_util.hh"
#include "kt84/geometry/CameraFree.hh"
#include "kt84/glfw_util.hh"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/utilities/elementary_geometry.h"
#include "args/args.hxx"

#include "examples/imgui_impl_glfw.h"
#include "examples/imgui_impl_opengl2.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreorder"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#include "imGuIZMOquat.h"
#pragma clang diagnostic pop

#define STB_IMAGE_IMPLEMENTATION
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"
#include "stb/stb_image.h"
#pragma clang diagnostic pop

#include "cit_resources.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#define CIT_ASSERT(CONDITION) GC_SAFETY_ASSERT((CONDITION), "");

#define LOG_TRACE(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::trace   , __VA_ARGS__)
#define LOG_DEBUG(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::debug   , __VA_ARGS__)
#define LOG_INFO(...)     spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::info    , __VA_ARGS__)
#define LOG_WARN(...)     spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::warn    , __VA_ARGS__)
#define LOG_ERROR(...)    spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::err     , __VA_ARGS__)
#define LOG_CRITICAL(...) spdlog::default_logger_raw()->log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION}, spdlog::level::critical, __VA_ARGS__)

#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace kt84::graphics_util;

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector4f;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::AlignedBox3d;
using kt84::MaxMinAverage;
using kt84::Camera;

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
MAKE_FORMATTABLE(Edge);
MAKE_FORMATTABLE(Halfedge);
MAKE_FORMATTABLE(Face);
MAKE_FORMATTABLE(Corner);
MAKE_FORMATTABLE(Vector2);
MAKE_FORMATTABLE(Vector3);
MAKE_FORMATTABLE(MaxMinAverage);



//-------------+
// Global data |
//-------------+

struct ModelData {
  std::string shortName;

  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::unique_ptr<CornerData<Vector2>> texCoords;

  CornerData<Vector2> cornerPositionsInFace;
  FaceData<double> faceAreas;

  AlignedBox3d bbox;
  kt84::TextureObject texture;
  kt84::CameraFree camera;
  kt84::DisplayList dispList;
} mdataA, mdataB;

size_t nV, nE, nF;

FaceData<Matrix2d> jacobianAtoB;
FaceData<Matrix2d> jacobianBtoA;
FaceData<double> energyPerFace;

Camera* camera_active = nullptr;
int window_width;
int window_height;
kt84::TextureObject colormapTexture;

const std::array<std::string, 3> fillModes = {"Distortion", "Texture", "Constant"};

struct Tweaks {
  struct {} i;
  struct {
    bool showFaces = true;
    bool showEdges = true;
    bool showBoundingBox = true;
  } b;
  struct {
    float edgeWidth = 1.f;
    float textureScale = 1.f;
    float logMaxEnergyDensity = 3.f;
  } f;
  struct {} d;
  struct {
    Vector3f faceColor = Vector3f::Constant(0.5f);
  } vec3;
  struct {
    std::string fillModeA = fillModes[0];
    std::string fillModeB = fillModes[0];
  } s;
} tweaks;

struct LightParam {
  bool enabled = false;
  float ambient = 0.f;
  float diffuse = 0.f;
  float specular = 0.f;
  Vector4f position = {0.f, 0.f, 1.f, 0.f};
};
std::array<LightParam, 8> lightParam;

struct MaterialParam {
  float ambient = 0.2f;
  float specular = 0.5f;
  float shininess = 20.f;
} materialParam;

//-----------+
// Utilities |
//-----------+

inline std::vector<Corner> getAdjacentCorners(Face f) {
  size_t d = f.degree();
  std::vector<Corner> res(d);
  size_t i = 0;
  for (Corner c : f.adjacentCorners())
    res[i++] = c;
  return res;
}

void initializeCamera(ModelData& mdata) {
  Vector3d center = 0.5 * (mdata.bbox.min() + mdata.bbox.max());
  mdata.camera.init(center + Vector3d(0,0,mdata.bbox.diagonal().norm() * 1.2), center, Vector3d::UnitY());
}

void setCameraMatrix(const ModelData& mdata) {
  if (&mdata == &mdataA)
    glViewport(0, 0, window_width / 2, window_height);
  else
    glViewport(window_width / 2, 0, window_width / 2, window_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double zNear = mdata.camera.center_to_eye().norm() * 0.1;
  double zFar  = zNear * 10 + mdata.bbox.diagonal().norm() * 10;
  double aspect_ratio = window_width / 2 / (double)window_height;
  gluPerspective(40, aspect_ratio, zNear, zFar);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Set light direction in the camera coordinate
  for (size_t i = 0; i < 8; ++i)
    glLightPosition4f(lightParam[i].position, i);

  gluLookAt(mdata.camera.get_eye(), mdata.camera.center, mdata.camera.get_up());
}

void invalidateDispList() {
  mdataA.dispList.invalidate();
  mdataB.dispList.invalidate();
}

//----------------+
// Core functions |
//----------------+

void readMesh(ModelData& mdata, const std::string& filename) {
  LOG_INFO("Reading mesh ... ({})", mdata.shortName);

  std::tie(mdata.mesh, mdata.geometry, mdata.texCoords) = readManifoldSurfaceMesh(filename);

  LOG_INFO("  Computing bounding box ...");
  for (Vertex v : mdata.mesh->vertices())
    mdata.bbox.extend(mdata.geometry->inputVertexPositions[v].castTo<Vector3d>());

  LOG_INFO("  Mesh stats:");
  LOG_INFO("    #V: {}", mdata.mesh->nVertices());
  LOG_INFO("    #E: {}", mdata.mesh->nEdges());
  LOG_INFO("    #F: {}", mdata.mesh->nFaces());
  LOG_INFO("    #B: {}", mdata.mesh->nBoundaryLoops());
  LOG_INFO("    Bounding box diagonal: {} (norm: {})", Vector3::castFrom(mdata.bbox.diagonal()), mdata.bbox.diagonal().norm());
  LOG_INFO("    Textured: {}", mdata.texCoords ? "yes" : "no");

  // Check degenerate edges
  mdata.geometry->requireEdgeLengths();
  mdata.geometry->requireFaceNormals();

  LOG_INFO("  Checking edge degenecary ...");
  for (Edge e : mdata.mesh->edges()) {
    if (mdata.geometry->edgeLengths[e] == 0.)
      LOG_ERROR("    Found degenerate edge: {} ({}, {})", e, e.halfedge().vertex(), e.halfedge().tipVertex());
  }
  LOG_INFO("  Done");
}

void checkTopologicalConsistency() {
  LOG_INFO("Checking topologycal consistency ...");
  bool good = true;

  if (mdataA.mesh->nVertices     () != mdataB.mesh->nVertices     ()) { LOG_ERROR("Inconsistent #V"); good = false; }
  if (mdataA.mesh->nEdges        () != mdataB.mesh->nEdges        ()) { LOG_ERROR("Inconsistent #E"); good = false; }
  if (mdataA.mesh->nFaces        () != mdataB.mesh->nFaces        ()) { LOG_ERROR("Inconsistent #F"); good = false; }
  if (mdataA.mesh->nBoundaryLoops() != mdataB.mesh->nBoundaryLoops()) { LOG_ERROR("Inconsistent #B"); good = false; }

  nV = mdataA.mesh->nVertices();
  nE = mdataA.mesh->nEdges();
  nF = mdataA.mesh->nFaces();

  for (size_t i = 0; i < nF; ++i) {
    Face fA = mdataA.mesh->face(i);
    Face fB = mdataB.mesh->face(i);

    if (fA.degree() != fB.degree()) {
      LOG_ERROR("Inconsistent face degree {} vs {} found in face {}", fA.degree(), fB.degree(), i);
      good = false;
      continue;
    }

    auto iA = fA.adjacentVertices().begin();
    auto iB = fB.adjacentVertices().begin();
    for (size_t j = 0; j < fA.degree(); ++j, ++iA, ++iB) {
      if ((*iA).getIndex() != (*iB).getIndex()) {
        LOG_ERROR("Inconsistent vertex id {} vs {} found in face {}", (*iA).getIndex(), (*iB).getIndex(), i);
        good = false;
      }
    }
  }
  LOG_INFO("  Done");
  CIT_ASSERT(good);
}

void computeCornerPositionsInFace(ModelData& mdata) {
  LOG_INFO("Computing 2D corner positions ({}) ...", mdata.shortName);

  mdata.cornerPositionsInFace = CornerData<Vector2>(*mdata.mesh, Vector2::undefined());
  MaxMinAverage mma;
  for (Face f : mdata.mesh->faces()) {
    size_t d = f.degree();
    std::vector<Corner> adjacentCorners = getAdjacentCorners(f);

    // Get 3D positions of face corners
    std::vector<Vector3> p(d);
    for (size_t i = 0; i < d; ++i)
      p[i] = mdata.geometry->inputVertexPositions[adjacentCorners[i].vertex()];

    // Get rotation matrix which brings p[0] to origin, p[1] to X axis, and p[2] on XY plane
    Vector3 dX = normalize(p[1] - p[0]);
    Vector3 dZ = normalize(cross(dX, p[2] - p[0]));
    Vector3 dY = cross(dZ, dX);

    if (!dX.isDefined() || !dZ.isDefined())
      continue;

    Matrix3d Minv;
    Minv << dX.castTo<Vector3d>(), dY.castTo<Vector3d>(), dZ.castTo<Vector3d>();
    Matrix3d M = Minv.inverse();

    mdata.cornerPositionsInFace[adjacentCorners[0]] = {0., 0.};
    for (size_t i = 1; i < d; ++i) {
      Vector3d q = M * (p[i] - p[0]).castTo<Vector3d>();
      mdata.cornerPositionsInFace[adjacentCorners[i]] = Vector2::castFrom(q);
      mma.update(std::fabs(q.z()));
    }
  }
  LOG_INFO("  Done");
  LOG_INFO("    Planarity error stats: {}", mma);
  // CIT_ASSERT(mma.max() < mdata.bbox.diagonal().norm() * 1.e-3);
}

void computeFaceAreas(ModelData& mdata) {
  LOG_INFO("Computing face areas ({}) ...", mdata.shortName);

  mdata.faceAreas = FaceData<double>(*mdata.mesh);
  for (Face f : mdata.mesh->faces()) {
    size_t d = f.degree();

    // If triangular, use elementary geometry
    if (d == 3) {
      double lA = mdata.geometry->edgeLengths[f.halfedge().edge()];
      double lB = mdata.geometry->edgeLengths[f.halfedge().next().edge()];
      double lC = mdata.geometry->edgeLengths[f.halfedge().next().next().edge()];
      mdata.faceAreas[f] = triangleArea(lA, lB, lC);

    // Otherwise, use 2D corner positions
    } else {
      std::vector<Corner> adjacentCorners = getAdjacentCorners(f);

      mdata.faceAreas[f] = 0.;
      for (size_t i = 2; i < d; ++i) {
        Corner c0 = adjacentCorners[i - 1];
        Corner c1 = adjacentCorners[i];

        Vector2 p0 = mdata.cornerPositionsInFace[c0];
        Vector2 p1 = mdata.cornerPositionsInFace[c1];
        CIT_ASSERT(p0.isDefined() && p1.isDefined());

        mdata.faceAreas[f] += cross(p0, p1) / 2.;
      }
    }

    if (mdata.faceAreas[f] == 0.) {
      std::string vid_str;
      for (Vertex v : f.adjacentVertices())
        vid_str += (vid_str.empty() ? "" : ", ") + std::to_string(v);
      LOG_WARN("  Found degenerate face: {} ({})", f, vid_str);
    }
  }
  LOG_INFO("  Done");
}

void normalizeSurfaceArea(ModelData& mdata) {
  LOG_INFO("Normalizing surface area ({}) ...", mdata.shortName);

  double surfaceArea = mdata.faceAreas.raw().sum();
  LOG_INFO("  Original surface area: {}", surfaceArea);

  for (Vertex v : mdata.mesh->vertices())
    mdata.geometry->inputVertexPositions[v] /= std::sqrt(surfaceArea);

  mdata.geometry->refreshQuantities();

  // Recompute face areas (with logging suppressed)
  LOG_INFO("  Recomputing face areas ...");
  spdlog::level::level_enum l = spdlog::get_level();
  spdlog::set_level(spdlog::level::level_enum::off);
  computeFaceAreas(mdata);
  spdlog::set_level(l);

  // Recompute bounding box
  LOG_INFO("  Recomputing bounding box ...");
  mdata.bbox.setEmpty();
  for (Vertex v : mdata.mesh->vertices())
    mdata.bbox.extend(mdata.geometry->inputVertexPositions[v].castTo<Vector3d>());

  LOG_INFO("  Done");
  LOG_INFO("    Normalized surface area: {}", mdata.faceAreas.raw().sum());
  LOG_INFO("    Bounding box diagonal after normalization: {} (norm: {})", Vector3::castFrom(mdata.bbox.diagonal()), mdata.bbox.diagonal().norm());
}

void computeJacobian() {
  LOG_INFO("Computing Jacobian ...");
  jacobianAtoB = FaceData<Matrix2d>(*mdataA.mesh);
  jacobianBtoA = FaceData<Matrix2d>(*mdataA.mesh);
  jacobianAtoB.setDefault(Matrix2d::Zero());
  jacobianBtoA.setDefault(Matrix2d::Zero());
  MaxMinAverage mma;
  for (size_t fid = 0; fid < nF; ++fid) {
    Face fA = mdataA.mesh->face(fid);
    Face fB = mdataB.mesh->face(fid);

    // Ignore degenerate faces
    if (mdataA.faceAreas[fA] == 0.) continue;
    if (mdataB.faceAreas[fB] == 0.) continue;

    std::vector<Corner> A_adjacentCorners = getAdjacentCorners(fA);
    std::vector<Corner> B_adjacentCorners = getAdjacentCorners(fB);

    size_t d = fA.degree();
    std::vector<Vector2> pA(d);
    std::vector<Vector2> pB(d);
    for (size_t i = 0; i < d; ++i) {
      pA[i] = mdataA.cornerPositionsInFace[A_adjacentCorners[i]];
      pB[i] = mdataB.cornerPositionsInFace[B_adjacentCorners[i]];
    }

    Vector2d A_d0 = (pA[1] - pA[0]).castTo<Vector2d>();
    Vector2d A_d1 = (pA[2] - pA[0]).castTo<Vector2d>();
    Vector2d B_d0 = (pB[1] - pB[0]).castTo<Vector2d>();
    Vector2d B_d1 = (pB[2] - pB[0]).castTo<Vector2d>();

    Matrix2d A_M;
    Matrix2d B_M;
    A_M << A_d0, A_d1;
    B_M << B_d0, B_d1;
    Matrix2d A_Minv = A_M.inverse();
    Matrix2d B_Minv = B_M.inverse();

    jacobianAtoB[fA] = B_M * A_Minv;
    jacobianBtoA[fA] = A_M * B_Minv;

    // For polygonal faces, verify that transforming fA's corner positions by this Jacobian agrees with fB's corner positions
    for (size_t i = 3; i < d; ++i) {
      Vector2 pBi_check = Vector2::castFrom<Vector2d>(jacobianAtoB[fA] * (pA[i] - pA[0]).castTo<Vector2d>()) + pB[0];
      mma.update(norm(pB[i] - pBi_check));
    }
  }
  LOG_INFO("  Done");
  if (mma.count())
    LOG_INFO("    Linear mapping error stats: {}", mma);
  // CIT_ASSERT(mma.max() < mdataB.bbox.diagonal().norm() * 0.005);
}

void computeEnergy() {
  LOG_INFO("Computing energy ...");
  energyPerFace = FaceData<double>(*mdataA.mesh, 0.);
  MaxMinAverage mma;
  for (size_t fid = 0; fid < nF; ++fid) {
    Face fA = mdataA.mesh->face(fid);
    Face fB = mdataB.mesh->face(fid);

    if (mdataA.faceAreas[fA] == 0.) continue;
    if (mdataB.faceAreas[fB] == 0.) continue;

    energyPerFace[fA] += mdataB.faceAreas[fB] * jacobianAtoB[fA].squaredNorm();
    energyPerFace[fA] += mdataA.faceAreas[fA] * jacobianBtoA[fA].squaredNorm();

    if (std::isnan(energyPerFace[fA])) {
      energyPerFace[fA] = 0.;
      continue;
    }

    // Energy density is obatined by dividing per-face energy by the average face area
    double energyDensity = energyPerFace[fA] / ((mdataA.faceAreas[fA] + mdataB.faceAreas[fB]) * 0.5);
    mma.update(energyDensity);
  }
  LOG_INFO("  Done");
  LOG_INFO("    Energy: {}", energyPerFace.raw().sum());
  LOG_INFO("    Energy density stats: {}, log(max)={}", mma, tweaks.f.logMaxEnergyDensity = std::log10(mma.max()));
}


//----------------+
// GLFW callbacks |
//----------------+

void mapCursorPos(GLFWwindow* window, int &x, int &y) {
  int wWidth, wHeight;
  glfwGetWindowSize(window, &wWidth, &wHeight);
  int fbWidth, fbHeight;
  glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
  x = x * fbWidth  / wWidth;
  y = y * fbHeight / wHeight;
}

void callback_framebuffersize(GLFWwindow* window, int width, int height) {
  window_width = width;
  window_height = height;
  mdataA.camera.reshape(width, height);
  mdataB.camera.reshape(width, height);
}

void callback_mousebutton(GLFWwindow* window, int button, int action, int mods) {
  auto modFlag = kt84::glfw_util::parseMods(mods);
  auto actionFlag = kt84::glfw_util::parseAction(action);
  auto mouse = kt84::glfw_util::getCursorPos(window);
  mapCursorPos(window, mouse.x, mouse.y);

  if (actionFlag.press) {
    bool clickedOnLeftSide = mouse.x < window_width / 2;
    if (modFlag.alt) {
      camera_active = clickedOnLeftSide ? &mdataA.camera : &mdataB.camera;
      camera_active->mouse_down(mouse.x, mouse.y, modFlag.shift ? Camera::DragMode::ZOOM : modFlag.ctrl ? Camera::DragMode::PAN : Camera::DragMode::ROTATE);
    } else if (modFlag.ctrl) {
    }
  } else if (actionFlag.release) {
    if (camera_active) {
      if (camera_active->drag_mode == Camera::DragMode::PAN) {
        double sx = (double)window_width * (camera_active == &mdataA.camera ? 0.25 : 0.75);
        double sy = (double)window_height * 0.5;
        double sz = (double)read_depth(sx, sy);
        if (sz < 1.) {
          setCameraMatrix(camera_active == &mdataA.camera ? mdataA : mdataB);
          camera_active->update_center(unproject(Vector3d{sx, sy, sz}));
        }
      }
      camera_active->mouse_up();
      camera_active = nullptr;
    }
  }
}

void callback_cursorpos(GLFWwindow* window, double xpos, double ypos) {
  int x = static_cast<int>(xpos);
  int y = static_cast<int>(ypos);
  mapCursorPos(window, x, y);

  if (camera_active) {
    camera_active->mouse_move(x, y);
  }
}

void callback_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
  auto actionFlag = kt84::glfw_util::parseAction(action);
  if (!actionFlag.press) return;
  if (key == GLFW_KEY_C) {
    initializeCamera(mdataA);
    initializeCamera(mdataB);
  }
}

//-----------+
// Rendering |
//-----------+

void drawViewportBorder() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glLineWidth(5);
  glColor3d(0,0,0);
  glBegin(GL_LINE_STRIP);
  glVertex2d(-1,-1);
  glVertex2d(1,-1);
  glVertex2d(1,1);
  glVertex2d(-1,1);
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void drawModel(const ModelData& mdata) {
  for (size_t i = 0; i < 8; ++i) {
    if (lightParam[i].enabled) {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0 + i);
    } else {
      glDisable(GL_LIGHT0 + i);
    }
  }

  // Faces
  if (tweaks.b.showFaces) {
    std::string fillMode = &mdata == &mdataA ? tweaks.s.fillModeA : tweaks.s.fillModeB;

    // Fall back to Constant when texCoords are unavailable
    if (fillMode == "Texture" && !mdata.texCoords)
      fillMode = "Constant";

    // Toggle texturing depending on fill mode
    if (fillMode == "Distortion" || fillMode == "Texture") {
      glEnable(GL_TEXTURE_2D);
      glColor3d(1,1,1);             // For texture modulation
      (fillMode == "Distortion" ? colormapTexture : mdata.texture).bind();
    }
    if (fillMode == "Constant") {
      glColor3f(tweaks.vec3.faceColor);
    }

    // Use GL_TRIANGLES for triangle mesh
    if (mdata.mesh->isTriangular()) {
      glBegin(GL_TRIANGLES);
      for (size_t fid = 0; fid < nF; ++fid) {
        Face fA = mdataA.mesh->face(fid);
        Face fB = mdataB.mesh->face(fid);
        Face f = &mdata == &mdataA ? fA : fB;

        if (mdataA.faceAreas[fA] == 0.) continue;
        if (mdataB.faceAreas[fB] == 0.) continue;
        if (energyPerFace[fA] == 0.) continue;

        glNormal3d(mdata.geometry->faceNormals[f]);

        if (fillMode == "Distortion") {
          // Energy density is obatined by dividing per-face energy by the average face area
          double energyDensity = energyPerFace[fA] / ((mdataA.faceAreas[fA] + mdataB.faceAreas[fB]) * 0.5);
          double t = std::log10(energyDensity / 4.) / (tweaks.f.logMaxEnergyDensity - std::log10(4.));
          CIT_ASSERT(t >= -1.e-10);
          t = std::max<double>(std::min<double>(t, 1.), 0.);
          glTexCoord2d(t, 0.5);
        }

        for (Corner c : f.adjacentCorners()) {
          if (fillMode == "Texture")
            glTexCoord2d(tweaks.f.textureScale * (*mdata.texCoords)[c]);
          glVertex3d(mdata.geometry->inputVertexPositions[c.vertex()]);
        }
      }
      glEnd();

    // Use GL_POLYGON for polygonal mesh
    } else {
      for (size_t fid = 0; fid < nF; ++fid) {
        Face fA = mdataA.mesh->face(fid);
        Face fB = mdataB.mesh->face(fid);
        Face f = &mdata == &mdataA ? fA : fB;

        if (mdataA.faceAreas[fA] == 0.) continue;
        if (mdataB.faceAreas[fB] == 0.) continue;
        if (energyPerFace[fA] == 0.) continue;

        glNormal3d(mdata.geometry->faceNormals[f]);

        if (fillMode == "Distortion") {
          // Energy density is obatined by dividing per-face energy by the average face area
          double energyDensity = energyPerFace[fA] / ((mdataA.faceAreas[fA] + mdataB.faceAreas[fB]) * 0.5);
          double t = std::log10(energyDensity / 4.) / (tweaks.f.logMaxEnergyDensity - std::log10(4.));
          CIT_ASSERT(t >= 0.);
          t = std::min<double>(t, 1.);
          glTexCoord2d(t, 0.5);
        }

        glBegin(GL_POLYGON);
        for (Corner c : (&mdata == &mdataA ? fA : fB).adjacentCorners()) {
          if (fillMode == "Texture")
            glTexCoord2d(tweaks.f.textureScale * (*mdata.texCoords)[c]);
          glVertex3d(mdata.geometry->inputVertexPositions[c.vertex()]);
        }
        glEnd();
      }
    }
  }

  // Drawing with no lighting and texturing
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);

  // Edges
  if (tweaks.b.showEdges) {
    glLineWidth(tweaks.f.edgeWidth);
    glColor3d(Vector3::constant(0.2));
    glBegin(GL_LINES);
    for (Edge e : mdata.mesh->edges()) {
      for (Vertex v : e.adjacentVertices())
        glVertex3d(mdata.geometry->inputVertexPositions[v]);
    }
    glEnd();
  }

  // Bounding box
  if (tweaks.b.showBoundingBox) {
    glLineWidth(5);
    glBegin(GL_LINES);
    using CornerType = AlignedBox3d::CornerType;
    glColor3d(1,0,0); glVertex3d(mdata.bbox.min()); glVertex3d(mdata.bbox.corner(CornerType::BottomRight));
    glColor3d(0,1,0); glVertex3d(mdata.bbox.min()); glVertex3d(mdata.bbox.corner(CornerType::TopLeft));
    glColor3d(0,0,1); glVertex3d(mdata.bbox.min()); glVertex3d(mdata.bbox.corner(CornerType::BottomLeftCeil));
    glEnd();
  }
}

void draw() {
  for (size_t i = 0; i < 8; ++i) {
    glLightAmbient3f (Vector3::constant(lightParam[i].ambient ), i);
    glLightDiffuse3f (Vector3::constant(lightParam[i].diffuse ), i);
    glLightSpecular3f(Vector3::constant(lightParam[i].specular), i);
  }
  glMaterialAmbient3f(Vector3::constant(materialParam.ambient));
  glMaterialSpecular3f(Vector3::constant(materialParam.specular));
  glMaterialShininess1f(materialParam.shininess);

  setCameraMatrix(mdataA);
  drawViewportBorder();
  mdataA.dispList.render([]() {
    drawModel(mdataA);
  });

  setCameraMatrix(mdataB);
  drawViewportBorder();
  mdataB.dispList.render([]() {
    drawModel(mdataB);
  });

  ImGui::Begin("Draw options (om_viewer)");

  if (ImGui::Checkbox("Show faces", &tweaks.b.showFaces)) invalidateDispList();
  if (ImGui::Checkbox("Show edges", &tweaks.b.showEdges)) invalidateDispList();
  if (ImGui::Checkbox("Show bounding box", &tweaks.b.showBoundingBox)) invalidateDispList();

  if (ImGui::TreeNode("Light parameters")) {
    static int lightIndex = 0;
    ImGui::SliderInt("Index", &lightIndex, 0, 7);
    if (ImGui::Checkbox("Enabled", &lightParam[lightIndex].enabled)) invalidateDispList();
    if (ImGui::SliderFloat("Ambient", &lightParam[lightIndex].ambient, 0.f, 1.f)) invalidateDispList();
    if (ImGui::SliderFloat("Diffuse", &lightParam[lightIndex].diffuse, 0.f, 1.f)) invalidateDispList();
    if (ImGui::SliderFloat("Specular", &lightParam[lightIndex].specular, 0.f, 1.f)) invalidateDispList();
    if (ImGui::gizmo3D("Direction", *(vec3*)&lightParam[lightIndex].position[0])) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Material parameters")) {
    if (ImGui::SliderFloat("Ambient", &materialParam.ambient, 0.f, 1.f)) invalidateDispList();
    if (ImGui::SliderFloat("Specular", &materialParam.specular, 0.f, 1.f)) invalidateDispList();
    if (ImGui::SliderFloat("Shininess", &materialParam.shininess, 0.f, 128.f)) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::ColorEdit3("Face color", &tweaks.vec3.faceColor[0])) invalidateDispList();

  if (ImGui::SliderFloat("Edge width", &tweaks.f.edgeWidth, 0.01f, 50.f)) invalidateDispList();
  if (ImGui::SliderFloat("Texture scale", &tweaks.f.textureScale, 0.f, 10.f)) invalidateDispList();
  if (ImGui::SliderFloat("Log max energy density", &tweaks.f.logMaxEnergyDensity, std::log10(4.1f), 4.f, "%.3f", 5.f)) invalidateDispList();

  ImGui::Separator();

  if (ImGui::BeginCombo("Fill mode A", tweaks.s.fillModeA.c_str())) {
    for (size_t i = 0; i < fillModes.size(); i++) {
      if (ImGui::Selectable(fillModes[i].c_str(), tweaks.s.fillModeA == fillModes[i])) {
        tweaks.s.fillModeA = fillModes[i];
        mdataA.dispList.invalidate();
      }
    }
    ImGui::EndCombo();
  }
  if (ImGui::BeginCombo("Fill mode B", tweaks.s.fillModeB.c_str())) {
    for (size_t i = 0; i < fillModes.size(); i++) {
      if (ImGui::Selectable(fillModes[i].c_str(), tweaks.s.fillModeB == fillModes[i])) {
        tweaks.s.fillModeB = fillModes[i];
        mdataB.dispList.invalidate();
      }
    }
    ImGui::EndCombo();
  }

  ImGui::End();
}

//------+
// Main |
//------+

int main(int argc, char** argv) {
  // Configure the argument parser
  // clang-format off
  args::ArgumentParser parser("Overlay Mesh Viewer");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<std::string> argTextureA(parser, "file", "Texture image file for A", {"textureA"});
  args::ValueFlag<std::string> argTextureB(parser, "file", "Texture image file for B", {"textureB"});
  std::string logLevelHelp;
  {
    std::ostringstream oss;
    oss << "Log level {";
    bool first = true;
    for (auto s : SPDLOG_LEVEL_NAMES) {
      if (!first)
        oss << ", ";
      first = false;
      oss << s;
    }
    oss << "}";
    logLevelHelp = oss.str();
  }
  args::ValueFlag<std::string> argLogLevel(parser, "level", logLevelHelp, {"log-level"});

  args::Positional<std::string> argMeshOnA(parser, "meshOnA", "Overlay mesh on model A");
  args::Positional<std::string> argMeshOnB(parser, "meshOnB", "Overlay mesh on model B");

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
  if (!argMeshOnB || !argMeshOnA) {
    std::cout << parser;
    return EXIT_FAILURE;
  }

  // Configure logging
  spdlog::set_pattern("%D %T %L %s:%# %v");
  if (argLogLevel) {
    int i = 0;
    for (auto s : SPDLOG_LEVEL_NAMES) {
      if (args::get(argLogLevel) == s) {
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
  spdlog::default_logger()->sinks().push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("om_viewer.log", true));

  mdataA.shortName = "A";
  mdataB.shortName = "B";

  readMesh(mdataA, args::get(argMeshOnA));
  readMesh(mdataB, args::get(argMeshOnB));

  checkTopologicalConsistency();

  computeCornerPositionsInFace(mdataA);
  computeCornerPositionsInFace(mdataB);

  computeFaceAreas(mdataA);
  computeFaceAreas(mdataB);

  normalizeSurfaceArea(mdataA);
  normalizeSurfaceArea(mdataB);

  computeJacobian();
  computeEnergy();

  // glewInit();

  // Initialize GLFW
  kt84::glfw_util::InitConfig initConfig;
  initConfig.window.title = "Overlay Mesh Viewer";
  initConfig.callback.framebuffersize = callback_framebuffersize;
  initConfig.callback.mousebutton = callback_mousebutton;
  initConfig.callback.cursorpos = callback_cursorpos;
  initConfig.callback.key = callback_key;
  kt84::glfw_util::LoopControl loopControl = kt84::glfw_util::init(initConfig);
  // Idle functions
  loopControl.idle_func[0] = []() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
  };
  loopControl.idle_func[1] = draw;
  loopControl.idle_func[2] = []() {
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
  };

  // OpenGL initial settings
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonOffset(1,1);
  glClearColor(1,1,1,1);

  lightParam[0].enabled = true;
  lightParam[0].diffuse = 1.f;
  lightParam[0].specular = 1.f;

  initializeCamera(mdataA);
  initializeCamera(mdataB);

  // Setup model texture
  stbi_set_flip_vertically_on_load(1);
  for (int i = 0; i < 2; ++i) {
    args::ValueFlag<std::string>* argTexture;
    ModelData* mdata;
    uint8_t* checkerboard_jpg;
    uint32_t* checkerboard_jpg_size;
    if (i == 0) {
      argTexture = &argTextureA;
      mdata = &mdataA;
      checkerboard_jpg = checkerboarda_jpg;
      checkerboard_jpg_size = &checkerboarda_jpg_size;
    } else {
      argTexture = &argTextureB;
      mdata = &mdataB;
      checkerboard_jpg = checkerboardb_jpg;
      checkerboard_jpg_size = &checkerboardb_jpg_size;
    }

    int x,y,n;
    unsigned char *data = *argTexture ?
      stbi_load(args::get(*argTexture).c_str(), &x, &y, &n, 4) :
      stbi_load_from_memory(checkerboard_jpg, *checkerboard_jpg_size, &x, &y, &n, 4);

    if (!data) {
      CIT_ASSERT(*argTexture);
      LOG_ERROR("Failed to load texture: {}", args::get(*argTexture));
      return EXIT_FAILURE;
    }

    mdata->texture.init();
    mdata->texture.bind();
    mdata->texture.set_default_param();
    mdata->texture.set_env(GL_MODULATE);
    mdata->texture.allocate(x, y);
    mdata->texture.copy_cpu2gpu(GL_UNSIGNED_BYTE, data);

    stbi_image_free(data);
  }

  // Setup colormap texture
  {
    int x,y,n;
    unsigned char *data = stbi_load_from_memory(colormap_jpg, colormap_jpg_size, &x, &y, &n, 4);
    if (!data) {
      LOG_ERROR("Failed to load colormap texture");
      return EXIT_FAILURE;
    }

    colormapTexture.init();
    colormapTexture.bind();
    colormapTexture.set_default_param();
    colormapTexture.set_env(GL_MODULATE);
    colormapTexture.set_wrap(GL_CLAMP_TO_EDGE);
    colormapTexture.allocate(x, y);
    colormapTexture.copy_cpu2gpu(GL_UNSIGNED_BYTE, data);

    stbi_image_free(data);
  }

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(loopControl.window, true);
  ImGui_ImplOpenGL2_Init();

  kt84::glfw_util::start_loop(loopControl);

  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  return EXIT_SUCCESS;
}
