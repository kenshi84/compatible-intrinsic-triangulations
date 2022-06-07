#include "cit.hpp"
#include "imgui.h"

#include "imgui/backends/imgui_impl_glfw.h"
#include "imgui/backends/imgui_impl_opengl2.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreorder"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#include "imGuIZMOquat.h"
#pragma clang diagnostic pop

#include "kt84/graphics/graphics_util.hh"
#include "kt84/glut_clone/geometry.hh"
using namespace kt84::graphics_util;

extern VectorXd dbg_descentDirection;
extern std::deque<VectorXsp> dbg_previousX;
extern std::deque<VectorXd> dbg_previousGradient;
extern std::deque<VertexData<std::vector<SurfacePoint>>> dbg_vertexPathA;
extern std::deque<VertexData<std::vector<SurfacePoint>>> dbg_vertexPathB;

namespace {

int initTweaks = [](){
  cit::tweaks.vec3["inputFaceColor"] = Vector3f::Constant(210.f) / 255.f;
  cit::tweaks.b["showAxis"] = true;
  cit::tweaks.b["showFaces"] = true;
  cit::tweaks.b["highQuality"] = false;
  cit::tweaks.i["seed"] = 0;

  // Vertex
  cit::tweaks.f["vertexRadius"];// = 0.003f;
  cit::tweaks.f["vertexSilhouetteRatio"] = 0.3f;
  cit::tweaks.f["vertexSilhouetteOffset"] = 0.3f;

  // Input edge
  cit::tweaks.f["inputEdgeRadius"];// = 0.002f;
  cit::tweaks.f["inputEdgeOffset"] = 0.f;

  // Intrinsic edge
  cit::tweaks.f["compatibleEdgeRadius"];// = 0.002f;
  cit::tweaks.f["incompatibleEdgeRadius"];// = 0.0025f;
  cit::tweaks.f["intrinsicEdgeOffset"] = 0.f;
  cit::tweaks.f["incompatibleEdgeBrightness"] = 0.4f;
  cit::tweaks.b["showIncompatiblePatchBoundary"] = true;

  // Mapped edge
  cit::tweaks.f["mappedEdgeRadius"];// = 0.001f;
  cit::tweaks.f["mappedEdgeOffset"] = 0.f;

  cit::tweaks.f["highlightedEdgeWidth"] = 20.f;
  cit::tweaks.f["vectorFieldLineWidth"];// = 15.f;
  cit::tweaks.f["vectorFieldLineScale"] = 10.f;
  cit::tweaks.f["logMaxEnergyDensity"] = 3.f;
  cit::tweaks.b["descentDirection"] = true;
  cit::tweaks.b["previousGradient"] = true;
  cit::tweaks.b["previousGradientAll"] = true;
  cit::tweaks.i["previousGradientSelected"] = 0;
  cit::tweaks.b["transportedGradient"] = true;
  cit::tweaks.b["smoothedGradient"] = true;

  // Texture
  cit::tweaks.f["textureScale"] = 1.f;
  cit::tweaks.vec3["textureOffset"] = {0.f, 0.f, 0.f};

  cit::tweaks.s["fillModeA"] = "Constant";
  cit::tweaks.s["fillModeB"] = "Constant";
  return 0;
}();
const std::array<std::string, 5> fillModes = {"TextureA", "TextureB", "Distortion", "fIntrinsic", "Constant"};
}

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

namespace cit {

void invalidateDispList() {
  mdataA.dispList.invalidate();
  mdataB.dispList.invalidate();
}

void dumpLightMaterialParams() {
  DLOG_INFO(0, "BEGIN dumpLightMaterialParams");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_INFO(0, "END   dumpLightMaterialParams"); });

  for (size_t i = 0; i < 8; ++i) {
    if (lightParam[i].enabled) {
      DLOG_INFO(1, "lightParam[{}].enabled = true;", i);
      DLOG_INFO(1, "lightParam[{}].ambient = {};", i, lightParam[i].ambient);
      DLOG_INFO(1, "lightParam[{}].diffuse = {};", i, lightParam[i].diffuse);
      DLOG_INFO(1, "lightParam[{}].specular = {};", i, lightParam[i].specular);
      DLOG_INFO(1, "lightParam[{}].position = Vector4f{{ {}, {}, {}, {} }};", i, lightParam[i].position[0], lightParam[i].position[1], lightParam[i].position[2], lightParam[i].position[3]);
    }
  }
  DLOG_INFO(1, "materialParam.ambient = {};"  , materialParam.ambient);
  DLOG_INFO(1, "materialParam.specular = {};" , materialParam.specular);
  DLOG_INFO(1, "materialParam.shininess = {};", materialParam.shininess);
}

void drawCoInTri(const CoInTri& cointriP, const CoInTri& cointriQ) {
  for (size_t i = 0; i < 8; ++i) {
    if (lightParam[i].enabled) {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0 + i);
    } else {
      glDisable(GL_LIGHT0 + i);
    }
  }

  const ModelData& mdataP = *cointriP.mdata;
  const ModelData& mdataQ = *cointriQ.mdata;

  // input faces
  if (tweaks.b["showFaces"]) {
    std::string fillMode = tweaks.s["fillMode" + mdataP.name];
    if (fillMode == "Texture" + mdataP.name && !mdataP.texCoords) fillMode = "Constant";
    if (fillMode == "Texture" + mdataQ.name && !mdataQ.texCoords) fillMode = "Constant";
    if (fillMode == "Constant" || fillMode == "fIntrinsic") {
      // Constant color, no texturing
      if (fillMode == "Constant") {
        glColor3f(tweaks.vec3["inputFaceColor"]);
        glBegin(GL_TRIANGLES);
        for (Face f : mdataP.mesh->faces()) {
          glNormal3d(mdataP.geometry->faceNormals[f]);
          for (Vertex v : f.adjacentVertices()) {
            glVertex3d(mdataP.geometry->vertexPositions[v]);
          }
        }
        glEnd();

      } else {
        for (const OverlayPolygon& polygon : currentConfig.overlayPolygons) {
          glColor3d(Vector3::constant(0.3) + 0.7 * getRandomColor(polygon.fIntrinsic.getIndex()));
          glBegin(GL_POLYGON);
          for (const OverlayWedge& wedge : polygon.wedges) {
            SurfacePoint facePointP = mdataP.name == "A" ? wedge.facePointA : wedge.facePointB;
            glNormal3d(mdataP.geometry->faceNormals[facePointP.face]);
            glVertex3d(facePointP.interpolate(mdataP.geometry->inputVertexPositions));
          }
          glEnd();
        }
      }

    } else {
      // Textured drawing
      glEnable(GL_TEXTURE_2D);
      if (fillMode == "Texture" + mdataP.name) {
        // Model P's original texture
        glColor3d(1,1,1);             // For texture modulation
        mdataP.texture.bind();
        glBegin(GL_TRIANGLES);
        for (Face f : mdataP.mesh->faces()) {
          glNormal3d(mdataP.geometry->faceNormals[f]);
          for (Corner c : f.adjacentCorners()) {
            glTexCoord2d(tweaks.f["textureScale"] * (*mdataP.texCoords)[c] + Vector2::castFrom(tweaks.vec3["textureOffset"]));
            glVertex3d(mdataP.geometry->vertexPositions[c.vertex()]);
          }
        }
        glEnd();

      } else if (fillMode == "Texture" + mdataQ.name && currentConfig.topologyValid) {
        // Model Q's texture transferred on P
        glColor3d(1,1,1);             // For texture modulation
        mdataQ.texture.bind();
        for (const OverlayPolygon& polygon : currentConfig.overlayPolygons) {
          glBegin(GL_POLYGON);
          for (const OverlayWedge& wedge : polygon.wedges) {
            SurfacePoint facePointP = mdataP.name == "A" ? wedge.facePointA : wedge.facePointB;
            SurfacePoint facePointQ = mdataP.name == "A" ? wedge.facePointB : wedge.facePointA;

            // Get interpolated texture coordinate at the wedge
            Halfedge Q_he = facePointQ.face.halfedge();
            Vector2 interpolated_texCoord{0., 0.};
            for (int i = 0; i < 3; ++i) {
              interpolated_texCoord += facePointQ.faceCoords[i] * (*mdataQ.texCoords)[Q_he.corner()];
              Q_he = Q_he.next();
            }
            CIT_ASSERT(Q_he == facePointQ.face.halfedge());

            glNormal3d(mdataP.geometry->faceNormals[facePointP.face]);
            glTexCoord2d(tweaks.f["textureScale"] * interpolated_texCoord + Vector2::castFrom(tweaks.vec3["textureOffset"]));
            glVertex3d(facePointP.interpolate(mdataP.geometry->inputVertexPositions));
          }
          glEnd();
        }

      } else if (fillMode == "Distortion") {
        // Distortion color map
        glColor3d(1,1,1);             // For texture modulation
        colormapTexture.bind();
        for (const OverlayPolygon& polygon : currentConfig.overlayPolygons) {
          double t = std::log(currentConfig.energyDensity[polygon.fIntrinsic] / 4.) / (tweaks.f["logMaxEnergyDensity"] - std::log(4.));
          CIT_ASSERT(t >= -1.e-10);
          t = std::max<double>(std::min<double>(t, 1.), 0.);
          glTexCoord2d(t, 0.5);
          glBegin(GL_POLYGON);
          for (const OverlayWedge& wedge : polygon.wedges) {
            SurfacePoint facePointP = mdataP.name == "A" ? wedge.facePointA : wedge.facePointB;
            glNormal3d(mdataP.geometry->faceNormals[facePointP.face]);
            glVertex3d(facePointP.interpolate(mdataP.geometry->inputVertexPositions));
          }
          glEnd();
        }
      }

      glDisable(GL_TEXTURE_2D);
    }
  }

  // vertices
  if (tweaks.f["vertexRadius"]) {
    glLineWidth(1);
    for (Vertex v : cointriP.signpostTri->intrinsicMesh->vertices()) {
      if (tweaks.b["vertexHideOriginal"] && v.getIndex() < mdataP.nV) continue;
      if (tweaks.b["vertexHideInserted"] && v.getIndex() >= mdataP.nV) continue;

      glPushMatrix();

      if (tweaks.b["vertexUseModelColor"]) {
        if (v.getIndex() < mdataP.nV)
          glColor3d(mdataP.inputEdgeColor);
        else
          glColor3d(mdataQ.inputEdgeColor);
      } else {
        glColor3d(getRandomColor(cointriP.uniqueID_per_Vertex[v]));
      }

      // Translate to the vertex position
      const SurfacePoint& p = cointriP.signpostTri->vertexLocations[v];
      Vector3 xyz = p.interpolate(mdataP.geometry->inputVertexPositions);
      glTranslated(xyz);

      // Determine actual radius
      double r = tweaks.f["vertexRadius"];
      Vector3 toEye = Vector3::castFrom(mdataP.camera.eye) - xyz;
      r *= norm(toEye);
      if (isVertexMerged(cointriP.uniqueID_per_Vertex[v]))
        r *= 2;

      const int N = tweaks.b["highQuality"] ? 30 : 6;
      kt84::glutSolidSphere(r, N, N / 2);

      // Draw silhouette if original
      if (v.getIndex() < mdataP.nV && tweaks.f["vertexSilhouetteRatio"]) {
        // Align Z axis to eye direction
        Vector3 unitZ{0,0,1};
        glRotated(angle(unitZ, toEye) * 180. / M_PI, cross(unitZ, toEye));

        // Slight offset along Z
        glTranslated(0, 0, r * tweaks.f["vertexSilhouetteOffset"]);

        // Draw disk
        glColor3d(0,0,0);
        glDisable(GL_LIGHTING);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2d(0,0);
        for (int i = 0; i <= N; ++i)
          glVertex2d(r * (1. + tweaks.f["vertexSilhouetteRatio"]) * Vector2::fromAngle(2. * M_PI * i / (double)N));
        glEnd();
        glEnable(GL_LIGHTING);
      }
      glPopMatrix();
    }
  }

  // input edges
  if (tweaks.f["inputEdgeRadius"]) {
    glColor3f(mdataP.inputEdgeColor);
    for (Edge e : mdataP.mesh->edges()) {
      glPushMatrix();

      Vector3 p0 = mdataP.geometry->vertexPositions[e.halfedge().vertex()];
      Vector3 p1 = mdataP.geometry->vertexPositions[e.halfedge().tipVertex()];

      // Determine actual radius
      double r = tweaks.f["inputEdgeRadius"];
      Vector3 toEye = Vector3::castFrom(mdataP.camera.eye) - 0.5 * (p0 + p1);
      r *= norm(toEye);

      // Translate
      glTranslated(p0 + tweaks.f["inputEdgeOffset"] * r * normalize(toEye));

      // Rotate
      Vector3 unitZ{0,0,1};
      glRotated(angle(unitZ, p1 - p0) * 180. / M_PI, cross(unitZ, p1 - p0));

      kt84::glutSolidCylinder(r, mdataP.geometry->edgeLengths[e], 6, 1);

      glPopMatrix();
    }
  }

  // intrinsic edges
  for (Edge e : cointriP.signpostTri->intrinsicMesh->edges()) {
    double r;
    Vector3 c;
    if (isEdgeCompatible(cointriP, e)) {
      r = tweaks.f["compatibleEdgeRadius"];
      if (tweaks.b["drawEdgeInBlack"])
        c = {0., 0., 0.};
      else if (tweaks.b["showIncompatiblePatchBoundary"] && isFaceCompatible(cointriP, e.halfedge().face()) + isFaceCompatible(cointriP, e.halfedge().twin().face()) == 1)
        c = {1.,0.,0.};
      else
        c = Vector3::constant(0.5) + 0.5 * getRandomColor(getUniqueEdgeID(cointriP, e));

    } else {
      r = tweaks.f["incompatibleEdgeRadius"];
      c = Vector3::constant(tweaks.f["incompatibleEdgeBrightness"]);
    }

    if (!r) continue;
    glColor3d(c);

    for (size_t i = 0; i < cointriP.intrinsicEdgePath[e].size() - 1; ++i) {
      glPushMatrix();

      Vector3 p0 = cointriP.intrinsicEdgePath[e][i    ].interpolate(mdataP.geometry->inputVertexPositions);
      Vector3 p1 = cointriP.intrinsicEdgePath[e][i + 1].interpolate(mdataP.geometry->inputVertexPositions);

      // Determine actual radius
      double r2 = r;
      Vector3 toEye = Vector3::castFrom(mdataP.camera.eye) - 0.5 * (p0 + p1);
      r2 *= norm(toEye);

      // Translate
      glTranslated(p0 + tweaks.f["intrinsicEdgeOffset"] * r2 * normalize(toEye));

      // Rotate
      Vector3 unitZ{0,0,1};
      glRotated(angle(unitZ, p1 - p0) * 180. / M_PI, cross(unitZ, p1 - p0));

      kt84::glutSolidCylinder(r2, norm(p1 - p0), 6, 1);

      glPopMatrix();
    }
  }

  // Q's input edges mapped onto P
  if (currentConfig.topologyValid && tweaks.f["mappedEdgeRadius"]) {
    glColor3f(mdataQ.inputEdgeColor);
    for (Edge Q_e : mdataQ.mesh->edges()) {
      for (size_t i = 0; i < cointriQ.inputEdgePathOnOtherInput[Q_e].size() - 1; ++i) {
        glPushMatrix();

        Vector3 p0 = cointriQ.inputEdgePathOnOtherInput[Q_e][i    ].interpolate(mdataP.geometry->inputVertexPositions);
        Vector3 p1 = cointriQ.inputEdgePathOnOtherInput[Q_e][i + 1].interpolate(mdataP.geometry->inputVertexPositions);

        // Determine actual radius
        double r = tweaks.f["mappedEdgeRadius"];
        Vector3 toEye = Vector3::castFrom(mdataP.camera.eye) - 0.5 * (p0 + p1);
        r *= norm(toEye);

        // Translate
        glTranslated(p0 + tweaks.f["mappedEdgeOffset"] * r * normalize(toEye));

        // Rotate
        Vector3 unitZ{0,0,1};
        glRotated(angle(unitZ, p1 - p0) * 180. / M_PI, cross(unitZ, p1 - p0));

        kt84::glutSolidCylinder(r, norm(p1 - p0), 6, 1);

        glPopMatrix();
      }
    }
  }

  // Drawing with no lighting
  glDisable(GL_LIGHTING);

  // bounding box
  if (tweaks.b["showAxis"]) {
    glLineWidth(5);
    glBegin(GL_LINES);
    using CornerType = AlignedBox3d::CornerType;
    glColor3d(1,0,0); glVertex3d(mdataP.bbox.min()); glVertex3d(mdataP.bbox.corner(CornerType::BottomRight));
    glColor3d(0,1,0); glVertex3d(mdataP.bbox.min()); glVertex3d(mdataP.bbox.corner(CornerType::TopLeft));
    glColor3d(0,0,1); glVertex3d(mdataP.bbox.min()); glVertex3d(mdataP.bbox.corner(CornerType::BottomLeftCeil));
    glEnd();
  }

  // edges around highlighted vertex
  if (tweaks.f["highlightedEdgeWidth"]) {
    glLineWidth(tweaks.f["highlightedEdgeWidth"]);
    if (cointriP.uniqueID_to_Vertex.count(highlightedVertexID)) {
      Vertex v = cointriP.uniqueID_to_Vertex.at(highlightedVertexID);
      for (Edge e : v.adjacentEdges()) {
        if (isEdgeCompatible(cointriP, e)) {
          glColor3d(Vector3::constant(0.5) + 0.5 * getRandomColor(getUniqueEdgeID(cointriP, e)));
        } else {
          glColor3d(0,0,0);
        }
        glBegin(GL_LINE_STRIP);
        for (const SurfacePoint& p : cointriP.intrinsicEdgePath[e]) {
          glVertex3d(p.interpolate(mdataP.geometry->inputVertexPositions));
        }
        glEnd();
      }
    }

    // highlighted intrinsic edge
    Edge highlightedEdge = getEdgeByUniqueID(cointriP, highlightedEdgeID);
    if (highlightedEdge != Edge()) {
      if (isEdgeCompatible(cointriP, highlightedEdge)) {
        glColor3d(Vector3::constant(0.5) + 0.5 * getRandomColor(highlightedEdgeID));
      } else {
        glColor3d(0,0,0);
      }
      glBegin(GL_LINE_STRIP);
      for (const SurfacePoint& p : cointriP.intrinsicEdgePath[highlightedEdge]) {
        glVertex3d(p.interpolate(mdataP.geometry->inputVertexPositions));
      }
      glEnd();
    }

    // highlighted intrinsic face
    if (0 <= highlightedFaceID && highlightedFaceID < (int)cointriP.signpostTri->intrinsicMesh->nFaces()) {
      Face f = cointriP.signpostTri->intrinsicMesh->face(highlightedFaceID);
      for (Edge e : f.adjacentEdges()) {
        if (isEdgeCompatible(cointriP, e)) {
          glColor3d(Vector3::constant(0.5) + 0.5 * getRandomColor(getUniqueEdgeID(cointriP, e)));
        } else {
          glColor3d(0,0,0);
        }
        glBegin(GL_LINE_STRIP);
        for (const SurfacePoint& p : cointriP.intrinsicEdgePath[e]) {
          glVertex3d(p.interpolate(mdataP.geometry->inputVertexPositions));
        }
        glEnd();
      }
    }
  }

#if 0
  // highlighted input edge
  glLineWidth(tweaks.f["inputEdgeRadius"] + 2);
  glColor3f(1,0,0);
  for (Edge e : mdataP.mesh->edges()) {
    Vertex v0 = e.halfedge().vertex();
    Vertex v1 = e.halfedge().twin().vertex();
    if (highlightedEdgeID != std::pair<int,int>(v0.getIndex(), v1.getIndex())) continue;
    glBegin(GL_LINES);
    glVertex3d(mdataP.geometry->vertexPositions[v0]);
    glVertex3d(mdataP.geometry->vertexPositions[v1]);
    glEnd();
    break;
  }
#endif

#if 0
  if (tweaks.f["highlightedEdgeWidth"] && sysParam.angleThreshold) {
    glLineWidth(tweaks.f["highlightedEdgeWidth"]);
    glColor3f(1,0,0);
    for (Halfedge he1 : cointriP.signpostTri->intrinsicMesh->halfedges()) {
      if (cointriP.signpostTri->cornerAngle(he1.corner()) < sysParam.angleThreshold) {
        Halfedge he0 = he1.next().next();
        for (Halfedge he : {he0, he1}) {
          glBegin(GL_LINE_STRIP);
          for (const SurfacePoint& p : cointriP.intrinsicEdgePath[he.edge()])
            glVertex3d(p.interpolate(mdataP.geometry->inputVertexPositions));
          glEnd();
        }
      }
    }
  }
#endif

  // The rest visualizes derivative-related data
  if (!currentConfig.derivativeComputed)
    return;

  // descent vectors
  if (tweaks.f["vectorFieldLineWidth"] && tweaks.b["descentDirection"]) {
    glLineWidth(tweaks.f["vectorFieldLineWidth"]);
    glBegin(GL_LINES);
    for (size_t i = 0; i < mdataQ.nV; ++i) {
      Vertex Q_v = cointriQ.signpostTri->intrinsicMesh->vertex(i);
      Vertex P_v = cointriQ.correspondingVertex[Q_v];
      SurfacePoint p_sp = cointriP.signpostTri->vertexLocations[P_v];

      if (p_sp.type == SurfacePointType::Vertex) continue;

      if (p_sp.type == SurfacePointType::Edge)
        p_sp = convertEdgePointToFacePoint(p_sp);

      Vector3 p_xyz = p_sp.interpolate(mdataP.geometry->inputVertexPositions);

      // extract barycentric coordinate delta
      SurfacePoint delta_sp{p_sp.face, {}};
      size_t idx = i;
      if (&mdataP == &mdataA)
        idx += mdataA.nV;
      for (size_t di = 0; di < 2; ++di) {
        delta_sp.faceCoords[di] = dbg_descentDirection[2*idx + di];
      }
      delta_sp.faceCoords[2] = -(delta_sp.faceCoords[0] + delta_sp.faceCoords[1]);
      Vector3 delta_xyz = delta_sp.interpolate(mdataP.geometry->inputVertexPositions);
      delta_xyz *= tweaks.f["vectorFieldLineScale"];
      glColor3d(0,0,1);
      glVertex3d(p_xyz);
      glColor3d(1,1,0);
      glVertex3d(p_xyz + delta_xyz);
    }
    glEnd();
  }

  const std::array<Vector3, 6> gradientColors = {
    Vector3{1,0,0},
    Vector3{1,1,0},
    Vector3{0,1,0},
    Vector3{0,1,1},
    Vector3{0,0,1},
    Vector3{1,0,1}
  };

  if (tweaks.f["vectorFieldLineWidth"] && tweaks.b["previousGradient"] && !dbg_previousX.empty()) {
    // vertices in the previous configuration
    const size_t offset = &mdataP == &mdataA ? mdataA.nV : 0;
    for (size_t i = 0; i < mdataQ.nV; ++i) {
      Vertex Q_v = cointriQ.signpostTri->intrinsicMesh->vertex(i);
      Vector3 c = getRandomColor(cointriQ.uniqueID_per_Vertex[Q_v]);
      glColor3d(c);
      for (size_t k = 0; k < dbg_previousX.size(); ++k) {
        const SurfacePoint& p = dbg_previousX[k][offset + i];
        Vector3 xyz = p.interpolate(mdataP.geometry->inputVertexPositions);
        Vector3 toEye = Vector3::castFrom(mdataP.camera.eye) - xyz;
        double s = norm(toEye);
        glPushMatrix();
        glTranslated(xyz);
        kt84::glutSolidSphere(tweaks.f["vertexRadius"] * s, 6, 6);
        glPopMatrix();
      }
    }

    // previous gradient vectors
    glLineWidth(tweaks.f["vectorFieldLineWidth"]);
    glBegin(GL_LINES);
    for (size_t i = 0; i < mdataQ.nV; ++i) {
      for (size_t k = 0; k <= dbg_previousX.size(); ++k) {
        if (!tweaks.b["previousGradientAll"] && tweaks.i["previousGradientSelected"] != (int)k)
          continue;

        SurfacePoint p_sp = (k == 0 ? currentConfig.x : dbg_previousX[k - 1])[offset + i];

        if (p_sp.type == SurfacePointType::Vertex) continue;

        if (p_sp.type == SurfacePointType::Edge)
          p_sp = convertEdgePointToFacePoint(p_sp);

        Vector3 p_xyz = p_sp.interpolate(mdataP.geometry->inputVertexPositions);

        // extract barycentric coordinate delta
        SurfacePoint delta_sp{p_sp.face, {}};
        size_t idx = i;
        if (&mdataP == &mdataA)
          idx += mdataA.nV;
        for (size_t di = 0; di < 2; ++di) {
          delta_sp.faceCoords[di] = -(k == 0 ? currentConfig.gradient : dbg_previousGradient[k - 1])[2*idx + di];
        }
        delta_sp.faceCoords[2] = -(delta_sp.faceCoords[0] + delta_sp.faceCoords[1]);
        Vector3 delta_xyz = delta_sp.interpolate(mdataP.geometry->inputVertexPositions);
        delta_xyz *= tweaks.f["vectorFieldLineScale"];
        glColor3d(0.5, 0.5, 0.5);
        glVertex3d(p_xyz);
        glColor3d(gradientColors[k]);
        glVertex3d(p_xyz + delta_xyz);
      }
    }

    // vertex path
    const auto& dbg_vertexPathQ = &mdataP == &mdataA ? dbg_vertexPathB : dbg_vertexPathA;
    for (Vertex Q_v : mdataQ.mesh->vertices()) {
      glLineWidth(3);
      for (size_t k = 0; k < dbg_previousX.size(); ++k) {
        auto& vertexPath = k == 0 ? cointriQ.vertexPath : dbg_vertexPathQ[k - 1];
        if (!vertexPath.getMesh()) {
          continue;
        }
        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < vertexPath[Q_v].size(); ++i) {
          const SurfacePoint& P_sp = vertexPath[Q_v][i];
          glColor3d(0,0,i / (vertexPath[Q_v].size() - 1));
          glVertex3d(P_sp.interpolate(mdataP.geometry->inputVertexPositions));
        }
        glEnd();
      }
    }
    glEnd();
  }

  if (tweaks.f["vectorFieldLineWidth"] && tweaks.b["transportedGradient"] && !currentHistory.empty()) {
    glLineWidth(tweaks.f["vectorFieldLineWidth"]);
    glBegin(GL_LINES);
    const size_t offset = &mdataP == &mdataA ? mdataA.nV : 0;
    std::vector<const VectorXd*> transportedGradient = { &currentConfig.gradient };
    for (const auto& history : currentHistory) {
      transportedGradient.push_back(&history.first);
    }
    for (size_t k = 0; k < transportedGradient.size(); ++k) {
      if (!tweaks.b["previousGradientAll"] && tweaks.i["previousGradientSelected"] != (int)k)
        continue;

      for (size_t i = 0; i < mdataQ.nV; ++i) {
        Vertex Q_v = cointriQ.signpostTri->intrinsicMesh->vertex(i);
        Vertex P_v = cointriQ.correspondingVertex[Q_v];
        SurfacePoint p_sp = cointriP.signpostTri->vertexLocations[P_v];

        if (p_sp.type == SurfacePointType::Vertex) continue;

        if (p_sp.type == SurfacePointType::Edge)
          p_sp = convertEdgePointToFacePoint(p_sp);

        Vector3 p_xyz = p_sp.interpolate(mdataP.geometry->inputVertexPositions);

        // extract barycentric coordinate delta
        SurfacePoint delta_sp{p_sp.face, {}};
        for (size_t di = 0; di < 2; ++di) {
          delta_sp.faceCoords[di] = -(*transportedGradient[k])[2 * (offset + i) + di];
        }
        delta_sp.faceCoords[2] = -(delta_sp.faceCoords[0] + delta_sp.faceCoords[1]);
        Vector3 delta_xyz = delta_sp.interpolate(mdataP.geometry->inputVertexPositions);
        delta_xyz *= tweaks.f["vectorFieldLineScale"];
        glColor3d(0.5, 0.5, 0.5);
        glVertex3d(p_xyz);
        glColor3d(gradientColors[k]);
        glVertex3d(p_xyz + delta_xyz);
      }
    }
    glEnd();
  }

  if (tweaks.f["vectorFieldLineWidth"] && tweaks.b["smoothedGradient"]) {
    const size_t offset = &mdataP == &mdataA ? mdataA.nV : 0;
    glLineWidth(tweaks.f["vectorFieldLineWidth"]);
    glBegin(GL_LINES);
    for (size_t i = 0; i < mdataQ.nV; ++i) {
      Vertex Q_v = cointriQ.signpostTri->intrinsicMesh->vertex(i);
      Vertex P_v = cointriQ.correspondingVertex[Q_v];
      SurfacePoint p_sp = cointriP.signpostTri->vertexLocations[P_v];

      if (p_sp.type == SurfacePointType::Vertex) continue;

      if (p_sp.type == SurfacePointType::Edge)
        p_sp = convertEdgePointToFacePoint(p_sp);

      Vector3 p_xyz = p_sp.interpolate(mdataP.geometry->inputVertexPositions);

      // extract barycentric coordinate delta
      SurfacePoint delta_sp{p_sp.face, {}};
      for (size_t di = 0; di < 2; ++di) {
        delta_sp.faceCoords[di] = -currentConfig.smoothedGradient[2 * (offset + i) + di];
      }
      delta_sp.faceCoords[2] = -(delta_sp.faceCoords[0] + delta_sp.faceCoords[1]);
      Vector3 delta_xyz = delta_sp.interpolate(mdataP.geometry->inputVertexPositions);
      delta_xyz *= tweaks.f["vectorFieldLineScale"];
      glColor3d(0.5, 0.5, 0.5);
      glVertex3d(p_xyz);
      glColor3d(1, 0.8, 0.5);
      glVertex3d(p_xyz + delta_xyz);
    }
    glEnd();
  }
}

} // namespace cit

void cit::draw() {
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
    drawCoInTri(currentConfig.cointriA, currentConfig.cointriB);
  });

  setCameraMatrix(mdataB);
  drawViewportBorder();
  mdataB.dispList.render([]() {
    drawCoInTri(currentConfig.cointriB, currentConfig.cointriA);
  });

  ImGui::Begin("Draw tweaks                  ");
  if (ImGui::ColorEdit3("Input face color", &tweaks.vec3["inputFaceColor"][0])) invalidateDispList();
  if (ImGui::Checkbox("Show axis", &tweaks.b["showAxis"])) invalidateDispList();
  if (ImGui::Checkbox("Show faces", &tweaks.b["showFaces"])) invalidateDispList();
  if (ImGui::ColorEdit3("A's input edge color", &mdataA.inputEdgeColor[0])) invalidateDispList();
  if (ImGui::ColorEdit3("B's input edge color", &mdataB.inputEdgeColor[0])) invalidateDispList();
  if (ImGui::Checkbox("High quality", &tweaks.b["highQuality"])) invalidateDispList();
  if (ImGui::SliderInt("Random seed", &tweaks.i["seed"], 0, 100)) invalidateDispList();

  if (ImGui::TreeNode("Vertex")) {
    if (ImGui::SliderFloat("Radius", &tweaks.f["vertexRadius"], 0.f, 0.02f, "%.5f")) invalidateDispList();
    if (ImGui::Checkbox("Use model color", &tweaks.b["vertexUseModelColor"])) invalidateDispList();
    if (ImGui::Checkbox("Hide original", &tweaks.b["vertexHideOriginal"])) invalidateDispList();
    if (ImGui::Checkbox("Hide inserted", &tweaks.b["vertexHideInserted"])) invalidateDispList();
    if (ImGui::SliderFloat("Silhouette ratio", &tweaks.f["vertexSilhouetteRatio"], 0.f, 2.f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Silhouette offset", &tweaks.f["vertexSilhouetteOffset"], 0.f, 2.f, "%.5f")) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Input edge")) {
    if (ImGui::SliderFloat("Radius", &tweaks.f["inputEdgeRadius"], 0.f, 0.01f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Offset", &tweaks.f["inputEdgeOffset"], 0.f, 5.f, "%.5f")) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Intrinsic edge")) {
    if (ImGui::SliderFloat("Radius (compatible)", &tweaks.f["compatibleEdgeRadius"], 0.f, 0.01f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Radius (incompatible)", &tweaks.f["incompatibleEdgeRadius"], 0.f, 0.01f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Offset", &tweaks.f["intrinsicEdgeOffset"], 0.f, 5.f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Incompatible edge brightness", &tweaks.f["incompatibleEdgeBrightness"], 0.f, 0.5f)) invalidateDispList();
    if (ImGui::Checkbox("Incomatible patch", &tweaks.b["showIncompatiblePatchBoundary"])) invalidateDispList();
    if (ImGui::Checkbox("Draw compatible edges in black", &tweaks.b["drawEdgeInBlack"])) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Mapped edge")) {
    if (ImGui::SliderFloat("Radius", &tweaks.f["mappedEdgeRadius"], 0.f, 0.01f, "%.5f")) invalidateDispList();
    if (ImGui::SliderFloat("Offset", &tweaks.f["mappedEdgeOffset"], 0.f, 5.f, "%.5f")) invalidateDispList();
    ImGui::TreePop();
  }

  if (ImGui::DragInt("Highlighted vertex ID", &cit::highlightedVertexID, 1.0f, -INT_MAX, INT_MAX)) invalidateDispList();
  if (ImGui::DragInt("Highlighted face ID", &cit::highlightedFaceID, 1.0f, -1, INT_MAX)) invalidateDispList();
  if (ImGui::SliderFloat("Highlighted edge width", &tweaks.f["highlightedEdgeWidth"], 0.f, 100.f)) invalidateDispList();
  if (ImGui::SliderFloat("Vector field line width", &tweaks.f["vectorFieldLineWidth"], 0.f, 100.f)) invalidateDispList();
  if (ImGui::SliderFloat("Vector field line scale", &tweaks.f["vectorFieldLineScale"], 0.f, 1000.f, "%.30f", 100.f)) invalidateDispList();
  if (ImGui::SliderFloat("Log max energy density", &tweaks.f["logMaxEnergyDensity"], std::log10(4.1), 4.f, "%.3f", 5.f)) invalidateDispList();

  ImGui::Spacing();

  if (ImGui::TreeNode("Texture")) {
    if (ImGui::SliderFloat("Scale", &tweaks.f["textureScale"], 0.f, 10.f)) invalidateDispList();
    if (ImGui::SliderFloat("Offset U", &tweaks.vec3["textureOffset"][0], 0.f, 1.f)) invalidateDispList();
    if (ImGui::SliderFloat("Offset V", &tweaks.vec3["textureOffset"][1], 0.f, 1.f)) invalidateDispList();
    ImGui::TreePop();
  }

  ImGui::Spacing();
  if (ImGui::BeginCombo("Fill mode A", tweaks.s["fillModeA"].c_str())) {
    for (size_t i = 0; i < fillModes.size(); i++) {
      if (ImGui::Selectable(fillModes[i].c_str(), tweaks.s["fillModeA"] == fillModes[i])) {
        tweaks.s["fillModeA"] = fillModes[i];
        mdataA.dispList.invalidate();
      }
    }
    ImGui::EndCombo();
  }
  if (ImGui::BeginCombo("Fill mode B", tweaks.s["fillModeB"].c_str())) {
    for (size_t i = 0; i < fillModes.size(); i++) {
      if (ImGui::Selectable(fillModes[i].c_str(), tweaks.s["fillModeB"] == fillModes[i])) {
        tweaks.s["fillModeB"] = fillModes[i];
        mdataB.dispList.invalidate();
      }
    }
    ImGui::EndCombo();
  }

  ImGui::Separator();
  if (ImGui::Checkbox("Descent direction", &tweaks.b["descentDirection"])) { mdataA.dispList.invalidate(); mdataB.dispList.invalidate(); }
  ImGui::Separator();
  if (ImGui::Checkbox("Previous gradient", &tweaks.b["previousGradient"])) { mdataA.dispList.invalidate(); mdataB.dispList.invalidate(); }
  if (ImGui::Checkbox("Show all", &tweaks.b["previousGradientAll"])) { mdataA.dispList.invalidate(); mdataB.dispList.invalidate(); }
  if (ImGui::SliderInt("Selected index in history", &tweaks.i["previousGradientSelected"], 0, sysParam.historyMaxSize)) invalidateDispList();
  ImGui::Separator();
  if (ImGui::Checkbox("Transported gradient", &tweaks.b["transportedGradient"])) { mdataA.dispList.invalidate(); mdataB.dispList.invalidate(); }
  if (ImGui::Checkbox("Smoothed gradient", &tweaks.b["smoothedGradient"])) { mdataA.dispList.invalidate(); mdataB.dispList.invalidate(); }

  // ImGui::Separator();
  // ImGui::SliderScalar("Angle threshold", ImGuiDataType_Double, &sysParam.angleThreshold, staticValuePtr(0.), staticValuePtr(3.), "%.7f", 2.0f);

  static char cameraParamBuf[512];
  {
    std::ostringstream oss;
    oss << " " << mdataA.camera.eye.transpose();
    oss << " " << mdataA.camera.center.transpose();
    oss << " " << mdataA.camera.up.transpose();
    oss << " " << mdataB.camera.eye.transpose();
    oss << " " << mdataB.camera.center.transpose();
    oss << " " << mdataB.camera.up.transpose();
    CIT_ASSERT(oss.str().size() < 512);
    strncpy(cameraParamBuf, oss.str().c_str(), oss.str().size());
    cameraParamBuf[oss.str().size()] = '\0';
  }
  if (ImGui::InputText("Camera parameters", cameraParamBuf, 512)) {
    std::vector<std::string> s;
    std::istringstream iss(cameraParamBuf);
    while (!iss.eof()) {
      s.emplace_back();
      iss >> s.back();
    }
    if (s.size() == 18) {
      std::array<double, 18> d;
      bool success = true;
      for (size_t i = 0; i < 18; ++i) {
        try {
          d[i] = std::stof(s[i]);
        } catch (...) {
          success = false;
          break;
        }
      }
      if (success) {
        size_t i = 0;
        mdataA.camera.eye    << d[i], d[i + 1], d[i + 2];   i += 3;
        mdataA.camera.center << d[i], d[i + 1], d[i + 2];   i += 3;
        mdataA.camera.up     << d[i], d[i + 1], d[i + 2];   i += 3;
        mdataB.camera.eye    << d[i], d[i + 1], d[i + 2];   i += 3;
        mdataB.camera.center << d[i], d[i + 1], d[i + 2];   i += 3;
        mdataB.camera.up     << d[i], d[i + 1], d[i + 2];   i += 3;
      }
    }
  }

  if (ImGui::TreeNode("Light parameters")) {
    static int lightIndex = 0;
    ImGui::SliderInt("Index", &lightIndex, 0, 7);
    if (ImGui::Checkbox("Enabled", &lightParam[lightIndex].enabled)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::SliderFloat("Ambient", &lightParam[lightIndex].ambient, 0.f, 1.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::SliderFloat("Diffuse", &lightParam[lightIndex].diffuse, 0.f, 1.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::SliderFloat("Specular", &lightParam[lightIndex].specular, 0.f, 1.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::gizmo3D("Direction", *(vec3*)&lightParam[lightIndex].position[0])) { invalidateDispList(); dumpLightMaterialParams(); }
    ImGui::TreePop();
  }

  if (ImGui::TreeNode("Material parameters")) {
    if (ImGui::SliderFloat("Ambient", &materialParam.ambient, 0.f, 1.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::SliderFloat("Specular", &materialParam.specular, 0.f, 1.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    if (ImGui::SliderFloat("Shininess", &materialParam.shininess, 0.f, 128.f)) { invalidateDispList(); dumpLightMaterialParams(); }
    ImGui::TreePop();
  }

  ImGui::End();
}
