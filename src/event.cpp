#include "cit.hpp"
#include "geometrycentral/surface/trace_geodesic.h"

#include "kt84/glfw_util.hh"
#include "kt84/graphics/graphics_util.hh"
using namespace kt84::graphics_util;

void mapCursorPos(GLFWwindow* window, int &x, int &y) {
    int wWidth, wHeight;
    glfwGetWindowSize(window, &wWidth, &wHeight);
    int fbWidth, fbHeight;
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    x = x * fbWidth  / wWidth;
    y = y * fbHeight / wHeight;
}

void cit::callback_framebuffersize(GLFWwindow* window, int width, int height) {
  window_width = width;
  window_height = height;
  mdataA.camera.reshape(width, height);
  mdataB.camera.reshape(width, height);
}

void cit::callback_mousebutton(GLFWwindow* window, int button, int action, int mods) {
  CoInTri& cointriA = currentConfig.cointriA;
  CoInTri& cointriB = currentConfig.cointriB;

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
      double sz = (double)read_depth(mouse.x, window_height - mouse.y);
      if (sz < 1.) {
        ModelData& mdata = clickedOnLeftSide ? mdataA : mdataB;
        CoInTri& cointri = clickedOnLeftSide ? cointriA : cointriB;
        setCameraMatrix(mdata);
        Vector3 p = Vector3::castFrom(unproject(Vector3d{(double)mouse.x, (double)(window_height - mouse.y), sz}));
        if (modFlag.shift) {
#if 0
          if (cointriA.uniqueID_to_Vertex.count(highlightedVertexID)) {
            // relocate highlighted vertex to clicked position
            Vertex A_v = cointriA.uniqueID_to_Vertex.at(highlightedVertexID);
            Vertex B_v = cointriB.uniqueID_to_Vertex.at(highlightedVertexID);
            if (clickedOnLeftSide) {
              if (A_v.getIndex() < mdataA.nV) {
                // A_v is fixed so cannot be relocated
                return;
              }
              // identify clicked position
              kt84::MinSelector<SurfacePoint> clickedPosition;
              for (Face A_f_input : cointriA.signpostTri->inputMesh.faces()) {
                std::array<Vector3, 3> q;
                int i = 0;
                for (Vertex A_v_input : A_f_input.adjacentVertices()) {
                  q[i++] = mdataA.geometry->inputVertexPositions[A_v_input];
                }

                // projection by least squares:
                //    q0*l0 + q1*l1 + q2*(1-l0-l1) ~= p
                //
                //    |q0-q2  q1-q2| * |l0| ~= |p-q2|
                //    |            |   |l1|    |    |
                Eigen::Matrix<double, 3, 2> M;
                M.col(0) = (q[0] - q[2]).castTo<Vector3d>();
                M.col(1) = (q[1] - q[2]).castTo<Vector3d>();
                Eigen::Matrix<double, 2, 3> Mt = M.transpose();
                Matrix2d MtM = Mt * M;
                Vector2d l = MtM.inverse() * Mt * (p - q[2]).castTo<Vector3d>();
                if (l[0] < 0 || l[1] < 0 || 1 - l[0] - l[1] < 0) {
                  // skip if outside the triangle
                  continue;
                }
                // pick the closest one
                SurfacePoint sp{A_f_input, Vector3{l[0], l[1], 1 - l[0] - l[1]}};
                Vector3 sp_xyz = sp.interpolate(mdataA.geometry->inputVertexPositions);
                clickedPosition.update(norm(p - sp_xyz), sp);
              }
              SurfacePoint pointOnIntrinsic = cointriA.signpostTri->equivalentPointOnIntrinsic(clickedPosition.value);
              auto minElement = std::min_element(&pointOnIntrinsic.faceCoords[0], &pointOnIntrinsic.faceCoords[0] + 3);
              if (*minElement < 0.1) {
                *minElement = 0.;
                pointOnIntrinsic.faceCoords /= sum(pointOnIntrinsic.faceCoords);
                pointOnIntrinsic = pointOnIntrinsic.reduced();
                CIT_ASSERT(pointOnIntrinsic.type == SurfacePointType::Edge);
                DLOG_INFO(0, "Relocating to an edge point {}", pointOnIntrinsic);
              } else {
                CIT_ASSERT(pointOnIntrinsic.type == SurfacePointType::Face);
                DLOG_INFO(0, "Relocating to a face point {}", pointOnIntrinsic);
              }
              if (cointriA.signpostTri->relocateInsertedVertex(A_v, pointOnIntrinsic)) {
                computeEdgePath(currentConfig);
              } else {
                DLOG_INFO(0, "Failed to relocate vertex {}", A_v);
              }
            } else {
              // TODO
            }
          }
#elif 1
          // select intrinsic edge
          kt84::MinSelector<Edge> nearestEdge;
          for (Edge e : cointri.signpostTri->intrinsicMesh->edges()) {
            Vertex v0 = e.halfedge().vertex();
            Vertex v1 = e.halfedge().twin().vertex();
            Vector3 q0 = cointri.signpostTri->vertexLocations[v0].interpolate(mdata.geometry->inputVertexPositions);
            Vector3 q1 = cointri.signpostTri->vertexLocations[v1].interpolate(mdata.geometry->inputVertexPositions);
            nearestEdge.update(norm(p - (q0 + q1) / 2.), e);
          }
          highlightedEdgeID = getUniqueEdgeID(cointri, nearestEdge.value);
#else
          // select input edge
          kt84::MinSelector<Edge> nearestEdge;
          for (Edge e : mdata.mesh->edges()) {
            Vertex v0 = e.halfedge().vertex();
            Vertex v1 = e.halfedge().twin().vertex();
            Vector3 q0 = mdata.geometry->inputVertexPositions[v0];
            Vector3 q1 = mdata.geometry->inputVertexPositions[v1];
            nearestEdge.update(norm(p - (q0 + q1) / 2.), e);
          }
          highlightedEdgeID = {nearestEdge.value.halfedge().vertex().getIndex(), nearestEdge.value.halfedge().twin().vertex().getIndex()};

#endif
        } else {
          kt84::MinSelector<Vertex> nearestVertex;
          for (Vertex v : cointri.signpostTri->intrinsicMesh->vertices()) {
            Vector3 q = cointri.signpostTri->vertexLocations[v].interpolate(mdata.geometry->inputVertexPositions);
            nearestVertex.update(norm(p - q), v);
          }
          highlightedVertexID = cointri.uniqueID_per_Vertex[nearestVertex.value];
        }
        mdataA.dispList.invalidate();
        mdataB.dispList.invalidate();
      }
    }
  } else if (actionFlag.release) {
    if (camera_active) {
      // Print camera parameters
      std::cout << "Camera param:" << std::endl
        << mdataA.camera.eye   .transpose() << " "
        << mdataA.camera.center.transpose() << " "
        << mdataA.camera.up    .transpose() << " "
        << mdataB.camera.eye   .transpose() << " "
        << mdataB.camera.center.transpose() << " "
        << mdataB.camera.up    .transpose() << std::endl;

      if (camera_active->drag_mode == Camera::DragMode::PAN) {
        double sx = (double)window_width * (camera_active == &mdataA.camera ? 0.25 : 0.75);
        double sy = (double)window_height * 0.5;
        double sz = (double)read_depth(sx, sy);
        if (sz < 1.) {
          setCameraMatrix(camera_active == &mdataA.camera ? mdataA : mdataB);
          camera_active->update_center(unproject(Vector3d{sx, sy, sz}));
        }
      }
      mdataA.dispList.invalidate();
      mdataB.dispList.invalidate();
      camera_active->mouse_up();
      camera_active = nullptr;
    }
  }
}

void cit::callback_cursorpos(GLFWwindow* window, double xpos, double ypos) {
  int x = static_cast<int>(xpos);
  int y = static_cast<int>(ypos);
  mapCursorPos(window, x, y);

  if (camera_active) {
    camera_active->mouse_move(x, y);
  }
}

void cit::callback_key(GLFWwindow* window, int key, int scancode, int action, int mods) {
  CoInTri& cointriA = currentConfig.cointriA;
  CoInTri& cointriB = currentConfig.cointriB;

  auto actionFlag = kt84::glfw_util::parseAction(action);
  if (!actionFlag.press) return;
  if (key == GLFW_KEY_K) {
    Vector3d A_center = 0.5 * (mdataA.bbox.min() + mdataA.bbox.max());
    Vector3d B_center = 0.5 * (mdataB.bbox.min() + mdataB.bbox.max());
    mdataA.camera.init(A_center + Vector3d(0,0,mdataA.bbox.diagonal().norm() * 1.2), A_center, Vector3d::UnitY());
    mdataB.camera.init(B_center + Vector3d(0,0,mdataB.bbox.diagonal().norm() * 1.2), B_center, Vector3d::UnitY());
  }
  if (key == GLFW_KEY_B) {
    if (!cointriA.uniqueID_to_Vertex.count(highlightedVertexID))
      return;
    // let P the model whose distance between camera and highlited vertex is smaller than the other
    // then set Q's camera such that it looks at the corresponding vertex
    auto getDistanceFromEyeToHighlightedVertex = [](const ModelData& mdata, const CoInTri& cointri) -> std::tuple<double, SurfacePoint> {
      SurfacePoint sp = cointri.signpostTri->vertexLocations[cointri.uniqueID_to_Vertex.at(highlightedVertexID)];
      Vector3 p = sp.interpolate(mdata.geometry->inputVertexPositions);
      double d = norm(p - Vector3::castFrom(mdata.camera.eye));
      return std::tie(d, sp);
    };
    auto lookAtHighlightedVertex = [](ModelData& mdata, SurfacePoint sp, double d) {
      Vector3 p = sp.interpolate(mdata.geometry->inputVertexPositions);
      Vector3 n = Vector3::zero();
      if (sp.type == SurfacePointType::Vertex) {
        n = mdata.geometry->vertexNormals[sp.vertex];
      } else if (sp.type == SurfacePointType::Edge) {
        n += mdata.geometry->faceNormals[sp.edge.halfedge().face()];
        n += mdata.geometry->faceNormals[sp.edge.halfedge().twin().face()];
        n = n.normalize();
      } else {
        n = mdata.geometry->faceNormals[sp.face];
      }
      mdata.camera.init((p + d * n).castTo<Vector3d>(), p.castTo<Vector3d>(), Vector3d(0,1,0));
      mdata.dispList.invalidate();
    };
    double A_dist, B_dist;
    SurfacePoint A_sp, B_sp;
    std::tie(A_dist, A_sp) = getDistanceFromEyeToHighlightedVertex(mdataA, cointriA);
    std::tie(B_dist, B_sp) = getDistanceFromEyeToHighlightedVertex(mdataB, cointriB);
    if (A_dist < B_dist) {
      lookAtHighlightedVertex(mdataB, B_sp, A_dist);
    } else {
      lookAtHighlightedVertex(mdataA, A_sp, B_dist);
    }
  }
  if (key == GLFW_KEY_F) {
    Edge A_e = getEdgeByUniqueID(cointriA, highlightedEdgeID);
    Edge B_e = getEdgeByUniqueID(cointriB, highlightedEdgeID);
    if (A_e != Edge() && B_e != Edge()) {
      if (isEdgeFlippable(cointriA, A_e) && isEdgeFlippable(cointriB, B_e)) {
        double energyBefore = 0;
        if (currentConfig.topologyValid)
          energyBefore = currentConfig.energy;

        cointriA.signpostTri->flipEdgeIfPossible(A_e);
        cointriB.signpostTri->flipEdgeIfPossible(B_e);

        highlightedEdgeID = getUniqueEdgeID(cointriA, A_e);

        updateCorrespondence(currentConfig);

        cointriA.signpostTri->refreshQuantities();
        cointriB.signpostTri->refreshQuantities();

        if (currentConfig.topologyValid) {
          double energyAfter = computeEnergy(currentConfig);
          double energyDelta = energyAfter - energyBefore;
          DLOG_INFO(0, "energy delta: {}", energyDelta);
        }

        computeEdgePath(currentConfig);
        computeOverlayPolygons(currentConfig);
        mdataA.dispList.invalidate();
        mdataB.dispList.invalidate();
      }
    }

  }
  if (key == GLFW_KEY_SPACE) {
    tweaks.b["showFaces"] = !tweaks.b["showFaces"];
    mdataA.dispList.invalidate();
    mdataB.dispList.invalidate();
  }
  if (key == GLFW_KEY_R) {
    std::mt19937 gen(time(nullptr));
    for (Edge A_e : cointriA.signpostTri->intrinsicMesh->edges()) {
      if (isEdgeCompatible(cointriA, A_e)) {
        Edge B_e = cointriA.correspondingEdge[A_e];

        if (!isEdgeFlippable(cointriA, A_e)) continue;
        if (!isEdgeFlippable(cointriB, B_e)) continue;

        if (std::uniform_real_distribution<>(0.0, 1.0)(gen) < 0.5) continue;

        cointriA.signpostTri->flipEdgeIfPossible(A_e);
        cointriB.signpostTri->flipEdgeIfPossible(B_e);
      }
    }

    cointriA.signpostTri->refreshQuantities();
    cointriB.signpostTri->refreshQuantities();

    updateCorrespondence(currentConfig);

    computeEdgePath(currentConfig);
    computeOverlayPolygons(currentConfig);

    computeEnergy(currentConfig);

    mdataA.dispList.invalidate();
    mdataB.dispList.invalidate();
  }
  if (key == GLFW_KEY_E) {
    flipToMinimizeEnergy(currentConfig);
    computeEdgePath(currentConfig);
    computeOverlayPolygons(currentConfig);
    mdataA.dispList.invalidate();
    mdataB.dispList.invalidate();
  }
  if (key == GLFW_KEY_D) {
    flipToMaximizeMinAngle(currentConfig);
    computeEdgePath(currentConfig);
    computeOverlayPolygons(currentConfig);
    mdataA.dispList.invalidate();
    mdataB.dispList.invalidate();
  }
}
