#include "cit.hpp"
#include "geometrycentral/surface/trace_geodesic.h"

void cit::insertVertices(const VectorXsp& P_x, const CoInTri& cointriP, CoInTri& cointriQ, spdlog::logger_ptr logger, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN insertVertices");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   insertVertices"); });

  LOG_INFO(logger, logDepth, "Inserting vertices (from {} to {})", cointriP.mdata->name, cointriQ.mdata->name);
  const size_t P_nV = cointriP.mdata->nV;
  const size_t Q_nV = cointriQ.mdata->nV;
  CIT_ASSERT(P_nV != Q_nV);
  CIT_ASSERT(P_x.size() == (int)P_nV);

  cointriQ.uniqueID_per_Vertex = VertexData<int>(*cointriQ.signpostTri->intrinsicMesh, 0);

  // store unique IDs of Q's vertices
  for (Vertex Q_v : cointriQ.signpostTri->intrinsicMesh->vertices()) {
    int uniqueID = (P_nV < Q_nV ? -1 : 1) * (Q_v.getIndex() + 1);
    cointriQ.uniqueID_per_Vertex[Q_v] = uniqueID;
  }

  // Insert P's vertices to Q
  kt84::MaxMinAverage mma;
  for (size_t i = 0; i < P_nV; ++i) {
    SurfacePoint pointOnInput = P_x[i];
    SurfacePoint pointOnIntrinsic = cointriQ.signpostTri->equivalentPointOnIntrinsic(pointOnInput);
    if (pointOnIntrinsic.type == SurfacePointType::Vertex) {
      int uniqueID = P_nV < Q_nV ?
        i + 1 + (P_nV + 1) * (P_x[i].vertex.getIndex() + 1) :
        P_x[i].vertex.getIndex() + 1 + (Q_nV + 1) * (i + 1);
      cointriQ.uniqueID_per_Vertex[pointOnIntrinsic.vertex] = uniqueID;
    } else {
      Vertex newV = cointriQ.signpostTri->insertVertex(pointOnIntrinsic);
      int uniqueID = (P_nV < Q_nV ? 1 : -1) * (i + 1);
      cointriQ.uniqueID_per_Vertex[newV] = uniqueID;
      // sanity check
      SurfacePoint P_x_check = cointriQ.signpostTri->equivalentPointOnInput(newV);
      Vector3 p0 = P_x[i].interpolate(cointriQ.mdata->geometry->inputVertexPositions);
      Vector3 p1 = P_x_check.interpolate(cointriQ.mdata->geometry->inputVertexPositions);
      mma.update(norm(p0 - p1));
    }
  }
  if (mma.max() < 1.e-5)
    LOG_DEBUG(logger, logDepth, "  Equivalent point sanity check: {}", mma.max());
  else
    LOG_WARN(logger, logDepth, "  Equivalent point sanity check: {}", mma.max());
}

void cit::updateCorrespondence(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN updateCorrespondence");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   updateCorrespondence"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  // Update uniqueID_to_Vertex
  auto updateCorrespondence_inner1 = [](CoInTri& cointri) {
    cointri.uniqueID_to_Vertex.clear();
    for (Vertex v : cointri.signpostTri->intrinsicMesh->vertices()) {
      int u = cointri.uniqueID_per_Vertex[v];
      CIT_ASSERT(!cointri.uniqueID_to_Vertex.count(u));
      cointri.uniqueID_to_Vertex[u] = v;
    }
  };
  updateCorrespondence_inner1(cointriA);
  updateCorrespondence_inner1(cointriB);
  CIT_ASSERT(cointriA.uniqueID_to_Vertex.size() == cointriB.uniqueID_to_Vertex.size());

  // Make vertex correspondence
  cointriA.correspondingVertex = VertexData<Vertex>(*cointriA.signpostTri->intrinsicMesh);
  cointriB.correspondingVertex = VertexData<Vertex>(*cointriB.signpostTri->intrinsicMesh);
  for (Vertex A_v : cointriA.signpostTri->intrinsicMesh->vertices()) {
    Vertex B_v = cointriB.uniqueID_to_Vertex.at(cointriA.uniqueID_per_Vertex[A_v]);
    cointriA.correspondingVertex[A_v] = B_v;
    cointriB.correspondingVertex[B_v] = A_v;
  }

  // Make edge correspondence
  cointriA.correspondingEdge = EdgeData<Edge>(*cointriA.signpostTri->intrinsicMesh);
  cointriB.correspondingEdge = EdgeData<Edge>(*cointriB.signpostTri->intrinsicMesh);
  size_t nCompatibleEdges = 0;
  for (Edge A_e : cointriA.signpostTri->intrinsicMesh->edges()) {
    std::pair<int, int> ue = getUniqueEdgeID(cointriA, A_e);
    Vertex B_v0 = cointriB.uniqueID_to_Vertex.at(ue.first);
    Vertex B_v1 = cointriB.uniqueID_to_Vertex.at(ue.second);
    Edge B_e = B_v0.connectingEdge(B_v1);
    if (B_e != Edge()) {
      cointriA.correspondingEdge[A_e] = B_e;
      cointriB.correspondingEdge[B_e] = A_e;
      ++nCompatibleEdges;
    }
  }
  if (config.nCompatibleEdges != nCompatibleEdges) {
    config.nCompatibleEdges = nCompatibleEdges;
    LOG_INFO(logger, logDepth, "  Compatible edges: {} / {}", nCompatibleEdges, config.cointriA.signpostTri->intrinsicMesh->nEdges());
  }

  // Make halfedge correspondence
  cointriA.correspondingHalfedge = HalfedgeData<Halfedge>(*cointriA.signpostTri->intrinsicMesh);
  cointriB.correspondingHalfedge = HalfedgeData<Halfedge>(*cointriB.signpostTri->intrinsicMesh);
  for (Edge A_e : cointriA.signpostTri->intrinsicMesh->edges()) {
    if (!isEdgeCompatible(cointriA, A_e))
      continue;

    Edge B_e = cointriA.correspondingEdge[A_e];
    Halfedge A_he = A_e.halfedge();
    Halfedge B_he = B_e.halfedge();
    if (cointriA.correspondingVertex[A_he.vertex()] != B_he.vertex()) {
      B_he = B_he.twin();
      CIT_ASSERT(cointriA.correspondingVertex[A_he.vertex()] == B_he.vertex());
    }

    cointriA.correspondingHalfedge[A_he] = B_he;
    cointriB.correspondingHalfedge[B_he] = A_he;
    cointriA.correspondingHalfedge[A_he.twin()] = B_he.twin();
    cointriB.correspondingHalfedge[B_he.twin()] = A_he.twin();
  }

  // Make face correspondence
  cointriA.correspondingFace = FaceData<Face>(*cointriA.signpostTri->intrinsicMesh);
  cointriB.correspondingFace = FaceData<Face>(*cointriB.signpostTri->intrinsicMesh);
  cointriA.reversedFace = FaceData<Face>(*cointriA.signpostTri->intrinsicMesh);
  cointriB.reversedFace = FaceData<Face>(*cointriB.signpostTri->intrinsicMesh);
  size_t nCompatibleFaces = 0;
  size_t nReversedFaces = 0;
  for (Face A_f : cointriA.signpostTri->intrinsicMesh->faces()) {
    Halfedge A_he0 = A_f.halfedge();
    Halfedge A_he1 = A_he0.next();
    Halfedge A_he2 = A_he1.next();

    if (!isEdgeCompatible(cointriA, A_he0.edge())) continue;
    if (!isEdgeCompatible(cointriA, A_he1.edge())) continue;
    if (!isEdgeCompatible(cointriA, A_he2.edge())) continue;

    Halfedge B_he0 = cointriA.correspondingHalfedge[A_he0];
    Halfedge B_he1 = cointriA.correspondingHalfedge[A_he1];
    Halfedge B_he2 = cointriA.correspondingHalfedge[A_he2];

    bool isCompatible = true;

    if (B_he0.next() != B_he1) isCompatible = false;
    if (B_he1.next() != B_he2) isCompatible = false;
    if (B_he2.next() != B_he0) isCompatible = false;

    if (isCompatible) {
      Face B_f = B_he0.face();
      cointriA.correspondingFace[A_f] = B_f;
      cointriB.correspondingFace[B_f] = A_f;
      ++nCompatibleFaces;
      continue;
    }

    B_he0 = B_he0.twin();
    B_he1 = B_he1.twin();
    B_he2 = B_he2.twin();

    bool isReversed = true;
    /*
      A_f: P->Q->R
      B_f: R->Q->P

      A_he0 : P->Q
      A_he1 : Q->R
      A_he2 : R->P

      B_he0 : Q->P
      B_he1 : R->Q
      B_he2 : P->R
    */
    if (B_he0.next() != B_he2) isReversed = false;
    if (B_he1.next() != B_he0) isReversed = false;
    if (B_he2.next() != B_he1) isReversed = false;

    if (isReversed) {
      Face B_f = B_he0.face();
      cointriA.reversedFace[A_f] = B_f;
      cointriB.reversedFace[B_f] = A_f;
      ++nReversedFaces;
    }
  }
  if (config.nCompatibleFaces != nCompatibleFaces) {
    config.nCompatibleFaces = nCompatibleFaces;
    LOG_INFO(logger, logDepth, "  Compatible (reversed) faces: {} ({}) / {}", nCompatibleFaces, nReversedFaces, config.cointriA.signpostTri->intrinsicMesh->nFaces());
  }
}

size_t cit::mergeNearbyVertexPairs(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN mergeNearbyVertexPairs");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   mergeNearbyVertexPairs"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  size_t res = 0;

  while (true) {
    bool changed = false;

    for (Edge A_e : cointriA.signpostTri->intrinsicMesh->edges()) {
      if (!isEdgeCompatible(cointriA, A_e))
        continue;

      std::array<Vertex, 2> A_v = A_e.adjacentVertices();

      // Skip if either endpoint is already a merged vertex
      if (isVertexMerged(cointriA, A_v[0])) continue;
      if (isVertexMerged(cointriA, A_v[1])) continue;

      // Skip if both endpoints are inserted or original
      if ((A_v[0].getIndex() < mdataA.nV) + (A_v[1].getIndex() < mdataA.nV) != 1)
        continue;

      // Ensure A_v[0] is original vertex
      if (A_v[0].getIndex() >= mdataA.nV)
        std::swap(A_v[0], A_v[1]);

      std::array<Vertex, 2> B_v = {
        cointriA.correspondingVertex[A_v[0]],
        cointriA.correspondingVertex[A_v[1]]
      };

      // Get inserted vertices' locations on A & B
      SurfacePoint A_sp = cointriA.signpostTri->vertexLocations[A_v[1]];
      SurfacePoint B_sp = cointriB.signpostTri->vertexLocations[B_v[0]];
      CIT_ASSERT(A_sp.type != SurfacePointType::Vertex);
      CIT_ASSERT(B_sp.type != SurfacePointType::Vertex);

      // If edge point, treat it as face point
      A_sp = A_sp.inSomeFace();
      B_sp = B_sp.inSomeFace();

      Face A_inputFace = A_sp.face;
      Face B_inputFace = B_sp.face;

      std::array<Vertex, 3> A_inputFaceVertices = adjacentVertices(A_inputFace);
      std::array<Vertex, 3> B_inputFaceVertices = adjacentVertices(B_inputFace);

      // Check if the largest face coord component is close enough to 1
      size_t A_maxIndex;
      size_t B_maxIndex;
      max(A_sp.faceCoords, &A_maxIndex);
      max(B_sp.faceCoords, &B_maxIndex);

      if (A_sp.faceCoords[A_maxIndex] < 1. - sysParam.snapThreshold &&
          B_sp.faceCoords[B_maxIndex] < 1. - sysParam.snapThreshold) continue;

      // Check if the closest vertex is the same as the other vertex on A_e (B_e)
      if (A_inputFaceVertices[A_maxIndex].getIndex() != A_v[0].getIndex()) continue;
      if (B_inputFaceVertices[B_maxIndex].getIndex() != B_v[1].getIndex()) continue;

      //------------------+
      // Merge this edge! |
      //------------------+

      // Check if the two edges are both collapsible
      Halfedge A_he = A_v[0].connectingHalfedge(A_v[1]);
      Halfedge B_he = B_v[1].connectingHalfedge(B_v[0]);
      CIT_ASSERT(cointriA.correspondingHalfedge[A_he.twin()] == B_he);
      if (!cointriA.signpostTri->collapseInteriorEdge(A_he, true)) continue;
      if (!cointriB.signpostTri->collapseInteriorEdge(B_he, true)) continue;

      // Prevent degree-two vertex from being generated after merge
      if (A_v[0].degree() + A_v[1].degree() - 4 < 3) continue;
      if (B_v[0].degree() + B_v[1].degree() - 4 < 3) continue;

      if (A_he.next().tipVertex().degree() <= 3) continue;
      if (B_he.next().tipVertex().degree() <= 3) continue;
      if (A_he.twin().next().tipVertex().degree() <= 3) continue;
      if (B_he.twin().next().tipVertex().degree() <= 3) continue;

      // Identify unique vertex ID before/after merge
      std::array<int, 2> uniqueVertexID = {
        cointriA.uniqueID_per_Vertex[A_v[0]],
        cointriA.uniqueID_per_Vertex[A_v[1]]
      };
      int uniqueVertexID_new;
      if (mdataA.nV < mdataB.nV) {
        CIT_ASSERT(uniqueVertexID[0] == (int)(A_v[0].getIndex() + 1));
        CIT_ASSERT(uniqueVertexID[1] == -(int)(B_v[1].getIndex() + 1));
        uniqueVertexID_new = A_v[0].getIndex() + 1 + (mdataA.nV + 1) * (B_v[1].getIndex() + 1);
      } else {
        CIT_ASSERT(uniqueVertexID[0] == -(int)(A_v[0].getIndex() + 1));
        CIT_ASSERT(uniqueVertexID[1] == (int)(B_v[1].getIndex() + 1));
        uniqueVertexID_new = B_v[1].getIndex() + 1 + (mdataB.nV + 1) * (A_v[0].getIndex() + 1);
      }

      LOG_DEBUG(logger, logDepth, "Collapsing edge ({}, {})", uniqueVertexID[0], uniqueVertexID[1]);

      // Perform edge collpase
      cointriA.signpostTri->collapseInteriorEdge(A_he, false);
      cointriB.signpostTri->collapseInteriorEdge(B_he, false);

      // Assign new unique vertex ID (representing the overlap of AB vertices) to the collapsed vertex
      cointriA.uniqueID_per_Vertex[A_v[0]] = uniqueVertexID_new;
      cointriB.uniqueID_per_Vertex[B_v[1]] = uniqueVertexID_new;

      updateCorrespondence(config, logger, logDepth + 1);
      changed = true;
      ++res;
      break;
    }

    // Repeat until no more collapse is possible
    if (!changed) break;
  }

  if (res)
    LOG_INFO(logger, logDepth, "  Merged {} nearby vertex pairs", res);
  return res;
}

size_t cit::simpleFlip(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN simpleFlip");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   simpleFlip"); });

  auto simpleFlip_inner = [logger, logDepth](CoInTri& cointriP, const CoInTri& cointriQ) {
    size_t res = 0;
    for (Edge P_e : cointriP.signpostTri->intrinsicMesh->edges()) {
      if (isEdgeCompatible(cointriP, P_e) || !isEdgeFlippable(cointriP, P_e))
        continue;
      auto P_v = getOppositeVertices(P_e);
      auto ue = getUniqueEdgeID(cointriP, P_v);
      Edge Q_e = getEdgeByUniqueID(cointriQ, ue);
      if (Q_e != Edge()) {
        cointriP.signpostTri->flipEdgeIfPossible(P_e);
        ++res;
      }
    }
    if (res)
      LOG_INFO(logger, logDepth, "* Flipped {} edges (mutable {}, fixed {})", res, cointriP.mdata->name, cointriQ.mdata->name);
    return res;
  };

  size_t res_total = simpleFlip_inner(config.cointriA, config.cointriB);
  res_total += simpleFlip_inner(config.cointriB, config.cointriA);
  return res_total;
}

size_t cit::simpleCoFlip(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN simpleCoFlip");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   simpleCoFlip"); });

  size_t res = 0;
  for (Edge A_e : config.cointriA.signpostTri->intrinsicMesh->edges()) {
    if (isEdgeCompatible(config.cointriA, A_e) || !isEdgeFlippable(config.cointriA, A_e))
      continue;
    auto A_v = getOppositeVertices(A_e);
    auto A_ue = getUniqueEdgeID(config.cointriA, {A_v[0], A_v[1]});
    for (Edge B_e : config.cointriB.signpostTri->intrinsicMesh->edges()) {
      if (isEdgeCompatible(config.cointriB, B_e) || !isEdgeFlippable(config.cointriB, B_e))
        continue;
      auto B_v = getOppositeVertices(B_e);
      auto B_ue = getUniqueEdgeID(config.cointriB, {B_v[0], B_v[1]});
      if (A_ue == B_ue) {
        config.cointriA.signpostTri->flipEdgeIfPossible(A_e);
        config.cointriB.signpostTri->flipEdgeIfPossible(B_e);
        ++res;
        break;
      }
    }
  }
  if (res)
    LOG_INFO(logger, logDepth, "** Co-flipped {} edges", res);
  return res;
}

size_t cit::doubleFlipCompatible(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN doubleFlipCompatible");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   doubleFlipCompatible"); });

/*
    A                            |        B
          ____*____              |          ____*____              
       __/ _//|\   \__           |       __/   /|\\_ \__           
    __/  _/ | | |     \__        |    __/     | | | \_  \__        
   /    /  /  |  \       \       |   /       /  |  \  \    \       
  *-e2-*  e0  |  e1       *      |  *      e0   |  e1  *-e3-*      
   \   | /    *    \     /       |   \     /    *    \ |   /       
    \  | |  _/ \_  |    /        |    \    |  _/ \_  | |  /        
     \ |/ _/     \_ \  /         |     \  / _/     \_ \| /         
      \||/         \| /          |      \ |/         \||/          
       *-------------*           |       *-------------*           

*/
  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  size_t res = 0;
  for (Edge A_e0 : cointriA.signpostTri->intrinsicMesh->edges()) {
    if (!isEdgeCompatible(cointriA, A_e0) || !isEdgeFlippable(cointriA, A_e0))
      continue;

    // Check edges adjacent to the one ring of its adjacent vertices
    std::set<Edge> edgesToCheck;
    for (Vertex A_v : {
      A_e0.halfedge().vertex(),
      A_e0.halfedge().tipVertex(),
      A_e0.halfedge().next().tipVertex(),
      A_e0.halfedge().twin().next().tipVertex()
    }) {
      for (Halfedge A_he : A_v.outgoingHalfedges()) {
        edgesToCheck.insert(A_he.edge());
        edgesToCheck.insert(A_he.next().edge());
      }
    }

    for (Edge A_e1 : edgesToCheck) {
      if (!isEdgeCompatible(cointriA, A_e1) || !isEdgeFlippable(cointriA, A_e1))
        continue;

      Edge B_e0 = cointriA.correspondingEdge[A_e0];
      Edge B_e1 = cointriA.correspondingEdge[A_e1];
      if (!isEdgeFlippable(cointriB, B_e0)) continue;
      if (!isEdgeFlippable(cointriB, B_e1)) continue;

      if (getUniqueEdgeID(cointriA, getOppositeVertices(A_e0)) != getUniqueEdgeID(cointriB, getOppositeVertices(B_e1)))
        continue;

      auto ue2 = getUniqueEdgeID(cointriB, getOppositeVertices(B_e0));
      Vertex A_v20 = cointriA.uniqueID_to_Vertex.at(ue2.first);
      Vertex A_v21 = cointriA.uniqueID_to_Vertex.at(ue2.second);
      Edge A_e2 = A_v20.connectingEdge(A_v21);
      if (!A_e2.getMesh()) continue;

      auto ue3 = getUniqueEdgeID(cointriA, getOppositeVertices(A_e1));
      Vertex B_v30 = cointriB.uniqueID_to_Vertex.at(ue3.first);
      Vertex B_v31 = cointriB.uniqueID_to_Vertex.at(ue3.second);
      Edge B_e3 = B_v30.connectingEdge(B_v31);
      if (!B_e3.getMesh()) continue;

      // All checks passed, double-flip edges
      cointriA.signpostTri->flipEdgeIfPossible(A_e0);
      cointriA.signpostTri->flipEdgeIfPossible(A_e1);
      cointriB.signpostTri->flipEdgeIfPossible(B_e0);
      cointriB.signpostTri->flipEdgeIfPossible(B_e1);
      updateCorrespondence(config, logger, logDepth + 1);
      ++res;
      break;
    }
  }
  if (res)
    LOG_INFO(logger, logDepth, "*** Double-flipped {} compatible edge pairs", res);
  return res;
}

size_t cit::flipCompatible(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipCompatible");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipCompatible"); });

  size_t res = 0;
  for (Edge A_e : config.cointriA.signpostTri->intrinsicMesh->edges()) {
    if (!isEdgeCompatible(config.cointriA, A_e))
      continue;
    Edge B_e = config.cointriA.correspondingEdge[A_e];

    if (!isEdgeFlippable(config.cointriA, A_e)) continue;
    if (!isEdgeFlippable(config.cointriB, B_e)) continue;

    std::array<Vertex, 2> A_v = getOppositeVertices(A_e);
    std::array<Vertex, 2> B_v = getOppositeVertices(B_e);

    std::array<Vertex, 2> A_cv; // "cv" stands for "corresponding vertex"
    std::array<Vertex, 2> B_cv;
    for (int i = 0; i < 2; ++i) {
      B_cv[i] = config.cointriA.correspondingVertex[A_v[i]];
      A_cv[i] = config.cointriB.correspondingVertex[B_v[i]];
    }

    auto flipCompatibleEdgesIfPairedWithIncompatible_inner1 = [](const CoInTri& cointri, const std::array<Vertex, 2>& target_v) -> bool {
      // ok if there exists an incompatible edge between target_v[0] & target_v[1]
      Edge e1 = target_v[0].connectingEdge(target_v[1]);
      if (e1.getMesh() && !isEdgeCompatible(cointri, e1))
        return true;

      // also ok if there exists a flippable incompatible edge whose opposite vertices are target_v[0] & target_v[1]
      for (Halfedge he : target_v[0].outgoingHalfedges()) {
        if (he.next().twin().next().tipVertex() == target_v[1]) {
          Edge e2 = he.next().edge();
          return !isEdgeCompatible(cointri, e2) && isEdgeFlippable(cointri, e2);
        }
      }

      return false;
    };
    if (!flipCompatibleEdgesIfPairedWithIncompatible_inner1(config.cointriA, A_cv)) continue;
    if (!flipCompatibleEdgesIfPairedWithIncompatible_inner1(config.cointriB, B_cv)) continue;

    config.cointriA.signpostTri->flipEdgeIfPossible(A_e);
    config.cointriB.signpostTri->flipEdgeIfPossible(B_e);
    ++res;
  }
  if (res)
    LOG_INFO(logger, logDepth, "**** Flipped {} compatible edges paired with incompatible edges", res);
  return res;
}

namespace cit {
namespace detail {

struct FlatPolygon {
  std::set<Edge> interiorEdges;
  std::vector<int> boundaryVertexUniqueIDs;
  std::map<Vertex, Vector2> vertexPositions;    // In a temporary 2D coordinate system
};

std::vector<FlatPolygon> getFlatPolygons(const CoInTri& cointri) {
  std::vector<FlatPolygon> flatPolygons;

  EdgeData<char> visited{*cointri.signpostTri->intrinsicMesh, false};
  for (Edge e0 : cointri.signpostTri->intrinsicMesh->edges()) {
    if (isEdgeCompatible(cointri, e0) || visited[e0])
      continue;

    FlatPolygon flatPolygon;

    // Find flat polygons by front propagation
    std::deque<Edge> edgesToProcess;
    edgesToProcess.push_back(e0);
    while (!edgesToProcess.empty()) {
      Edge e1 = edgesToProcess.front();
      CIT_ASSERT(!visited[e1]);
      flatPolygon.interiorEdges.insert(e1);
      edgesToProcess.pop_front();
      visited[e1] = true;
      for (Halfedge he2 : {
        e1.halfedge().next(),
        e1.halfedge().next().next(),
        e1.halfedge().twin().next(),
        e1.halfedge().twin().next().next()
      }) {
        Edge e2 = he2.edge();
        if (std::find(edgesToProcess.begin(), edgesToProcess.end(), e2) != edgesToProcess.end())
          continue;
        if (!isFaceCompatible(cointri, he2.twin().face()) && !visited[e2])
          edgesToProcess.push_back(e2);
      }
    }

    // Ensure each vertex lies at the boundary of compatible/incompatible faces. otherwise, skip
    std::set<Vertex> vertexSet;
    for (Edge e : flatPolygon.interiorEdges) {
      vertexSet.insert(e.halfedge().vertex());
      vertexSet.insert(e.halfedge().next().vertex());
      vertexSet.insert(e.halfedge().next().next().vertex());
      vertexSet.insert(e.halfedge().twin().next().next().vertex());
    }
    bool found = false;
    for (Vertex v : vertexSet) {
      bool isInterior = true;
      for (Face f : v.adjacentFaces()) {
        if (isFaceCompatible(cointri, f)) {
          isInterior = false;
          break;
        }
      }
      if (isInterior) {
        found = true;
        break;
      }
    }
    if (found)
      continue;

    // Ensure every vertex has a half-disk topology in the flat polygon (i.e. boundary doesn't touch itself)
    found = false;
    for (Vertex v : vertexSet) {
      size_t count = 0;
      for (Halfedge he : v.outgoingHalfedges()) {
        count += isFaceCompatible(cointri, he.face()) ^ isFaceCompatible(cointri, he.twin().face());
      }
      CIT_ASSERT(count % 2 == 0);
      if (count > 2) {
        found = true;
        break;
      }
    }
    if (found)
      continue;

    // Collect counterclockwise ordered vertex list
    std::vector<Vertex> vertexList;
    Halfedge heStart;
    for (Edge e : flatPolygon.interiorEdges) {
      for (Halfedge he : {
        e.halfedge().next(),
        e.halfedge().next().next(),
        e.halfedge().twin().next(),
        e.halfedge().twin().next().next()
      }) {
        if (isFaceCompatible(cointri, he.twin().face())) {
          heStart = he;
          break;
        }
      }
      if (heStart.getMesh())
        break;
    }
    CIT_ASSERT(heStart.getMesh());
    Halfedge he = heStart;
    do {
      vertexList.push_back(he.vertex());
      he = he.next();
      while (!isFaceCompatible(cointri, he.twin().face())) {
        he = he.twin().next();
      }
    } while (he != heStart);

    // Skip if this face set doesn't have a disc topology
    if (vertexList.size() != vertexSet.size())
      continue;

    // Fill in boundaryVertexUniqueIDs
    for (Vertex v : vertexList) {
      flatPolygon.boundaryVertexUniqueIDs.push_back(cointri.uniqueID_per_Vertex[v]);
    }
    // Make it unique by rotating the smallest to front
    std::rotate(
      flatPolygon.boundaryVertexUniqueIDs.begin(),
      std::min_element(flatPolygon.boundaryVertexUniqueIDs.begin(), flatPolygon.boundaryVertexUniqueIDs.end()),
      flatPolygon.boundaryVertexUniqueIDs.end()
    );

    // Lay out vertices to a temporary 2D coordinate system
    Edge eFirst = *flatPolygon.interiorEdges.begin();
    flatPolygon.vertexPositions[eFirst.halfedge().vertex()] = {0., 0.};
    flatPolygon.vertexPositions[eFirst.halfedge().tipVertex()] = {cointri.signpostTri->edgeLengths[eFirst], 0.};
    while (true) {
      bool changed = false;

      for (Edge e : flatPolygon.interiorEdges) {
        std::array<Vertex, 2> v = e.adjacentVertices();
        std::array<Vertex, 2> vOpp = getOppositeVertices(e);

        if (!flatPolygon.vertexPositions.count(v[0])) continue;
        if (!flatPolygon.vertexPositions.count(v[1])) continue;

        Vector2 pA = flatPolygon.vertexPositions.at(v[0]);
        Vector2 pB = flatPolygon.vertexPositions.at(v[1]);

        if (!flatPolygon.vertexPositions.count(vOpp[0])) {
          double lBC = cointri.signpostTri->edgeLengths[e.halfedge().next().edge()];
          double lCA = cointri.signpostTri->edgeLengths[e.halfedge().next().next().edge()];
          flatPolygon.vertexPositions[vOpp[0]] = layoutTriangleVertexFromLength(pA, pB, lBC, lCA);
          changed = true;
        }

        if (!flatPolygon.vertexPositions.count(vOpp[1])) {
          double lAD = cointri.signpostTri->edgeLengths[e.halfedge().twin().next().edge()];
          double lDB = cointri.signpostTri->edgeLengths[e.halfedge().twin().next().next().edge()];
          flatPolygon.vertexPositions[vOpp[1]] = layoutTriangleVertexFromLength(pB, pA, lAD, lDB);
          changed = true;
        }
      }

      if (!changed)
        break;
    }
    CIT_ASSERT(flatPolygon.vertexPositions.size() == flatPolygon.boundaryVertexUniqueIDs.size());

    flatPolygons.push_back(flatPolygon);
  }
  return flatPolygons;
}

bool connectAllEdgesToVertex_check(const CoInTri& cointri, const FlatPolygon& flatPolygon, int uniqueVertexID) {
  Vertex v0 = cointri.uniqueID_to_Vertex.at(uniqueVertexID);
  Vector2 p0 = flatPolygon.vertexPositions.at(v0);

  const size_t n = flatPolygon.vertexPositions.size();
  for (size_t i = 0; i < n; ++i) {
    Vertex v1 = cointri.uniqueID_to_Vertex.at(flatPolygon.boundaryVertexUniqueIDs[i]);
    Vertex v2 = cointri.uniqueID_to_Vertex.at(flatPolygon.boundaryVertexUniqueIDs[(i + 1) % n]);

    if (v1 == v0 || v2 == v0)
      continue;

    Vector2 p1 = flatPolygon.vertexPositions.at(v1);
    Vector2 p2 = flatPolygon.vertexPositions.at(v2);

    if (cross(normalize(p1 - p0), normalize(p2 - p0)) <= 1.e-5)
      return false;
  }
  return true;
}

bool connectAllEdgesToVertex(CoInTri& cointri, const FlatPolygon& flatPolygon, int uniqueVertexID) {
  Vertex vTarget = cointri.uniqueID_to_Vertex.at(uniqueVertexID);
  while (true) {
    size_t nGood = 0;
    bool changed = false;
    for (Edge e : flatPolygon.interiorEdges) {
      if (e.halfedge().vertex() == vTarget || e.halfedge().tipVertex() == vTarget) {
        ++nGood;
        continue;
      }

      if (!isEdgeFlippable(cointri, e))
        continue;

      for (Vertex vOpp : getOppositeVertices(e)) {
        if (vOpp == vTarget) {
          cointri.signpostTri->flipEdgeIfPossible(e);
          ++nGood;
          changed = true;
          break;
        }
      }
    }
    if (nGood == flatPolygon.interiorEdges.size())
      return true;
    if (!changed)
      return false;
  }
}

} // namespace cit
} // namespace detail

size_t cit::flipFlatPolygon(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipFlatPolygon");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipFlatPolygon"); });

  size_t res = 0;

  std::vector<detail::FlatPolygon> A_flatPolygons = detail::getFlatPolygons(config.cointriA);
  std::vector<detail::FlatPolygon> B_flatPolygons = detail::getFlatPolygons(config.cointriB);

  // process each flat polygon
  for (size_t i = 0; i < A_flatPolygons.size(); ++i) {
    detail::FlatPolygon* A_flatPolygon = &A_flatPolygons[i];
    detail::FlatPolygon* B_flatPolygon = nullptr;
    for (size_t j = 0; j < B_flatPolygons.size(); ++j) {
      if (A_flatPolygon->boundaryVertexUniqueIDs == B_flatPolygons[j].boundaryVertexUniqueIDs) {
        B_flatPolygon = &B_flatPolygons[j];
        break;
      }
    }
    if (!B_flatPolygon)
      continue;

    for (int uniqueVertexID : A_flatPolygon->boundaryVertexUniqueIDs) {
      if (detail::connectAllEdgesToVertex_check(config.cointriA, *A_flatPolygon, uniqueVertexID) &&
          detail::connectAllEdgesToVertex_check(config.cointriB, *B_flatPolygon, uniqueVertexID))
      {
        bool A_success = detail::connectAllEdgesToVertex(config.cointriA, *A_flatPolygon, uniqueVertexID);
        bool B_success = detail::connectAllEdgesToVertex(config.cointriB, *B_flatPolygon, uniqueVertexID);
        if (!A_success) continue;        // CIT_ASSERT(A_success);    // Not sure why this doesn't hold sometimes
        if (!B_success) continue;        // CIT_ASSERT(B_success);
        res += A_flatPolygon->interiorEdges.size();
        break;
      }
    }
  }

  if (res)
    LOG_INFO(logger, logDepth, "***** Flipped {} incompatible edges in flat polygon", res);
  return res;
}

size_t cit::flipToCompatible(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipToCompatible");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipToCompatible"); });

  std::vector<Eigen::Matrix<size_t, 6, 1>> nFlipsArr;
  while (true) {
    // Sanity check about vertex degree
    checkVertexDegree(config.cointriA);
    checkVertexDegree(config.cointriB);
    auto scopeExit = kt84::make_ScopeExit([&config](){
      checkVertexDegree(config.cointriA);
      checkVertexDegree(config.cointriB);
    });

    updateCorrespondence(config, logger, logDepth + 1);

    nFlipsArr.push_back({});
    auto& nFlips = nFlipsArr.back();
    nFlips.setZero();
    if ((nFlips[0] = simpleFlip(config,           logger, logDepth + 1))) continue;
    if ((nFlips[2] = simpleCoFlip(config,         logger, logDepth + 1))) continue;
    if ((nFlips[3] = flipCompatible(config,       logger, logDepth + 1))) continue;
    if ((nFlips[4] = flipFlatPolygon(config,      logger, logDepth + 1))) continue;
    if ((nFlips[5] = doubleFlipCompatible(config, logger, logDepth + 1))) continue;
    break;
  }
  size_t res = 0;
  for (const auto& nFlips : nFlipsArr)
    res += nFlips.sum();

  config.cointriA.signpostTri->refreshQuantities();
  config.cointriB.signpostTri->refreshQuantities();

  return res;
}

size_t cit::flipPatchInteriorCompatibleEdgesToExploreVariations(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipPatchInteriorCompatibleEdgesToExploreVariations");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipPatchInteriorCompatibleEdgesToExploreVariations"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  // To not flip the same edge again
  std::unordered_set<Edge> flippedEdges;

  while (true) {
    // Collect flippable compatible edges inside incompatible patches
    bool flipped = false;
    for (PatchPair& pp : extractInconsistentPatches(config, false)) {
      Patch& A_patch = pp.first;
      analyzeInconsistentPatch(cointriA, A_patch);
      for (Edge A_e : A_patch.interiorEdges) {
        if (flippedEdges.count(A_e))
          continue;

        if (!isEdgeCompatible(cointriA, A_e))
          continue;
        Edge B_e = cointriA.correspondingEdge[A_e];

        if (!isEdgeFlippable(cointriA, A_e)) continue;
        if (!isEdgeFlippable(cointriB, B_e)) continue;

        // Check if flipping A_e/B_e preserves the number of compatible edges

        std::array<Vertex, 2> A_vOpp = getOppositeVertices(A_e);
        std::array<Vertex, 2> B_vOpp = getOppositeVertices(B_e);

        auto A_ue = getUniqueEdgeID(cointriA, A_vOpp);
        auto B_ue = getUniqueEdgeID(cointriB, B_vOpp);

        // When A_ue == B_ue, A_e/B_e have compatible edge diamonds, so they can be flipped

        // When A_ue != B_ue, A_e/B_e can still be flipped if there exists a corresponding incompatible edge with edge ID A_ue/B_ue
        if (A_ue != B_ue) {
          std::array<Vertex, 2> A_v = {
            cointriA.uniqueID_to_Vertex.at(B_ue.first),
            cointriA.uniqueID_to_Vertex.at(B_ue.second)
          };
          std::array<Vertex, 2> B_v = {
            cointriB.uniqueID_to_Vertex.at(A_ue.first),
            cointriB.uniqueID_to_Vertex.at(A_ue.second)
          };
          Halfedge A_he = A_v[0].connectingHalfedge(A_v[1]);
          Halfedge B_he = B_v[0].connectingHalfedge(B_v[1]);
          if (A_he == Halfedge()) continue;
          if (B_he == Halfedge()) continue;

          CIT_ASSERT(!isEdgeCompatible(cointriA, A_he.edge()));
          CIT_ASSERT(!isEdgeCompatible(cointriB, B_he.edge()));
        }

        cointriA.signpostTri->flipEdgeIfPossible(A_e);
        cointriB.signpostTri->flipEdgeIfPossible(B_e);

        LOG_DEBUG(logger, logDepth, "  Flipped compatible edge {}/{} in an incompatible patch (nV={}, nF={})", A_e, B_e, A_patch.vertices.size(), A_patch.faces.size());
        flippedEdges.insert(A_e);

        if (flipToCompatible(config, logger, logDepth + 1)) {
          LOG_DEBUG(logger, logDepth, "    Flipping did improve compatibility! Repeat the process");
          flipped = true;
          break;
        } else {
          LOG_DEBUG(logger, logDepth, "    Flipping didn't improve compatibility. Move on to the next edge");
        }
      }
      if (flipped)
        break;
    }
    if (!flipped)
      break;
  }

  if (!flippedEdges.empty())
    LOG_INFO(logger, logDepth, "****** Flipped {} compatible edges within incompatible patches to explore variations", flippedEdges.size());
  return flippedEdges.size();
}

size_t cit::convexifyConcaveIncompatiblePatch(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN convexifyConcaveIncompatiblePatch");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   convexifyConcaveIncompatiblePatch"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  auto convexifyConcaveIncompatiblePatch_inner = [&logger, &logDepth](CoInTri& cointri, const Patch& patch) -> std::set<int> {
    if (patch.boundaryHalfedges.size() > 1)
      return {}; // Complicated case with multiple boundary loops, skip

    std::set<int> res;
    for (size_t i = 0; i < patch.boundaryHalfedges[0].size(); ++i) {
      std::array<Halfedge, 2> he = {
        patch.boundaryHalfedges[0][(i + patch.boundaryHalfedges[0].size() - 1) % patch.boundaryHalfedges[0].size()],
        patch.boundaryHalfedges[0][i]
      };
      std::array<Vertex, 3> v = {
        he[0].vertex(),
        he[1].vertex(),
        he[1].tipVertex()
      };
      CIT_ASSERT(v[1] == he[0].tipVertex());
      std::array<int, 3> u;
      for (int j = 0; j < 3; ++j)
        u[j] = cointri.uniqueID_per_Vertex[v[j]];

      if (patch.vertexAngleSums.at(v[1]) < M_PI) continue;
      LOG_DEBUG(logger, logDepth, "Model {}'s inconsitent patch boundary vertex {} has large angle sum {}, trying to convexify...", cointri.mdata->name, u[1], patch.vertexAngleSums.at(v[1]));

      if (isVertexMovable(cointri, v[1])) {
        // Simply move v[1] to the midpoint of the line connecting v[0] & v[2]
        Vector2 p0 = cointri.signpostTri->halfedgeVectorsInVertex[he[0].twin()];
        Vector2 p2 = cointri.signpostTri->halfedgeVectorsInVertex[he[1]];
        Vector2 traceVec = 0.5 * (p0 + p2);
        traceVec += 0.001 * (p2 - p0) * Vector2{0, -1};   // Nudge outwards
        TraceGeodesicResult traceResult = traceGeodesic(*cointri.signpostTri, {v[1]}, traceVec);
        if (cointri.signpostTri->relocateInsertedVertex(v[1], traceResult.endPoint)) {
          LOG_DEBUG(logger, logDepth, "Was able to relocate the free vertex {}", u[1]);
          res.insert(u[1]);
        }
      }
    }
    return res;
  };

  std::vector<PatchPair> patchPairs = extractInconsistentPatches(config, false);
  for (PatchPair& pp: patchPairs) {
    Patch& A_patch = pp.first;
    Patch& B_patch = pp.second;
    analyzeInconsistentPatch(cointriA, A_patch);
    analyzeInconsistentPatch(cointriB, B_patch);
  }

  std::vector<std::set<int>> A_res;
  std::vector<std::set<int>> B_res;

  for (const PatchPair& pp : patchPairs) {
    Patch A_patch, B_patch;
    std::tie(A_patch, B_patch) = pp;
    if (!isPatchSizeConsistent(A_patch, B_patch)) continue;
    std::set<int> A_res_inner = convexifyConcaveIncompatiblePatch_inner(cointriA, A_patch);
    std::set<int> B_res_inner = convexifyConcaveIncompatiblePatch_inner(cointriB, B_patch);
    if (!A_res_inner.empty()) A_res.push_back(A_res_inner);
    if (!B_res_inner.empty()) B_res.push_back(B_res_inner);
  }

  size_t res = 0;
  if (!A_res.empty()) {
    LOG_INFO(logger, logDepth, "===== Convexified {} incompatible patches on A =====", A_res.size());
    for (size_t i = 0; i < A_res.size(); ++i) {
      std::ostringstream oss;
      for (int u : A_res[i])
        oss << u << " ";
      LOG_INFO(logger, logDepth, "======= processed concave vertices in patch #{}: {} =======", i, oss.str());
      res += A_res[i].size();
    }
  }
  if (!B_res.empty()) {
    LOG_INFO(logger, logDepth, "===== Convexified {} incompatible patches on B =====", B_res.size());
    for (size_t i = 0; i < B_res.size(); ++i) {
      std::ostringstream oss;
      for (int u : B_res[i])
        oss << u << " ";
      LOG_INFO(logger, logDepth, "======= processed concave vertices in patch #{}: {} =======", i, oss.str());
      res += B_res[i].size();
    }
  }
  if (res)
    LOG_INFO(logger, logDepth, "======= processed {} concave vertices in total =======", res);

  return res;
}

size_t cit::mergeInconsistentVertexPairs(Configuration& config, bool useVertexAdjacency, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN mergeInconsistentVertexPairs");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   mergeInconsistentVertexPairs"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  size_t res = 0;
  std::set<std::pair<int, int>> uniqueIDs;

  std::vector<PatchPair> patchPairs = extractInconsistentPatches(config, useVertexAdjacency);

  // bunch of template matchings
  auto mergeInconsistentVertexPairs_template1 = [&useVertexAdjacency, &cointriA](const Patch& A_patch, const Patch& B_patch, std::vector<std::vector<Halfedge>>& A_heCollapse) -> bool {
    if (useVertexAdjacency) return false;

    A_heCollapse.clear();
    A_heCollapse.resize(1);

    if (!A_patch.interiorVertices.empty()) return false;

    std::vector<Vertex> A_vConcave;
    for (Vertex A_v : A_patch.boundaryVertices) {
      if (A_patch.vertexAngleSums.at(A_v) >= M_PI)
        A_vConcave.push_back(A_v);
    }
    std::vector<Vertex> B_vConcave;
    for (Vertex B_v : B_patch.boundaryVertices) {
      if (B_patch.vertexAngleSums.at(B_v) >= M_PI)
        B_vConcave.push_back(B_v);
    }
    if (A_vConcave.size() != 1) return false;
    if (B_vConcave.size() != 1) return false;

    for (Halfedge he : A_vConcave[0].outgoingHalfedges()) {
      if (cointriA.correspondingVertex[he.tipVertex()] == B_vConcave[0] && A_patch.boundaryEdges.count(he.edge())) {
        A_heCollapse[0] = {he};
        return true;
      }
    }
    return false;
  };
  auto mergeInconsistentVertexPairs_template2 = [&useVertexAdjacency, &cointriA, &cointriB](const Patch& A_patch, const Patch& B_patch, std::vector<std::vector<Halfedge>>& A_heCollapse) -> bool {
    if (useVertexAdjacency) return false;

    A_heCollapse.clear();
    A_heCollapse.resize(1);

    if (A_patch.interiorVertices.size() != 1) return false;

    Vertex A_vInterior = *A_patch.interiorVertices.begin();

    Face A_fReversed;
    for (Face A_f : A_patch.faces) {
      if (isFaceReversed(cointriA, A_f)) {
        if (A_fReversed.getMesh())
          return false;
        A_fReversed = A_f;
      }
    }
    if (!A_fReversed.getMesh()) return false;

    std::vector<Halfedge> B_heCollapse;
    for (Halfedge A_he : A_fReversed.adjacentHalfedges()) {
      if (A_he.vertex() == A_vInterior || A_he.tipVertex() == A_vInterior) {
        A_heCollapse[0].push_back(A_he);
        B_heCollapse.push_back(cointriA.correspondingHalfedge[A_he]);
      }
    }
    if (A_heCollapse[0].size() != 2) return false;

    // prefer one of the two whose relative edge length sum is smaller
    std::array<double, 2> A_edgeLen;
    std::array<double, 2> B_edgeLen;
    for (size_t i = 0; i < 2; ++i) {
      A_edgeLen[i] = cointriA.signpostTri->edgeLengths[A_heCollapse[0][i].edge()];
      B_edgeLen[i] = cointriB.signpostTri->edgeLengths[B_heCollapse[i].edge()];
    }
    double A_edgeLenSum = A_edgeLen[0] + A_edgeLen[1];
    double B_edgeLenSum = B_edgeLen[0] + B_edgeLen[1];

    if (A_edgeLen[0]/A_edgeLenSum + B_edgeLen[0]/B_edgeLenSum > A_edgeLen[1]/A_edgeLenSum + B_edgeLen[1]/B_edgeLenSum)
      std::swap(A_heCollapse[0][0], A_heCollapse[0][1]);

    return true;
  };
  auto mergeInconsistentVertexPairs_template3 = [&useVertexAdjacency, &config, &cointriA, &cointriB, &logger, &logDepth](const Patch& A_patch, const Patch& B_patch, std::vector<std::vector<Halfedge>>& A_heCollapse) -> bool {
    if (useVertexAdjacency) return false;

    A_heCollapse.clear();

    if (A_patch.interiorVertices.size() % 2 != 0) return false;

    for (Edge A_e : A_patch.interiorEdges) {
      if (!isEdgeCompatible(cointriA, A_e)) continue;
      Edge B_e = cointriA.correspondingEdge[A_e];

      // consider interior edge only when both of its adjacent faces are reversed
      if (!isFaceReversed(cointriA, A_e.halfedge().face()) || !isFaceReversed(cointriA, A_e.halfedge().twin().face()))
        continue;

      // if both of the endpoints are interior vertices, then collapse it
      std::array<Vertex, 2> A_v = A_e.adjacentVertices();
      std::array<Vertex, 2> B_v = B_e.adjacentVertices();
      if (A_patch.interiorVertices.count(A_v[0]) && A_patch.interiorVertices.count(A_v[1]) &&
          B_patch.interiorVertices.count(B_v[0]) && B_patch.interiorVertices.count(B_v[1]))
      {
        A_heCollapse.push_back({A_e.halfedge()});
        continue;
      }

      // if the edge is flippable, and its opposite vertices are compatible and are both interior, then flip and collapse it
      if (!isEdgeFlippable(cointriA, A_e)) continue;
      if (!isEdgeFlippable(cointriB, B_e)) continue;

      std::array<Vertex, 2> A_vOpp = getOppositeVertices(A_e);
      std::array<Vertex, 2> B_vOpp = getOppositeVertices(B_e);

      if (getUniqueEdgeID(cointriA, A_vOpp) != getUniqueEdgeID(cointriB, B_vOpp))
        continue;

      if (A_patch.interiorVertices.count(A_vOpp[0]) && A_patch.interiorVertices.count(A_vOpp[1]) &&
          B_patch.interiorVertices.count(B_vOpp[0]) && B_patch.interiorVertices.count(B_vOpp[1]))
      {
        cointriA.signpostTri->flipEdgeIfPossible(A_e);
        cointriB.signpostTri->flipEdgeIfPossible(B_e);
        updateCorrespondence(config, logger, logDepth + 1);
        A_heCollapse.push_back({A_e.halfedge()});
      }
    }
    return !A_heCollapse.empty();
  };
  auto mergeInconsistentVertexPairs_template4 = [&useVertexAdjacency](const Patch& A_patch, const Patch& B_patch, std::vector<std::vector<Halfedge>>& A_heCollapse) -> bool {
    if (!useVertexAdjacency) return false;

    A_heCollapse.clear();
    A_heCollapse.resize(1);

    if (A_patch.faces.size() != 4) return false;
    if (A_patch.interiorEdges.size() != 3) return false;
    if (A_patch.boundaryEdges.size() != 6) return false;
    if (A_patch.interiorVertices.size() != 1) return false;
    if (A_patch.boundaryVertices.size() != 5) return false;

    Vertex A_vInterior = *A_patch.interiorVertices.begin();
    CIT_ASSERT(A_vInterior.degree() == 3);

    Vertex A_vHourglass;
    for (Vertex A_v : A_patch.boundaryVertices) {
      int count = 0;
      for (Edge A_e : A_v.adjacentEdges()) {
        count += A_patch.boundaryEdges.count(A_e);
      }
      if (count == 4)
        A_vHourglass = A_v;
    }
    CIT_ASSERT(A_vHourglass.getMesh());

    for (Halfedge A_heOut : A_vInterior.outgoingHalfedges()) {
      if (A_heOut.tipVertex() == A_vHourglass) {
        A_heCollapse[0] = {A_heOut};
        return true;
      }
    }
    return false;
  };
  auto mergeInconsistentVertexPairs_template5 = [&useVertexAdjacency, &config, &cointriA, &cointriB, &logger, &logDepth](const Patch& A_patch, const Patch& B_patch, std::vector<std::vector<Halfedge>>& A_heCollapse) -> bool {
    if (useVertexAdjacency) return false;

    A_heCollapse.clear();
    A_heCollapse.resize(1);

    if (A_patch.interiorVertices.size() != 0) return false;

    // Find a compatible interior edge
    for (Edge A_e : A_patch.interiorEdges) {
      Edge B_e;
      if (isEdgeCompatible(cointriA, A_e)) {
        // Easy case, the edge is already compatible
        B_e = cointriA.correspondingEdge[A_e];

      } else {
        if (isEdgeFlippable(cointriA, A_e)) {
          // For flippable A_e, see if flipping makes it compatible
          auto ue = getUniqueEdgeID(cointriA, getOppositeVertices(A_e));
          for (Edge B_e2 : B_patch.interiorEdges) {
            if (ue == getUniqueEdgeID(cointriB, B_e2)) {
              CIT_ASSERT(!isEdgeCompatible(cointriB, B_e2));
              B_e = B_e2;
              cointriA.signpostTri->flipEdgeIfPossible(A_e);
              updateCorrespondence(config, logger, logDepth + 1);
              break;
            }
          }
        }
        if (!B_e.getMesh()) {
          // Find fippable B_e that becomes compatible when flipped
          auto ue = getUniqueEdgeID(cointriA, A_e);
          for (Edge B_e2 : B_patch.interiorEdges) {
            if (!isEdgeCompatible(cointriB, B_e2) && isEdgeFlippable(cointriB, B_e2) && ue == getUniqueEdgeID(cointriB, getOppositeVertices(B_e2))) {
              B_e = B_e2;
              cointriB.signpostTri->flipEdgeIfPossible(B_e);
              updateCorrespondence(config, logger, logDepth + 1);
              break;
            }
          }
        }
      }

      if (!B_e.getMesh()) continue;

      // the edge to be collapsed must have concave vertices at both of its endpoints, for both A and B
      std::array<Vertex, 2> A_v = A_e.adjacentVertices();
      std::array<Vertex, 2> B_v = B_e.adjacentVertices();
      if (A_patch.vertexAngleSums.at(A_v[0]) < M_PI) continue;
      if (A_patch.vertexAngleSums.at(A_v[1]) < M_PI) continue;
      if (B_patch.vertexAngleSums.at(B_v[0]) < M_PI) continue;
      if (B_patch.vertexAngleSums.at(B_v[1]) < M_PI) continue;

      // if there is more than one such edges, then bad
      if (!A_heCollapse[0].empty()) return false;

      A_heCollapse[0] = {A_e.halfedge()};
    }

    return !A_heCollapse[0].empty();
  };

  // apply templates to each of patchPairs
  for (size_t k = 0; k < patchPairs.size(); ++k) {
    PatchPair& pp = patchPairs[k];
    Patch& A_patch = pp.first;
    Patch& B_patch = pp.second;
    analyzeInconsistentPatch(cointriA, A_patch);
    analyzeInconsistentPatch(cointriB, B_patch);

    if (!isPatchSizeConsistent(A_patch, B_patch)) continue;

    if (!useVertexAdjacency)
      if (A_patch.boundaryEdges.size() != A_patch.boundaryVertices.size()) continue;

    // check if the interior vertices for A & B have the same set of unique IDs
    std::set<int> A_interiorVertexUniqueIDs;
    std::set<int> B_interiorVertexUniqueIDs;
    for (Vertex A_v : A_patch.interiorVertices) A_interiorVertexUniqueIDs.insert(cointriA.uniqueID_per_Vertex[A_v]);
    for (Vertex B_v : B_patch.interiorVertices) B_interiorVertexUniqueIDs.insert(cointriB.uniqueID_per_Vertex[B_v]);
      if (A_interiorVertexUniqueIDs != B_interiorVertexUniqueIDs) continue;

    std::vector<std::vector<Halfedge>> A_heCollapse;    // outer loop represents multiple edges to collapse, inner loop represents preference over which edge to collapse
    if (mergeInconsistentVertexPairs_template1(A_patch, B_patch, A_heCollapse) ||
        mergeInconsistentVertexPairs_template2(A_patch, B_patch, A_heCollapse) ||
        mergeInconsistentVertexPairs_template3(A_patch, B_patch, A_heCollapse) ||
        mergeInconsistentVertexPairs_template4(A_patch, B_patch, A_heCollapse) ||
        mergeInconsistentVertexPairs_template5(A_patch, B_patch, A_heCollapse))
    {
      CIT_ASSERT(!A_heCollapse.empty());
      for (size_t i = 0; i < A_heCollapse.size(); ++i) {
        CIT_ASSERT(!A_heCollapse[i].empty());
        for (Halfedge A_he : A_heCollapse[i]) {
          CIT_ASSERT(A_he.getMesh());

          if (A_he.isDead()) continue;

          CIT_ASSERT(isEdgeCompatible(cointriA, A_he.edge()));

          Halfedge B_he = cointriA.correspondingHalfedge[A_he].twin();
          CIT_ASSERT(!B_he.isDead());

          std::array<Vertex, 2> A_v = { A_he.vertex(), A_he.twin().vertex() };
          std::array<Vertex, 2> B_v = { B_he.vertex(), B_he.twin().vertex() };

          // if either of the vertices already went through another merge, skip
          if (isVertexMerged(cointriA.uniqueID_per_Vertex[A_v[0]])) continue;
          if (isVertexMerged(cointriA.uniqueID_per_Vertex[A_v[1]])) continue;

          // Prevent degree-two vertex from being generated after merge
          // -- Check the central vertex
          if (A_v[0].degree() + A_v[1].degree() - 4 < 3) continue;
          if (B_v[0].degree() + B_v[1].degree() - 4 < 3) continue;
          // -- Check the opposite vertices
          std::array<Vertex, 2> A_vOpp = getOppositeVertices(A_he.edge());
          std::array<Vertex, 2> B_vOpp = getOppositeVertices(B_he.edge());
          if (A_vOpp[0].degree() < 4 || A_vOpp[1].degree() < 4) continue;
          if (B_vOpp[0].degree() < 4 || B_vOpp[1].degree() < 4) continue;

          // Edge is collapsible only when its one end is original and the other is inserted
          if ((A_v[0].getIndex() < mdataA.nV) + (A_v[1].getIndex() < mdataA.nV) != 1) continue;
          CIT_ASSERT((B_v[0].getIndex() < mdataB.nV) + (B_v[1].getIndex() < mdataB.nV) == 1);

          // Ensure the first vertex is fixed
          if (A_v[0].getIndex() >= mdataA.nV) {
            A_he = A_he.twin();
            B_he = B_he.twin();
            std::swap(A_v[0], A_v[1]);
            std::swap(B_v[0], B_v[1]);
          }
          CIT_ASSERT(A_v[0].getIndex() < mdataA.nV);
          CIT_ASSERT(A_v[1].getIndex() >= mdataA.nV);
          CIT_ASSERT(B_v[0].getIndex() < mdataB.nV);
          CIT_ASSERT(B_v[1].getIndex() >= mdataB.nV);

          // both of its adjacent faces cannot be compatible
          Face A_f0 = A_he.face();
          Face A_f1 = A_he.twin().face();
          // CIT_ASSERT(!isFaceCompatible(cointriA, A_f0) || !isFaceCompatible(cointriA, A_f1));    // TODO: figure out why this assumption is wrong

          // some preparation before merge:
          // -- identify unique vertex/edge ID before/after merge
          std::array<int, 2> A_uniqueVertexID = {
            cointriA.uniqueID_per_Vertex[A_v[0]],
            cointriA.uniqueID_per_Vertex[A_v[1]],
          };
          int uniqueVertexID_new;
          if (mdataA.nV < mdataB.nV) {
            CIT_ASSERT(A_uniqueVertexID[0] == (int)(A_v[0].getIndex() + 1));
            CIT_ASSERT(A_uniqueVertexID[1] == -(int)(B_v[0].getIndex() + 1));
            uniqueVertexID_new = A_v[0].getIndex() + 1 + (mdataA.nV + 1) * (B_v[0].getIndex() + 1);
          } else {
            CIT_ASSERT(A_uniqueVertexID[0] == -(int)(A_v[0].getIndex() + 1));
            CIT_ASSERT(A_uniqueVertexID[1] == (int)(B_v[0].getIndex() + 1));
            uniqueVertexID_new = B_v[0].getIndex() + 1 + (mdataB.nV + 1) * (A_v[0].getIndex() + 1);
          }

          // perform edge collapse if possible
          bool A_collapsible = cointriA.signpostTri->collapseInteriorEdge(A_he, true);
          bool B_collapsible = cointriB.signpostTri->collapseInteriorEdge(B_he, true);
          if (!A_collapsible || !B_collapsible)
            continue;

          uniqueIDs.insert({cointriA.uniqueID_per_Vertex[A_he.vertex()], cointriA.uniqueID_per_Vertex[A_he.tipVertex()]});

          cointriA.signpostTri->collapseInteriorEdge(A_he, false);
          cointriB.signpostTri->collapseInteriorEdge(B_he, false);

          // assign new unique vertex ID (representing the overlap of AB vertices) to the collapsed vertex
          cointriA.uniqueID_per_Vertex[A_v[0]] = uniqueVertexID_new;
          cointriB.uniqueID_per_Vertex[B_v[0]] = uniqueVertexID_new;

          updateCorrespondence(config, logger, logDepth + 1);

          ++res;
          break;
        }
      }
    }
  }

  if (res) {
    LOG_INFO(logger, logDepth, "===== Merged {} vertex pairs in inconsistent positions =====", res);
    for (std::pair<int, int> u : uniqueIDs)
      LOG_INFO(logger, logDepth, "======= processed vertex pair: ({}, {}) =======", u.first, u.second);
  }
  return res;
}

bool cit::splitMergedVertices(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN splitMergedVertices");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   splitMergedVertices"); });

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  CIT_ASSERT(config.nCompatibleEdges == cointriA.signpostTri->intrinsicMesh->nEdges());

  std::set<std::array<int, 3>> res;

  for (Vertex A_vOld : cointriA.signpostTri->intrinsicMesh->vertices()) {
    int uniqueVertexID = cointriA.uniqueID_per_Vertex[A_vOld];
    if (!isVertexMerged(uniqueVertexID)) continue;;

    Vertex B_vOld = cointriA.correspondingVertex[A_vOld];

    if (mdataA.anchorVertexIDs.count(A_vOld.getIndex())) {
      CIT_ASSERT(mdataB.anchorVertexIDs.count(B_vOld.getIndex()));
      continue;
    }

    LOG_DEBUG(logger, logDepth, "====== Processing merged vertex {} ({} in A:, {} in B) ... ======", uniqueVertexID, A_vOld.getIndex(), B_vOld.getIndex());

    // unique ID for the new vertices
    const size_t A_nV = mdataA.nV;
    const size_t B_nV = mdataB.nV;
    int A_u, B_u;
    if (A_nV < B_nV) {
      A_u = uniqueVertexID % (A_nV + 1);
      B_u = uniqueVertexID / (A_nV + 1);
      B_u *= -1;
    } else {
      B_u = uniqueVertexID % (B_nV + 1);
      A_u = uniqueVertexID / (B_nV + 1);
      A_u *= -1;
    }

    CIT_ASSERT(!cointriA.uniqueID_to_Vertex.count(A_u));
    CIT_ASSERT(!cointriA.uniqueID_to_Vertex.count(B_u));
    CIT_ASSERT(!cointriB.uniqueID_to_Vertex.count(A_u));
    CIT_ASSERT(!cointriB.uniqueID_to_Vertex.count(B_u));

    res.insert({uniqueVertexID, A_u, B_u});

    // find a pair of edges around the vertex whose angle difference is closest to 180 degrees
    std::array<Halfedge, 2> A_halfedgePair;
    std::array<Halfedge, 2> B_halfedgePair;
    std::set<std::pair<int, int>> alreadyChecked;
    double matchingCost = std::numeric_limits<double>::infinity();
    for (Halfedge A_h0 : A_vOld.outgoingHalfedges()) {
      for (Halfedge A_h1 : A_vOld.outgoingHalfedges()) {
        if (A_h0 == A_h1)
          continue;
        std::pair<int, int> p{A_h0.getIndex(), A_h1.getIndex()};
        if (p.first > p.second)
          std::swap(p.first, p.second);
        if (alreadyChecked.count(p))
          continue;
        alreadyChecked.insert(p);

        double A_angle0 = cointriA.signpostTri->halfedgeVectorsInVertex[A_h0].arg();
        double A_angle1 = cointriA.signpostTri->halfedgeVectorsInVertex[A_h1].arg();
        double A_matchingCost = M_PI - angleAbsDifference(A_angle0, A_angle1);
        CIT_ASSERT(A_matchingCost >= 0);

        Halfedge B_h0 = cointriA.correspondingHalfedge[A_h0];
        Halfedge B_h1 = cointriA.correspondingHalfedge[A_h1];
        CIT_ASSERT(B_h0.vertex() == B_vOld);
        CIT_ASSERT(B_h1.vertex() == B_vOld);

        double B_angle0 = cointriB.signpostTri->halfedgeVectorsInVertex[B_h0].arg();
        double B_angle1 = cointriB.signpostTri->halfedgeVectorsInVertex[B_h1].arg();
        double B_matchingCost = M_PI - angleAbsDifference(B_angle0, B_angle1);
        CIT_ASSERT(B_matchingCost >= 0);

        if (A_matchingCost + B_matchingCost < matchingCost) {
          matchingCost = A_matchingCost + B_matchingCost;
          A_halfedgePair = {A_h0, A_h1};
          B_halfedgePair = {B_h0, B_h1};
        }
      }
    }

    // Avoid going through the exact middle of the two halfedge directions, as this may coinside with an existing edge if the input mesh is very regular
    const double angleNudge = 0.05;

    // Perform split on A; new vertex corresponds to B's vertex inserted on A
    std::array<double, 2> A_anglePair;
    for (int i = 0; i < 2; ++i) {
      A_anglePair[i] = cointriA.signpostTri->halfedgeVectorsInVertex[A_halfedgePair[i]].arg();
    }
    if (A_anglePair[0] > A_anglePair[1]) {
      A_anglePair[1] += 2. * M_PI;
    }
    double A_traceVecAngle = std::fmod((0.5 + angleNudge) * A_anglePair[0] + (0.5 - angleNudge) * A_anglePair[1], 2. * M_PI);
    double A_traceVecLen = std::numeric_limits<double>::infinity();
    for (Halfedge A_he : A_vOld.outgoingHalfedges()) {
      double triangleHeight = 2. * cointriA.signpostTri->area(A_he.face()) / cointriA.signpostTri->intrinsicEdgeLengths[A_he.next().edge()];
      A_traceVecLen = std::min<double>(A_traceVecLen, triangleHeight);
    }
    A_traceVecLen *= 0.1;
    Vector2 A_traceVec = Vector2::fromAngle(A_traceVecAngle) * A_traceVecLen;
    Halfedge A_heNew;
    try {
      A_heNew = cointriA.signpostTri->splitVertexAlongTwoEdges(A_halfedgePair[0], A_halfedgePair[1], A_traceVec);
    } catch (const std::runtime_error& e) {
      LOG_DEBUG(logger, logDepth + 1, "Failed to perform split on A: what = {}", e.what());
      return false;
    }
    Vertex A_vNew = A_heNew.tipVertex();
    cointriA.uniqueID_per_Vertex[A_vOld] = A_u;
    cointriA.uniqueID_per_Vertex[A_vNew] = B_u;

    // Perform split on B; new vertex corresponds to A's vertex inserted on B
    std::array<double, 2> B_anglePair;
    for (int i = 0; i < 2; ++i) {
      B_anglePair[i] = cointriB.signpostTri->halfedgeVectorsInVertex[B_halfedgePair[i]].arg();
    }
    if (B_anglePair[0] > B_anglePair[1]) {
      B_anglePair[1] += 2. * M_PI;
    }
    double B_traceVecAngle = std::fmod((0.5 + angleNudge) * B_anglePair[0] + (0.5 - angleNudge) * B_anglePair[1] + M_PI, 2. * M_PI);  // note the opposite direction
    double B_traceVecLen = std::numeric_limits<double>::infinity();
    for (Halfedge B_he : B_vOld.outgoingHalfedges()) {
      double triangleHeight = 2. * cointriB.signpostTri->area(B_he.face()) / cointriB.signpostTri->intrinsicEdgeLengths[B_he.next().edge()];
      B_traceVecLen = std::min<double>(B_traceVecLen, triangleHeight);
    }
    B_traceVecLen *= 0.1;
    Vector2 B_traceVec = Vector2::fromAngle(B_traceVecAngle) * B_traceVecLen;
    Halfedge B_heNew;
    try {
      B_heNew = cointriB.signpostTri->splitVertexAlongTwoEdges(B_halfedgePair[0], B_halfedgePair[1], B_traceVec);
    } catch (const std::runtime_error& e) {
      LOG_DEBUG(logger, logDepth + 1, "Failed to perform split on B: what = {}", e.what());
      return false;
    }
    Vertex B_vNew = B_heNew.tipVertex();
    cointriB.uniqueID_per_Vertex[B_vOld] = B_u;
    cointriB.uniqueID_per_Vertex[B_vNew] = A_u;

    updateCorrespondence(config, logger, logDepth + 1);
  }

  if (res.size()) {
    LOG_INFO(logger, logDepth, "===== Split {} merged vertices =====", res.size());
    for (std::array<int, 3> r : res)
      LOG_INFO(logger, logDepth, "======= Merged vertex {} split into {} & {} =======", r[0], r[1], r[2]);
  }
  return true;
}

cit::OverlayVertex cit::getOverlayVertexAtEdgeIntersection(const Configuration& config, std::array<Edge, 2> e) {
  OverlayVertex overlayVertex;
  for (int i = 0; i < 2; ++i) {
    SurfaceMesh* mesh = e[i].getMesh();
    if (mesh == mdataA.mesh.get()) {
      overlayVertex.eA = e[i];
    } else if (mesh == mdataB.mesh.get()) {
      overlayVertex.eB = e[i];
    } else {
      // Convert B_eIntrinsic to A_eIntrinsic, following convention (see comment in OverlayVertex definition)
      if (mesh == config.cointriB.signpostTri->intrinsicMesh.get())
        e[i] = config.cointriB.correspondingEdge[e[i]];
      CIT_ASSERT(e[i].getMesh() == config.cointriA.signpostTri->intrinsicMesh.get());
      overlayVertex.eIntrinsic = e[i];
    }
  }
  return overlayVertex;
}

void cit::storeOverlayWedgesAtEdgeIntersection(Configuration& config, std::array<SurfacePoint, 2> edgePoint, size_t logDepth) {
#if 0 // This clutters the log a lot
  DLOG_TRACE(logDepth, "BEGIN storeOverlayWedgesAtEdgeIntersection");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   storeOverlayWedgesAtEdgeIntersection"); });
#endif

  CIT_ASSERT(edgePoint[0].type == SurfacePointType::Edge);
  CIT_ASSERT(edgePoint[1].type == SurfacePointType::Edge);

  const CoInTri& cointriA = config.cointriA;
  const CoInTri& cointriB = config.cointriB;
  ManifoldSurfaceMesh* inputMeshA = mdataA.mesh.get();
  ManifoldSurfaceMesh* inputMeshB = mdataB.mesh.get();
  ManifoldSurfaceMesh* intrinsicMeshA = cointriA.signpostTri->intrinsicMesh.get();
  ManifoldSurfaceMesh* intrinsicMeshB = cointriB.signpostTri->intrinsicMesh.get();

  // Assume edge[0] is from an input mesh (A or B)
  CIT_ASSERT(edgePoint[0].getMesh() == inputMeshA || edgePoint[0].getMesh() == inputMeshB);

  // >>>>> Ensure edge[1] is oriented to the left relative to edge[0] >>>>>
  {
    const ModelData& mdataP = edgePoint[0].getMesh() == inputMeshA ? mdataA : mdataB;
    const ModelData& mdataQ = edgePoint[0].getMesh() == inputMeshA ? mdataB : mdataA;
    const CoInTri& cointriP = &mdataP == &mdataA ? cointriA : cointriB;
    const CoInTri& cointriQ = &mdataP == &mdataA ? cointriB : cointriA;
    ManifoldSurfaceMesh* inputMeshP = mdataP.mesh.get();
    ManifoldSurfaceMesh* inputMeshQ = mdataQ.mesh.get();
    ManifoldSurfaceMesh* intrinsicMeshP = &mdataP == &mdataA ? intrinsicMeshA : intrinsicMeshB;
    ManifoldSurfaceMesh* intrinsicMeshQ = &mdataP == &mdataA ? intrinsicMeshB : intrinsicMeshA;

    // Assume given input edge point is never preserved in the intrinsic mesh
    CIT_ASSERT(!cointriP.signpostTri->isInputEdgePointPreserved(edgePoint[0]));
    if (edgePoint[1].getMesh() == inputMeshQ)
      CIT_ASSERT(!cointriQ.signpostTri->isInputEdgePointPreserved(edgePoint[1]));

    // Look at the edge path on P corresponding to edgePoint[1].edge
    CIT_ASSERT(edgePoint[1].getMesh() != inputMeshP);
    CIT_ASSERT(edgePoint[1].getMesh() != intrinsicMeshQ);
    const std::vector<SurfacePoint>& P_edgePath = edgePoint[1].getMesh() == intrinsicMeshP
      ? cointriP.intrinsicEdgePath[edgePoint[1].edge]
      : cointriQ.inputEdgePathOnOtherInput[edgePoint[1].edge];

    // On the edge path, find one previous to edgePoint[0]
    // P_edgePath can have multiple intersections with edgePoint[0].edge if both edges are from input meshes, so select the closest match
    kt84::MinSelector<decltype(P_edgePath.end())> found(P_edgePath.end());
    for (auto sp = P_edgePath.begin(); sp != P_edgePath.end(); ++sp) {
      if (sp->edge == edgePoint[0].edge)
        found.update(std::fabs(sp->tEdge - edgePoint[0].tEdge), sp);
    }
    CIT_ASSERT(found.count);
    if (found.score >= 0.02) {    // Is such a large error like this really okay?
      DLOG_ERROR(logDepth, "  Edge point matching error exceeded threshold: min={}, count={}, average={}", found.score, found.count, found.average_score());
      CIT_ASSERT(false);
    }
    CIT_ASSERT(found.value != P_edgePath.end());
    CIT_ASSERT(found.value != P_edgePath.begin());      // The intersection should be somewhere inbetween

    SurfacePoint P_spPrev = *--found.value;

    Face P_fInput = sharedFace(edgePoint[0], P_spPrev);
    CIT_ASSERT(P_fInput != Face());

    Vector3 bcCross = edgePoint[0].inFace(P_fInput).faceCoords;
    Vector3 bcFrom_edge0 = SurfacePoint(edgePoint[0].edge.halfedge().vertex()).inFace(P_fInput).faceCoords;
    Vector3 bcFrom_edge1 = P_spPrev.inFace(P_fInput).faceCoords;

    Vector2 d0 = Vector2::castFrom(bcCross - bcFrom_edge0);
    Vector2 d1 = Vector2::castFrom(bcCross - bcFrom_edge1);

    // Check the orientation condition, swap if applicable
    if (cross(d0, d1) < 0)
      std::swap(edgePoint[0], edgePoint[1]);
  }
  // <<<<< Ensure edge[1] is oriented to the left relative to edge[0] <<<<<

  // Identify which mesh each edge belongs to (A, B, or intrinsic)
  std::array<std::string, 2> edgeType;
  Edge A_eInput;
  Edge B_eInput;
  double A_tEdge = -1.;
  double B_tEdge = -1.;
  for (int i = 0; i < 2; ++i) {
    SurfaceMesh* mesh = edgePoint[i].getMesh();

    // Output of isIntrinsicEdgePartiallyOriginal
    double tEdgeMin, tEdgeMax;
    bool reversed;

    if (mesh == inputMeshA) {
      edgeType[i] = "A";

    } else if (mesh == inputMeshB) {
      edgeType[i] = "B";

    } else if (mesh == intrinsicMeshA && cointriA.signpostTri->isIntrinsicEdgePartiallyOriginal(edgePoint[i].edge, &A_eInput, &tEdgeMin, &tEdgeMax, &reversed)) {
      edgeType[i] = "A";
      A_tEdge = (1. - edgePoint[i].tEdge) * (reversed ? tEdgeMax : tEdgeMin) + edgePoint[i].tEdge * (reversed ? tEdgeMin : tEdgeMax);

    } else if (mesh == intrinsicMeshB && cointriB.signpostTri->isIntrinsicEdgePartiallyOriginal(edgePoint[i].edge, &B_eInput, &tEdgeMin, &tEdgeMax, &reversed)) {
      edgeType[i] = "B";
      B_tEdge = (1. - edgePoint[i].tEdge) * (reversed ? tEdgeMax : tEdgeMin) + edgePoint[i].tEdge * (reversed ? tEdgeMin : tEdgeMax);

    } else if (mesh == intrinsicMeshA && cointriB.signpostTri->isIntrinsicEdgePartiallyOriginal(cointriA.correspondingEdge[edgePoint[i].edge], &B_eInput, &tEdgeMin, &tEdgeMax, &reversed)) {
      edgeType[i] = "B";
      B_tEdge = (1. - edgePoint[i].tEdge) * (reversed ? tEdgeMax : tEdgeMin) + edgePoint[i].tEdge * (reversed ? tEdgeMin : tEdgeMax);

    } else if (mesh == intrinsicMeshB && cointriA.signpostTri->isIntrinsicEdgePartiallyOriginal(cointriB.correspondingEdge[edgePoint[i].edge], &A_eInput, &tEdgeMin, &tEdgeMax, &reversed)) {
      edgeType[i] = "A";
      A_tEdge = (1. - edgePoint[i].tEdge) * (reversed ? tEdgeMax : tEdgeMin) + edgePoint[i].tEdge * (reversed ? tEdgeMin : tEdgeMax);
    }
  }

  // For each edge, if it belongs to A or B, figure out face points on its two adjacent faces
  std::array<SurfacePoint, 2> facePointA;
  std::array<SurfacePoint, 2> facePointB;
  for (int i = 0; i < 2; ++i) {
    SurfaceMesh* mesh = edgePoint[i].getMesh();
    if (edgeType[i] == "A") {
      if (mesh == inputMeshA) {
        facePointA[0] = edgePoint[i].inFace(edgePoint[i].edge.halfedge().face());
        facePointA[1] = edgePoint[i].inFace(edgePoint[i].edge.halfedge().twin().face());
      } else {
        CIT_ASSERT(A_eInput.getMesh());
        CIT_ASSERT(A_tEdge != -1.);
        facePointA[0] = SurfacePoint{A_eInput, A_tEdge}.inFace(A_eInput.halfedge().face());
        facePointA[1] = SurfacePoint{A_eInput, A_tEdge}.inFace(A_eInput.halfedge().twin().face());
      }
    } else if (edgeType[i] == "B") {
      if (mesh == inputMeshB) {
        facePointB[0] = edgePoint[i].inFace(edgePoint[i].edge.halfedge().face());
        facePointB[1] = edgePoint[i].inFace(edgePoint[i].edge.halfedge().twin().face());
      } else {
        CIT_ASSERT(B_eInput.getMesh());
        CIT_ASSERT(B_tEdge != -1.);
        facePointB[0] = SurfacePoint{B_eInput, B_tEdge}.inFace(B_eInput.halfedge().face());
        facePointB[1] = SurfacePoint{B_eInput, B_tEdge}.inFace(B_eInput.halfedge().twin().face());
      }
    }
  }

  std::vector<OverlayWedge> wedges(4);

  // Set common overlay vertex
  OverlayVertex overlayVertex = getOverlayVertexAtEdgeIntersection(config, {edgePoint[0].edge, edgePoint[1].edge});

  // If this is an input-input intersection, set tEdgeA/tEdgeB field to make distinction between potentially more than one intersections
  if ((edgePoint[0].getMesh() == inputMeshA || edgePoint[0].getMesh() == inputMeshB) &&
      (edgePoint[1].getMesh() == inputMeshA || edgePoint[1].getMesh() == inputMeshB))
  {
    for (int i = 0; i < 2; ++i)
      if (edgePoint[i].getMesh() == inputMeshA)
        overlayVertex.tEdgeA = edgePoint[i].tEdge;
      else
        overlayVertex.tEdgeB = edgePoint[i].tEdge;
  }

  for (OverlayWedge& wedge : wedges)
    wedge.overlayVertex = overlayVertex;

  double tEdge0 = edgePoint[0].tEdge;
  double tEdge1 = edgePoint[1].tEdge;
  Halfedge he0 = edgePoint[0].edge.halfedge();
  Halfedge he1 = edgePoint[1].edge.halfedge();

  // Convert intrinsic halfedge on B to its corresponding one on A, following convention (see comment in OverlayVertex definition)
  if (he0.getMesh() == intrinsicMeshB) he0 = cointriB.correspondingHalfedge[he0];
  if (he1.getMesh() == intrinsicMeshB) he1 = cointriB.correspondingHalfedge[he1];

  // Fill in connectivity & location info appropriately
  /*
           *
           |
          ^||
      he0 ||v
  e0   --> | -->
  o--------x--------*
       <-- | <--
          ^||
      he1 ||v
           |
           o e1
  */

  wedges[0].halfedgeIn = he0;
  wedges[0].halfedgeOut = he1;
  wedges[0].tHalfedgeIn = tEdge0;
  wedges[0].tHalfedgeOut = tEdge1;
  if (edgeType[0] == "A") wedges[0].facePointA = facePointA[0];
  if (edgeType[0] == "B") wedges[0].facePointB = facePointB[0];
  if (edgeType[1] == "A") wedges[0].facePointA = facePointA[0];
  if (edgeType[1] == "B") wedges[0].facePointB = facePointB[0];

  wedges[1].halfedgeIn = he1;
  wedges[1].halfedgeOut = he0.twin();
  wedges[1].tHalfedgeIn = tEdge1;
  wedges[1].tHalfedgeOut = 1 - tEdge0;
  if (edgeType[0] == "A") wedges[1].facePointA = facePointA[1];
  if (edgeType[0] == "B") wedges[1].facePointB = facePointB[1];
  if (edgeType[1] == "A") wedges[1].facePointA = facePointA[0];
  if (edgeType[1] == "B") wedges[1].facePointB = facePointB[0];

  wedges[2].halfedgeIn = he0.twin();
  wedges[2].halfedgeOut = he1.twin();
  wedges[2].tHalfedgeIn = 1 - tEdge0;
  wedges[2].tHalfedgeOut = 1 - tEdge1;
  if (edgeType[0] == "A") wedges[2].facePointA = facePointA[1];
  if (edgeType[0] == "B") wedges[2].facePointB = facePointB[1];
  if (edgeType[1] == "A") wedges[2].facePointA = facePointA[1];
  if (edgeType[1] == "B") wedges[2].facePointB = facePointB[1];

  wedges[3].halfedgeIn = he1.twin();
  wedges[3].halfedgeOut = he0;
  wedges[3].tHalfedgeIn = 1 - tEdge1;
  wedges[3].tHalfedgeOut = tEdge0;
  if (edgeType[0] == "A") wedges[3].facePointA = facePointA[0];
  if (edgeType[0] == "B") wedges[3].facePointB = facePointB[0];
  if (edgeType[1] == "A") wedges[3].facePointA = facePointA[1];
  if (edgeType[1] == "B") wedges[3].facePointB = facePointB[1];

  for (size_t i = 0; i < 4; ++i)
    CIT_ASSERT(wedges[i].tHalfedgeIn > 0.);

  config.overlayWedgesPerOverlayVertex[overlayVertex] = wedges;
}

// Internal function needed by computeEdgePath_inner2
namespace cit {
namespace detail {
void getNextTraceVec(const CoInTri& cointri, const std::array<SurfacePoint, 3> sp, double& traceVecAngleDelta, double& traceVecLen) {
  CIT_ASSERT(sp[1].type == SurfacePointType::Edge);

  std::array<Vector2, 4> Q_diamondP = getLayoutDiamond(*cointri.signpostTri, sp[1].edge.halfedge());
  std::array<Halfedge,4> Q_diamondHE = sp[1].edge.diamondBoundary();
  std::array<Vector2, 3> p;

  p[1] = (1 - sp[1].tEdge) * Q_diamondP[2] + sp[1].tEdge * Q_diamondP[0];

  for (int i : {0, 2}) {
    if (sp[i].type == SurfacePointType::Vertex) {
      // the vertex should be one of the opposite vertices of sp[1].edge
      int found = 1;
      if (Q_diamondHE[found].vertex() != sp[i].vertex)
        found = 3;
      CIT_ASSERT(Q_diamondHE[found].vertex() == sp[i].vertex);

      p[i] = Q_diamondP[found];

    } else {
      CIT_ASSERT(sp[i].type == SurfacePointType::Edge);

      // find where the edge is in the diamond, also adjust tEdge depending on edge orientation
      int found = -1;
      double t = sp[i].tEdge;
      for (int j = 0; j < 4; ++j) {
        if (Q_diamondHE[j].edge() == sp[i].edge) {
          found = j;
          if (Q_diamondHE[j] != sp[i].edge.halfedge()) {
            CIT_ASSERT(Q_diamondHE[j] == sp[i].edge.halfedge().twin());
            t = 1 - t;
          }
          break;
        }
      }
      CIT_ASSERT(found >= 0);
      p[i] = (1 - t) * Q_diamondP[found] + t * Q_diamondP[(found + 1) % 4];
    }
  }

  traceVecAngleDelta = (p[2] - p[1]).arg() - (p[1] - p[0]).arg();
  traceVecLen = norm(p[2] - p[1]);
};
}   // namespace detail
}   // namespace cit

void cit::computeEdgePath(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeEdgePath");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeEdgePath"); });

  config.overlayWedgesPerOverlayVertex = {};

  // Step 1: trace intrinsic edge over input mesh
  auto computeEdgePath_inner1 = [&logDepth](CoInTri& cointri) {
    DLOG_INFO(logDepth, "Tracing intrinsic edges on {} ... ", cointri.mdata->name);
    cointri.intrinsicEdgePath = cointri.signpostTri->traceEdges();
    DLOG_INFO(logDepth, "  Done");
  };
  computeEdgePath_inner1(config.cointriA);
  computeEdgePath_inner1(config.cointriB);

  // Step 2: rearrange intrinsic edge's path points on input mesh to get input edge's path points on intrinsic mesh
  auto computeEdgePath_inner2 = [&config, &logDepth](CoInTri& cointri) {
    DLOG_INFO(logDepth, "Rearranging traced path points to fill input edge path on {} ...", cointri.mdata->name);

    // use std::map to sort points on each input edge intersected by intrinsic edges
    EdgeData<std::map<double, SurfacePoint>> inputEdgePathOnIntrinsic_map(cointri.signpostTri->inputMesh);

    cointri.edgePointMap_input2intrinsic.clear();
    cointri.edgePointMap_intrinsic2input.clear();

    // First, register each input edge's endpoints to the sorted list
    for (Edge inputE : cointri.signpostTri->inputMesh.edges()) {
      std::array<Vertex, 2> inputV = inputE.adjacentVertices();
      inputEdgePathOnIntrinsic_map[inputE][0.] = SurfacePoint(cointri.signpostTri->intrinsicMesh->vertex(inputV[0].getIndex()));
      inputEdgePathOnIntrinsic_map[inputE][1.] = SurfacePoint(cointri.signpostTri->intrinsicMesh->vertex(inputV[1].getIndex()));
    }

    // Also, register all the intrinsic vertices whose locations are edge points
    for (Vertex intrinsicV : cointri.signpostTri->intrinsicMesh->vertices()) {
      SurfacePoint inputSP = cointri.signpostTri->vertexLocations[intrinsicV];
      if (inputSP.type == SurfacePointType::Edge)
        inputEdgePathOnIntrinsic_map[inputSP.edge][inputSP.tEdge] = SurfacePoint(intrinsicV);
    }

    // Then, process each intrinsic edge
    for (Edge intrinsicE : cointri.signpostTri->intrinsicMesh->edges()) {
      // If this intrinsic edge is partially original, skip it because both of its endpoints have already been registered above
      if (cointri.signpostTri->isIntrinsicEdgePartiallyOriginal(intrinsicE))
        continue;

      // The usual case where the intrinsic edge is non-original, possibly intersecting with input edges
      const std::vector<SurfacePoint>& inputSPs = cointri.intrinsicEdgePath[intrinsicE];
      const size_t n = inputSPs.size();
      CIT_ASSERT(n >= 2);

      // determine parameter values t[i] along intrinsicE
      std::vector<double> t(n);
      std::vector<Vector3> p(n);
      for (size_t i = 0; i < n; ++i) {
        p[i] = inputSPs[i].interpolate(cointri.mdata->geometry->inputVertexPositions);
        t[i] = i > 0 ? (t[i - 1] + norm(p[i] - p[i - 1]) / cointri.signpostTri->edgeLengths[intrinsicE]) : 0;
      }

      // insert each intrinsic surface point to sorted list on each input edge
      for (size_t i = 1; i < n - 1; ++i) {
        SurfacePoint intrinsicSP(intrinsicE, t[i]);

        const SurfacePoint& inputSP = inputSPs[i];
        CIT_ASSERT(inputSP.type == SurfacePointType::Edge);

        inputEdgePathOnIntrinsic_map[inputSP.edge][inputSP.tEdge] = intrinsicSP;

        cointri.edgePointMap_intrinsic2input[intrinsicSP] = inputSP;
        cointri.edgePointMap_input2intrinsic[inputSP] = intrinsicSP;

        if (config.topologyValid)
          storeOverlayWedgesAtEdgeIntersection(config, {inputSP, intrinsicSP});
      }
    }

    // store sorted result
    cointri.inputEdgePathOnIntrinsic = EdgeData<std::vector<SurfacePoint>>(cointri.signpostTri->inputMesh);
    for (Edge inputE : cointri.signpostTri->inputMesh.edges()) {
      for (auto p : inputEdgePathOnIntrinsic_map[inputE])
        cointri.inputEdgePathOnIntrinsic[inputE].push_back(p.second);
    }
    DLOG_INFO(logDepth, "  Done");
  };
  computeEdgePath_inner2(config.cointriA);
  computeEdgePath_inner2(config.cointriB);

  // Store maxSignpostError
  ResultInfo& ri = resultInfoTable[config.id];
  ri.A_maxSignpostError = config.cointriA.signpostTri->maxSignpostError(*mdataA.geometry, config.cointriA.intrinsicEdgePath);
  ri.B_maxSignpostError = config.cointriB.signpostTri->maxSignpostError(*mdataB.geometry, config.cointriB.intrinsicEdgePath);

  // Step 3: map each input edge path on intrinsic mesh from P to Q, trace it over Q_inputMesh
  auto computeEdgePath_inner3 = [&config, &logDepth](CoInTri& cointriP, const CoInTri& cointriQ, std::vector<std::array<SurfacePoint, 2>>& inputEdgeMutualIntersections) {
    const ModelData& mdataP = *cointriP.mdata;
    const ModelData& mdataQ = *cointriQ.mdata;

    DLOG_INFO(logDepth, "Tracing input edge path of {} over {} ... ", mdataP.name, mdataQ.name);

    cointriP.inputEdgePathOnOtherInput = EdgeData<std::vector<SurfacePoint>>(cointriP.signpostTri->inputMesh);

    for (Edge P_inputE : cointriP.signpostTri->inputMesh.edges()) {
      Vertex P_inputVStart = P_inputE.halfedge().vertex();
      Vertex P_inputVEnd   = P_inputE.halfedge().tipVertex();

      Vertex P_intrinsicVStart = cointriP.signpostTri->intrinsicMesh->vertex(P_inputVStart.getIndex());
      Vertex P_intrinsicVEnd   = cointriP.signpostTri->intrinsicMesh->vertex(P_inputVEnd  .getIndex());

      Vertex Q_intrinsicVStart = cointriP.correspondingVertex[P_intrinsicVStart];
      Vertex Q_intrinsicVEnd   = cointriP.correspondingVertex[P_intrinsicVEnd];

      SurfacePoint Q_inputSPStart = cointriQ.signpostTri->vertexLocations[Q_intrinsicVStart];
      SurfacePoint Q_inputSPEnd   = cointriQ.signpostTri->vertexLocations[Q_intrinsicVEnd];

      size_t n = cointriP.inputEdgePathOnIntrinsic[P_inputE].size();
      CIT_ASSERT(n >= 2);

      // Push the first point
      cointriP.inputEdgePathOnOtherInput[P_inputE].push_back(Q_inputSPStart);

      SurfacePoint Q_intrinsicSPPrev;     // Needed for getNextTraceVec

      TraceOptions options;
      options.includePath = true;
      TraceGeodesicResult traceResult;    // Result of i-th round may be used by (i+1)-th round

      // We keep track of the parametric position on P_inputE in order to identify overlay vertex of
      // potentially more than one input-input intersections between the same pair of input edges.
      // Note that this parameterization is piecewise arc-length (one piece per (i,i+1) segment)
      double P_tEdge = 0.;

      for (size_t i = 0; i < n - 1; ++i) {
        SurfacePoint P_intrinsicSP     = cointriP.inputEdgePathOnIntrinsic[P_inputE][i];
        SurfacePoint P_intrinsicSPNext = cointriP.inputEdgePathOnIntrinsic[P_inputE][i + 1];

        // These two points correspond to a line segment on P's input face, so this segment's length is simply the Euclidean distance
        double P_segmentLength;
        {
          SurfacePoint P_inputSP     = P_intrinsicSP    .type == SurfacePointType::Vertex ? cointriP.signpostTri->vertexLocations[P_intrinsicSP    .vertex] : cointriP.edgePointMap_intrinsic2input.at(P_intrinsicSP    );
          SurfacePoint P_inputSPNext = P_intrinsicSPNext.type == SurfacePointType::Vertex ? cointriP.signpostTri->vertexLocations[P_intrinsicSPNext.vertex] : cointriP.edgePointMap_intrinsic2input.at(P_intrinsicSPNext);
          CIT_ASSERT(checkAdjacent(P_inputSP, P_inputSPNext));

          Vector3 p0 = P_inputSP    .interpolate(mdataP.geometry->inputVertexPositions);
          Vector3 p1 = P_inputSPNext.interpolate(mdataP.geometry->inputVertexPositions);

          P_segmentLength = norm(p1 - p0);
        }

        bool isSegmentPreserved = false;    // Is this segment (i,i+1) on P_inputE preserved in the intrinsic mesh?

        // This point is an intersection between P_inputE and P's intrinsic edge --> trace geodesic
        if (P_intrinsicSP.type == SurfacePointType::Edge) {
          SurfacePoint Q_intrinsicSP = getCorrespondingEdgePoint(cointriP, P_intrinsicSP);

          // Before doing the next tracing, store facePointB(A) info for OverlayVertex corresponding to this inputA(B)-intrinsic intersection
          if (!cointriQ.signpostTri->isIntrinsicEdgePartiallyOriginal(Q_intrinsicSP.edge)) {
            OverlayVertex overlayVertex = getOverlayVertexAtEdgeIntersection(config, {P_inputE, Q_intrinsicSP.edge});
            for (OverlayWedge& wedge : config.overlayWedgesPerOverlayVertex.at(overlayVertex)) {
              if (mdataP.name == "A")
                wedge.facePointB = traceResult.endPoint;
              else
                wedge.facePointA = traceResult.endPoint;
            }
          }

          // Get the next tracing vector
          SurfacePoint Q_intrinsicSPNext = P_intrinsicSPNext.type == SurfacePointType::Edge
            ? getCorrespondingEdgePoint(cointriP, P_intrinsicSPNext)
            : SurfacePoint(cointriP.correspondingVertex[P_intrinsicSPNext.vertex]);
          double traceVecAngleDelta;
          double traceVecLen;
          detail::getNextTraceVec(cointriQ, {Q_intrinsicSPPrev, Q_intrinsicSP, Q_intrinsicSPNext}, traceVecAngleDelta, traceVecLen);
          double traceVecAngle = traceResult.endingDir.arg() + traceVecAngleDelta;
          Vector2 traceVec = Vector2::fromAngle(traceVecAngle) * traceVecLen;

          // Trace geodesic
          traceResult = traceGeodesic(cointriQ.signpostTri->inputGeom, traceResult.endPoint, traceVec, options);

          // Needed for getNextTraceVec in the next round
          Q_intrinsicSPPrev = Q_intrinsicSP;

          // If P_intrinsicSPNext is a vertex point, trimming may be necessary
          if (P_intrinsicSPNext.type == SurfacePointType::Vertex) {
            Vertex P_intrinsicVNext = P_intrinsicSPNext.vertex;
            SurfacePoint P_inputSPNext = cointriP.signpostTri->vertexLocations[P_intrinsicVNext];

            // Some sanity checks
            if (P_inputSPNext.type == SurfacePointType::Vertex) {
              // The next point is the last end of P_inputE
              CIT_ASSERT(i == n - 2);
              CIT_ASSERT(P_inputSPNext.vertex == P_inputVEnd);
            } else {
              // The next point is an edge point on P_inputE
              CIT_ASSERT(P_inputSPNext.type == SurfacePointType::Edge);
              CIT_ASSERT(P_inputSPNext.edge == P_inputE);
            }

            Vertex Q_intrinsicVNext = cointriP.correspondingVertex[P_intrinsicVNext];
            SurfacePoint Q_inputSPNext = cointriQ.signpostTri->vertexLocations[Q_intrinsicVNext];

            // Trimming is necessary when Q_inputSPNext is not a face point.
            if (Q_inputSPNext.type != SurfacePointType::Face) {
              // Some sanity checks
              if (i < n - 2) {
                // Before the last endpoint, this can happen only when P_intrinsicVNext is an edge point
                CIT_ASSERT(P_inputSPNext.type == SurfacePointType::Edge);
                CIT_ASSERT(Q_inputSPNext.type == SurfacePointType::Vertex);

              } else if (mdataP.anchorVertexIDs.count(P_inputVEnd.getIndex())) {
                // The last endpoint is an anchor vertex
                CIT_ASSERT(mdataQ.anchorVertexIDs.count(Q_inputSPNext.vertex.getIndex()));
                CIT_ASSERT(P_inputSPNext.vertex == P_inputVEnd);
                CIT_ASSERT(Q_inputSPNext.vertex == Q_inputSPEnd.vertex);

              } else {
                // The last endpoint is an edge point on Q
                CIT_ASSERT(Q_inputSPNext.type == SurfacePointType::Edge);
              }

              bool success = trimTraceResult(traceResult, Q_inputSPNext);
              CIT_ASSERT(success);
              if (success)
                traceResult.pathPoints.push_back(Q_inputSPNext);
            }
          }

        // This point is either the first endpoint of P_inputE or an edge point on P_inputE
        } else {
          CIT_ASSERT(P_intrinsicSP.type == SurfacePointType::Vertex);
          Vertex P_intrinsicV = P_intrinsicSP.vertex;
          SurfacePoint P_inputSP = cointriP.signpostTri->vertexLocations[P_intrinsicV];

          // Some sanity checks
          if (P_inputSP.type == SurfacePointType::Vertex) {
            // This point is the first end of P_inputE
            CIT_ASSERT(i == 0);
            CIT_ASSERT(P_inputSP.vertex == P_inputVStart);
          } else {
            // This point is an edge point on P_inputE
            CIT_ASSERT(P_inputSP.type == SurfacePointType::Edge);
            CIT_ASSERT(P_inputSP.edge == P_inputE);
          }

          Vertex Q_intrinsicV = cointriP.correspondingVertex[P_intrinsicV];
          SurfacePoint Q_inputSP = cointriQ.signpostTri->vertexLocations[Q_intrinsicV];

          // The next point is an intrinsic edge point --> trace geodesic
          if (P_intrinsicSPNext.type == SurfacePointType::Edge) {
            SurfacePoint Q_intrinsicSPNext = getCorrespondingEdgePoint(cointriP, P_intrinsicSPNext);

            // Calculate the tracing direction by working on Q's intrinsic face shared by Q_intrinsicV and Q_intrinsicSPNext
            Face Q_intrinsicF = sharedFace(SurfacePoint(Q_intrinsicV), Q_intrinsicSPNext);

            // Get halfedge adjacent to Q_intrinsicF which starts from Q_intrinsicV
            Halfedge Q_intrinsicHE = Q_intrinsicF.halfedge();
            while (Q_intrinsicHE.vertex() != Q_intrinsicV)
              Q_intrinsicHE = Q_intrinsicHE.next();

            double t = Q_intrinsicSPNext.tEdge;
            if (Q_intrinsicSPNext.edge.halfedge() != Q_intrinsicHE.next()) {
              CIT_ASSERT(Q_intrinsicSPNext.edge.halfedge().twin() == Q_intrinsicHE.next());
              t = 1. - t;
            }

            // Obtain a trace vector in the local coordinate system w.r.t. Q_intrinsicHE
            Vector2 traceVec = cointriQ.signpostTri->halfedgeVectorsInFace[Q_intrinsicHE];
            traceVec += t * cointriQ.signpostTri->halfedgeVectorsInFace[Q_intrinsicHE.next()];
            traceVec /= cointriQ.signpostTri->halfedgeVectorsInFace[Q_intrinsicHE].normalize();

            // Transform the vector to a canonical coordinate system w.r.t. Q_intrinsicV
            double angleScaling = 2. * M_PI /  cointriQ.signpostTri->intrinsicVertexAngleSums[Q_intrinsicV];
            double traceVecAngle = angleScaling * (traceVec.arg() + cointriQ.signpostTri->intrinsicHalfedgeDirections[Q_intrinsicHE]);
            double traceVecLen = traceVec.norm();
            traceVec = Vector2::fromAngle(traceVecAngle) * traceVecLen;

            // Trace geodesic
            traceResult = traceGeodesic(cointriQ.signpostTri->inputGeom, Q_inputSP, traceVec, options);

            // Needed for getNextTraceVec in the next round
            Q_intrinsicSPPrev = SurfacePoint(Q_intrinsicV);

          // The next point is also an intrinsic vertex; this segment (i,i+1) on P_inputE is preserved in the intrinsic mesh,
          // so simply copy the intrinsic edge's path
          } else {
            CIT_ASSERT(P_intrinsicSPNext.type == SurfacePointType::Vertex);
            Vertex P_intrinsicVNext = P_intrinsicSPNext.vertex;
            Vertex Q_intrinsicVNext = cointriP.correspondingVertex[P_intrinsicVNext];

            Halfedge P_intrinsicHE = P_intrinsicV.connectingHalfedge(P_intrinsicVNext);
            Halfedge Q_intrinsicHE = Q_intrinsicV.connectingHalfedge(Q_intrinsicVNext);
            CIT_ASSERT(P_intrinsicHE != Halfedge());
            CIT_ASSERT(Q_intrinsicHE != Halfedge());
            CIT_ASSERT(cointriP.correspondingHalfedge[P_intrinsicHE] == Q_intrinsicHE);

            Edge P_intrinsicE = P_intrinsicHE.edge();
            Edge Q_intrinsicE = Q_intrinsicHE.edge();

            traceResult.pathPoints = cointriQ.intrinsicEdgePath[Q_intrinsicE];    // (Repurpose traceResult.pathPoints for the sake of code simplicity)

            // The orientation of P_inputE and P_intrinsicE can be inconsistent even after calling makeIntrinsicEdgesConsistentlyOriented()
            // due to rare cases where P_inputE coincides with an input edge in Q
            if (P_intrinsicE.halfedge() != P_intrinsicHE)
              std::reverse(traceResult.pathPoints.begin(), traceResult.pathPoints.end());

            isSegmentPreserved = true;
          }
        }

        // Compute this segment's length on Q
        double Q_segmentLength = 0.;
        for (size_t j = 0; j < traceResult.pathPoints.size() - 1; ++j) {
          Vector3 p0 = traceResult.pathPoints[j    ].interpolate(mdataQ.geometry->inputVertexPositions);
          Vector3 p1 = traceResult.pathPoints[j + 1].interpolate(mdataQ.geometry->inputVertexPositions);
          Q_segmentLength += norm(p1 - p0);
        }

        for (size_t j = 1; j < traceResult.pathPoints.size(); ++j) {
          cointriP.inputEdgePathOnOtherInput[P_inputE].push_back(traceResult.pathPoints[j]);

          // Get the parametric position on P_inputE corresponding to this path point
          Vector3 p0 = traceResult.pathPoints[j - 1].interpolate(mdataQ.geometry->inputVertexPositions);
          Vector3 p1 = traceResult.pathPoints[j    ].interpolate(mdataQ.geometry->inputVertexPositions);
          P_tEdge += (norm(p1 - p0) / Q_segmentLength) * P_segmentLength / mdataP.geometry->edgeLengths[P_inputE];

          // If both input edges on P & Q are non-preserved, store wedges at the intersection
          if (!isSegmentPreserved && j != traceResult.pathPoints.size() - 1) {
            SurfacePoint P_inputSP = SurfacePoint(P_inputE, P_tEdge);
            SurfacePoint Q_inputSP = traceResult.pathPoints[j];
            CIT_ASSERT(Q_inputSP.type == SurfacePointType::Edge);

            if (!cointriQ.signpostTri->isInputEdgePointPreserved(Q_inputSP))
              inputEdgeMutualIntersections.push_back({P_inputSP, Q_inputSP});    // Queue it for later processing
          }
        }
      }

      // Sanity check for the last endpoint
      SurfacePoint Q_inputSPEnd_check = cointriP.inputEdgePathOnOtherInput[P_inputE].back();
      if (Q_inputSPEnd.type == SurfacePointType::Face) {
        CIT_ASSERT(Q_inputSPEnd.face == Q_inputSPEnd_check.face);
        if (norm(Q_inputSPEnd.faceCoords - Q_inputSPEnd_check.faceCoords) >= 5.e-3) {
          DLOG_ERROR(logDepth, "Sanity check error exceeded threshold: {}", norm(Q_inputSPEnd.faceCoords - Q_inputSPEnd_check.faceCoords));
          CIT_ASSERT(false);
        }

        // Replace it with the one stored in signpost->vertexLocations, for consistency
        cointriP.inputEdgePathOnOtherInput[P_inputE].back() = Q_inputSPEnd;

      } else {
        CIT_ASSERT(Q_inputSPEnd == Q_inputSPEnd_check); // This should hold even for an edge point thanks to trimming
      }
    }
    DLOG_INFO(logDepth, "  Done");
  };
  if (config.topologyValid) {
    std::vector<std::array<SurfacePoint, 2>> inputEdgeMutualIntersections_AtoB;
    std::vector<std::array<SurfacePoint, 2>> inputEdgeMutualIntersections_BtoA;
    computeEdgePath_inner3(config.cointriA, config.cointriB, inputEdgeMutualIntersections_AtoB);
    computeEdgePath_inner3(config.cointriB, config.cointriA, inputEdgeMutualIntersections_BtoA);

    DLOG_INFO(logDepth, "Storing wedges for input-input edge intersections ...");
    // Store wedges for the A->B direction
    for (const auto& i : inputEdgeMutualIntersections_AtoB) {
      storeOverlayWedgesAtEdgeIntersection(config, i);
    }

#if 0
    DLOG_INFO(2, "Checking consistency ...");      // Takes a while
    // Do consistency check for the other direction
    for (const auto& i : inputEdgeMutualIntersections_BtoA) {
      Edge this_eA = i[1].edge;
      Edge this_eB = i[0].edge;
      double this_tEdgeA = i[1].tEdge;
      double this_tEdgeB = i[0].tEdge;

      // Check the existence of a numerically identical overlay vertex in the list
      kt84::MaxMinAverage mma;
      for (const auto& p : config.overlayWedgesPerOverlayVertex) {
        Edge other_eA = p.first.eA;
        Edge other_eB = p.first.eB;
        double other_tEdgeA = p.first.tEdgeA;
        double other_tEdgeB = p.first.tEdgeB;
        if (other_eA == Edge() || other_eA != this_eA) continue;
        if (other_eB == Edge() || other_eB != this_eB) continue;
        mma.update(norm(Vector2{this_tEdgeA - other_tEdgeA, this_tEdgeB - other_tEdgeB}));
      }
      CIT_ASSERT(mma.count() && mma.min() < 1.e-4);
    }
#endif
    DLOG_INFO(logDepth, "  Done");
  }

  config.edgePathComputed = true;
  mdataA.dispList.invalidate();
  mdataB.dispList.invalidate();
}

namespace cit {
namespace detail {

struct HalfedgeFacePair {
  Halfedge halfedge;
  double tHalfedge = -1.;
  Face inputFace;          // Input face to which this halfedge belongs
};
// Sort outgoing input/intrinsic halfedges according to signpost angle
std::map<double, HalfedgeFacePair> getSortedHalfedges(const CoInTri& cointriP, const CoInTri& cointriQ, Vertex P_vIntrinsic) {
  const ModelData& mdataP = *cointriP.mdata;

  std::map<double, HalfedgeFacePair> sortedHalfedges;

  // When registering a portion of a halfedge for an edge point, we need to nudge the tHalfedge value
  // so that the opposite wedges can be robustly separated.
  const double nudgeFactor = 1.e-10;

  Vertex Q_vIntrinsic = cointriP.correspondingVertex[P_vIntrinsic];

  SurfacePoint P_spInput = cointriP.signpostTri->vertexLocations[P_vIntrinsic];
  SurfacePoint Q_spInput = cointriQ.signpostTri->vertexLocations[Q_vIntrinsic];

  // Regardless of P_spInput.type, always register each of intrinsic edges adjacent to P_vIntrinsic
  for (Halfedge P_heIntrinsic : P_vIntrinsic.outgoingHalfedges()) {
    double angle = cointriP.signpostTri->intrinsicHalfedgeDirections[P_heIntrinsic] * (2. * M_PI / cointriP.signpostTri->intrinsicVertexAngleSums[P_vIntrinsic]);

    // If P_heIntrinsic is partially original, set the inputFace field to the associated input face
    Edge P_eInput;
    double tEdgeMin, tEdgeMax;
    bool reversed;
    if (cointriP.signpostTri->isIntrinsicEdgePartiallyOriginal(P_heIntrinsic.edge(), &P_eInput, &tEdgeMin, &tEdgeMax, &reversed)) {
      Face P_fInput;
      if (P_heIntrinsic.getIndex() % 2 == 0) {
        P_fInput = (reversed ? P_eInput.halfedge().twin() : P_eInput.halfedge()).face();
      } else {
        P_fInput = (reversed ? P_eInput.halfedge() : P_eInput.halfedge().twin()).face();
      }
      sortedHalfedges[angle] = {P_heIntrinsic, 0., P_fInput};

    // Otherwise, we'll handle inputFace later
    } else {
      sortedHalfedges[angle] = {P_heIntrinsic, 0., Face()};
    }
  }

  // ==================================================================
  // >>>>> Register Q's non-preserved input halfedges to the list >>>>>
  // ==================================================================
  auto registerQsInputHalfedge = [&](Halfedge Q_heInput, double tHalfedge) {
    // Get the corresponding edge path and find P_spInput in it
    const std::vector<SurfacePoint>& P_edgePath = cointriQ.inputEdgePathOnOtherInput[Q_heInput.edge()];
    auto found = std::find(P_edgePath.begin(), P_edgePath.end(), P_spInput);
    CIT_ASSERT(found != P_edgePath.end());

    SurfacePoint P_spInputNext = *(found + (Q_heInput.getIndex() % 2 == 0 ? 1 : -1));

    // Calculate the vector from P_spInput to P_spInputNext by working on the shared face
    Face P_fInput = sharedFace(P_spInput, P_spInputNext);
    CIT_ASSERT(P_fInput != Face());

    Vector3 faceCoordsDelta = P_spInputNext.inFace(P_fInput).faceCoords - P_spInput.inFace(P_fInput).faceCoords;

    const std::array<Vector2, 3> vertCoords = {{
      {0., 0.},
      mdataP.geometry->halfedgeVectorsInFace[P_fInput.halfedge()],
      -mdataP.geometry->halfedgeVectorsInFace[P_fInput.halfedge().next().next()]
    }};
    Vector2 vectorInface = {0., 0.};
    for (int i = 0; i < 3; ++i)
      vectorInface += faceCoordsDelta[i] * vertCoords[i];

    double angle = vectorInface.arg();

    // Adjust angle as necessary, depending on P_spInput.type
    if (P_spInput.type == SurfacePointType::Face) {
      // No adjustment needed

    } else if (P_spInput.type == SurfacePointType::Edge) {
      Edge P_eInput = P_spInput.edge;

      Halfedge P_heInput = P_eInput.halfedge();
      if (P_heInput.face() != P_fInput)
        P_heInput = P_heInput.twin();
      CIT_ASSERT(P_heInput.face() == P_fInput);

      // Make the angle relative to P_heInput
      angle -= mdataP.geometry->halfedgeVectorsInFace[P_heInput].arg();

      // If P_fInput is not adjacent to P_eInput.halfedge(), add M_PI
      if (P_heInput != P_eInput.halfedge())
        angle += M_PI;

    } else {
      Vertex P_vInput = P_spInput.vertex;

      // Find outgoing halfedge adjacent to P_fInput
      Halfedge P_heInput = P_vInput.halfedge();
      while (P_heInput.face() != P_fInput)
        P_heInput = P_heInput.twin().next();

      // Make the angle relative to P_heInput, followed by scaling & offsetting
      angle -= mdataP.geometry->halfedgeVectorsInFace[P_heInput].arg();
      angle *= 2. * M_PI / cointriP.signpostTri->intrinsicVertexAngleSums[P_vIntrinsic];
      angle += mdataP.geometry->halfedgeVectorsInVertex[P_heInput].arg();
    }

    // Register Q_heInput to the list at appropriate angle
    angle = standardizeAngle(angle);
    sortedHalfedges[angle] = {Q_heInput, tHalfedge, P_fInput};
  };

  // Q_vIntrinsic corresponds to Q's input vertex --> register each of non-preserved input edges adjacent to Q_vInput
  if (Q_spInput.type == SurfacePointType::Vertex) {
    Vertex Q_vInput = Q_spInput.vertex;
    for (Halfedge Q_heInput : Q_vInput.outgoingHalfedges()) {
      // Is Q_heInput preserved?
      double tEdgeTemp = Q_heInput.edge().halfedge() == Q_heInput ? 1.e-15 : 1. - 1.e-15;
      if (!cointriQ.signpostTri->isInputEdgePointPreserved(SurfacePoint(Q_heInput.edge(), tEdgeTemp)))
        registerQsInputHalfedge(Q_heInput, 0.);
    }

  // Q_vIntrinsic corresponds to an edge point on Q --> register each of opposite input halfedges if non-preserved
  } else if (Q_spInput.type == SurfacePointType::Edge) {
    CIT_ASSERT(P_spInput.type == SurfacePointType::Vertex);
    Edge Q_eInput = Q_spInput.edge;

    // Forward direction
    if (!cointriQ.signpostTri->isInputEdgePointPreserved(SurfacePoint(Q_eInput, Q_spInput.tEdge + 1.e-15)))
      registerQsInputHalfedge(Q_eInput.halfedge(), Q_spInput.tEdge + nudgeFactor);

    // Backward direction
    if (!cointriQ.signpostTri->isInputEdgePointPreserved(SurfacePoint(Q_eInput, Q_spInput.tEdge - 1.e-15)))
      registerQsInputHalfedge(Q_eInput.halfedge().twin(), 1. - Q_spInput.tEdge + nudgeFactor);
  }
  // ==================================================================
  // <<<<< Register Q's non-preserved input halfedges to the list <<<<<
  // ==================================================================

  // ==================================================================
  // >>>>> Register P's non-preserved input halfedges to the list >>>>>
  // ==================================================================
  auto registerPsInputHalfedge = [&](Halfedge P_heInput, Halfedge P_heIntrinsic, double angle = -1., double tHalfedge = -1.) {
    Face P_fInput = P_heInput.face();

    // This input halfedge is not preserved
    if (P_heIntrinsic == Halfedge()) {
      CIT_ASSERT(angle != -1. && tHalfedge != -1.);
      sortedHalfedges[angle] = {P_heInput, tHalfedge, P_fInput};

    // This input halfedge is preserved and thus has already been added to the list above
    // Find the corresponding entry in the list and set inputFace
    } else {
      CIT_ASSERT(angle == -1. && tHalfedge == -1.);

      bool found = false;
      for (auto& p : sortedHalfedges) {
        // We can't just use Halfedge::operator==() because it doesn't look at the referenced mesh
        if (strictlyEqual(p.second.halfedge, P_heIntrinsic)) {
          p.second.inputFace = P_fInput;
          found = true;
          break;
        }
      }
      CIT_ASSERT(found);
    }
  };

  // P_vIntrinsic corresponds to P's input vertex --> register each of non-preserved input edges adjacent to P_vInput
  if (P_spInput.type == SurfacePointType::Vertex) {
    Vertex P_vInput = P_spInput.vertex;
    for (Halfedge P_heInput : P_vInput.outgoingHalfedges()) {
      // Is P_heInput preserved?
      double tEdgeTemp = P_heInput.edge().halfedge() == P_heInput ? 1.e-15 : 1. - 1.e-15;
      SurfacePoint P_spIntrinsic;
      // If preserved, set the corresponding entry's inputFace
      if (cointriP.signpostTri->isInputEdgePointPreserved(SurfacePoint(P_heInput.edge(), tEdgeTemp), &P_spIntrinsic)) {
        CIT_ASSERT(P_spIntrinsic.type == SurfacePointType::Edge);
        Halfedge P_heIntrinsic = P_spIntrinsic.edge.halfedge();

        if (P_heIntrinsic.vertex() != P_vIntrinsic)
          P_heIntrinsic = P_heIntrinsic.twin();
        CIT_ASSERT(P_heIntrinsic.vertex() == P_vIntrinsic);

        registerPsInputHalfedge(P_heInput, P_heIntrinsic);

      // Otherwise, register P_heInput to the list
      } else {
        double angle = mdataP.geometry->halfedgeVectorsInVertex[P_heInput].arg();
        angle = standardizeAngle(angle);
        registerPsInputHalfedge(P_heInput, Halfedge(), angle, 0.);
      }
    }

  // P_vIntrinsic corresponds to an edge point on P --> register each of opposite input halfedges if non-preserved
  } else if (P_spInput.type == SurfacePointType::Edge) {
    CIT_ASSERT(Q_spInput.type == SurfacePointType::Vertex);
    Edge     P_eInput  = P_spInput.edge;
    Halfedge P_heInput = P_eInput.halfedge();

    // Forward direction
    SurfacePoint P_spIntrinsic;
    if (cointriP.signpostTri->isInputEdgePointPreserved(SurfacePoint(P_eInput, P_spInput.tEdge + 1.e-15), &P_spIntrinsic)) {
      CIT_ASSERT(P_spIntrinsic.type == SurfacePointType::Edge);
      Halfedge P_heIntrinsic = P_spIntrinsic.edge.halfedge();

      if (P_heIntrinsic.vertex() != P_vIntrinsic)
        P_heIntrinsic = P_heIntrinsic.twin();
      CIT_ASSERT(P_heIntrinsic.vertex() == P_vIntrinsic);

      registerPsInputHalfedge(P_heInput, P_heIntrinsic);

    } else {
      registerPsInputHalfedge(P_heInput, Halfedge(), 0., P_spInput.tEdge + nudgeFactor);
    }

    // Backward direction
    if (cointriP.signpostTri->isInputEdgePointPreserved(SurfacePoint(P_eInput, P_spInput.tEdge - 1.e-15), &P_spIntrinsic)) {
      CIT_ASSERT(P_spIntrinsic.type == SurfacePointType::Edge);
      Halfedge P_heIntrinsic = P_spIntrinsic.edge.halfedge();

      if (P_heIntrinsic.vertex() != P_vIntrinsic)
        P_heIntrinsic = P_heIntrinsic.twin();
      CIT_ASSERT(P_heIntrinsic.vertex() == P_vIntrinsic);

      registerPsInputHalfedge(P_heInput.twin(), P_heIntrinsic);

    } else {
      registerPsInputHalfedge(P_heInput.twin(), Halfedge(), M_PI, 1. - P_spInput.tEdge + nudgeFactor);
    }
  }
  // ==================================================================
  // <<<<< Register P's non-preserved input halfedges to the list <<<<<
  // ==================================================================

  // If P_vIntrinsic lies on P's input face, all the halfedges should be associated with that face
  if (P_spInput.type == SurfacePointType::Face) {
    for (auto& p : sortedHalfedges)
      p.second.inputFace = P_spInput.face;

  // Otherwise, fill in missing inputFace field for entries corresponding to non-original intrinsic edges
  } else {
    for (auto i = sortedHalfedges.begin(); i != sortedHalfedges.end(); ++i) {
      if (i->second.inputFace == Face()) {
        // Walk clockwise to find entry with inputFace filled
        auto j = i;
        do {
          if (j == sortedHalfedges.begin())
            j = sortedHalfedges.end();
          --j;
        } while (j->second.inputFace == Face());

        i->second.inputFace = j->second.inputFace;
      }
    }
  }

  return sortedHalfedges;
}

struct HalfedgeFaceTriplet {
  Halfedge halfedge;
  double tHalfedge = -1.;
  Face A_inputFace;        // Input face in A to which this halfedge belongs
  Face B_inputFace;        // Input face in B to which this halfedge belongs
};
std::map<double, HalfedgeFaceTriplet> mergeSortedHalfedges(
  const Configuration& config,
  const std::map<double, HalfedgeFacePair>& A_sortedHalfedges,
  const std::map<double, HalfedgeFacePair>& B_sortedHalfedges
) {
  std::map<double, HalfedgeFaceTriplet> AB_sortedHalfedges;

  CIT_ASSERT(A_sortedHalfedges.size() == B_sortedHalfedges.size());

  auto iA = A_sortedHalfedges.begin();
  auto iB = B_sortedHalfedges.begin();

  // Align iB to iA so that they point to the same halfedge
  for (; iB != B_sortedHalfedges.end(); ++iB) {
    Halfedge A_halfedge = iA->second.halfedge;
    Halfedge B_halfedge = iB->second.halfedge;

    // We can't just use Halfedge::operator==() because it doesn't look at the referenced mesh
    if (strictlyEqual(A_halfedge, B_halfedge))
      break;
  }
  CIT_ASSERT(iB != B_sortedHalfedges.end());

  for (; iA != A_sortedHalfedges.end(); ++iA, ++iB) {
    // iB wraps around
    if (iB == B_sortedHalfedges.end())
      iB = B_sortedHalfedges.begin();

    double angle = iA->first;

    HalfedgeFacePair A_hfp = iA->second;
    HalfedgeFacePair B_hfp = iB->second;

    // The two lists should have the same set of halfedges in the same order
    CIT_ASSERT(A_hfp.halfedge == B_hfp.halfedge);
    CIT_ASSERT(A_hfp.tHalfedge == B_hfp.tHalfedge);

    // The inputFace field should be non-null for both
    CIT_ASSERT(A_hfp.inputFace.getMesh() == mdataA.mesh.get());
    CIT_ASSERT(B_hfp.inputFace.getMesh() == mdataB.mesh.get());

    AB_sortedHalfedges[angle] = {A_hfp.halfedge, A_hfp.tHalfedge, A_hfp.inputFace, B_hfp.inputFace};
  }

  return AB_sortedHalfedges;
}

struct OverlayWedgesPerHalfedge {
  HalfedgeData<std::map<double, OverlayWedge>> A, B, intrinsic;
};

OverlayWedge getNextWedge(const OverlayWedgesPerHalfedge& overlayWedgesPerHalfedge, const OverlayWedge& wedge) {
  SurfaceMesh* mesh = wedge.halfedgeOut.getMesh();

  const std::map<double, OverlayWedge>& wedges = mesh == mdataA.mesh.get()
    ? overlayWedgesPerHalfedge.A[wedge.halfedgeOut]
    : mesh == mdataB.mesh.get()
      ? overlayWedgesPerHalfedge.B[wedge.halfedgeOut]
      : overlayWedgesPerHalfedge.intrinsic[wedge.halfedgeOut];

  for (const auto& p : wedges) {
    CIT_ASSERT(p.first == p.second.tHalfedgeIn);
    if (p.first > wedge.tHalfedgeOut)
      return p.second;
  }
  throw std::logic_error("Couldn't find the next overlay wedge");
}

OverlayWedge getOppositeWedge(const OverlayWedgesPerHalfedge& overlayWedgesPerHalfedge, const OverlayWedge& wedge) {
  SurfaceMesh* mesh = wedge.halfedgeOut.getMesh();

  const std::map<double, OverlayWedge>& wedges = mesh == mdataA.mesh.get()
    ? overlayWedgesPerHalfedge.A[wedge.halfedgeOut.twin()]
    : mesh == mdataB.mesh.get()
      ? overlayWedgesPerHalfedge.B[wedge.halfedgeOut.twin()]
      : overlayWedgesPerHalfedge.intrinsic[wedge.halfedgeOut.twin()];

  for (const auto& p : wedges) {
    if (std::fabs(p.first - (1. - wedge.tHalfedgeOut)) < 1.e-10)
      return p.second;
  }
  throw std::logic_error("Couldn't find the opposite overlay wedge");
}

} // namespace detail
} // namespace cit

void cit::computeOverlayPolygons(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeOverlayPolygons");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeOverlayPolygons"); });

  if (!config.topologyValid)
    return;

  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  config.overlayPolygons = {};

  // Step 1: For each intrinsic vertex, generate an overlay mesh vertex and its adjacent wedges
  DLOG_DEBUG(logDepth + 1, "Step 1");
  for (Vertex A_vIntrinsic : cointriA.signpostTri->intrinsicMesh->vertices()) {
    Vertex B_vIntrinsic = cointriA.correspondingVertex[A_vIntrinsic];

    SurfacePoint A_spInput = cointriA.signpostTri->vertexLocations[A_vIntrinsic];
    SurfacePoint B_spInput = cointriB.signpostTri->vertexLocations[B_vIntrinsic];

    // Generate an overlay vertex in the relevant category
    OverlayVertex overlayVertex;
    switch(A_spInput.type) {
      case SurfacePointType::Vertex : overlayVertex.vA = A_spInput.vertex; break;
      case SurfacePointType::Edge   : overlayVertex.eA = A_spInput.edge  ; break;
      case SurfacePointType::Face   : break;
    }
    switch(B_spInput.type) {
      case SurfacePointType::Vertex : overlayVertex.vB = B_spInput.vertex; break;
      case SurfacePointType::Edge   : overlayVertex.eB = B_spInput.edge  ; break;
      case SurfacePointType::Face   : break;
    }

    // Assertion for anchor vertex
    if (overlayVertex.vA != Vertex() && overlayVertex.vB != Vertex()) {
      CIT_ASSERT(mdataA.anchorVertexIDs.count(A_spInput.vertex.getIndex()));
      CIT_ASSERT(mdataB.anchorVertexIDs.count(B_spInput.vertex.getIndex()));
    }

    // Generate a sorted list of outgoing halfedges for each of A & B
    std::map<double, detail::HalfedgeFacePair> A_sortedHalfedges = detail::getSortedHalfedges(cointriA, cointriB, A_vIntrinsic);
    std::map<double, detail::HalfedgeFacePair> B_sortedHalfedges = detail::getSortedHalfedges(cointriB, cointriA, B_vIntrinsic);

    // Convert intrinsic halfedge on B to its corresponding one on A, following convention (see comment in OverlayVertex definition)
    for (auto& p : B_sortedHalfedges) {
      auto& B_hfp = p.second;
      if (B_hfp.halfedge.getMesh() == config.cointriB.signpostTri->intrinsicMesh.get())
        B_hfp.halfedge = config.cointriB.correspondingHalfedge[B_hfp.halfedge];
    }

    // Merge the two lists
    std::map<double, detail::HalfedgeFaceTriplet> AB_sortedHalfedges = detail::mergeSortedHalfedges(config, A_sortedHalfedges, B_sortedHalfedges);

    // Make wedges
    for (auto j0 = AB_sortedHalfedges.begin(); j0 != AB_sortedHalfedges.end(); ++j0) {
      auto j1 = j0;
      ++j1;
      if (j1 == AB_sortedHalfedges.end())
        j1 = AB_sortedHalfedges.begin();

      Halfedge he0 = j0->second.halfedge;
      Halfedge he1 = j1->second.halfedge;

      Face A_fInput = j0->second.A_inputFace;
      Face B_fInput = j0->second.B_inputFace;

      OverlayWedge wedge;
      wedge.overlayVertex = overlayVertex;
      wedge.halfedgeIn = he1.twin();
      wedge.halfedgeOut = he0;
      wedge.tHalfedgeIn = 1. - j1->second.tHalfedge;
      wedge.tHalfedgeOut = j0->second.tHalfedge;
      wedge.facePointA = A_spInput.inFace(A_fInput);
      wedge.facePointB = B_spInput.inFace(B_fInput);

      config.overlayWedgesPerOverlayVertex[overlayVertex].push_back(wedge);
    }
  }

  // Step 2: Group wedges by halfedgeIn, sort each group by tHalfedgeIn
  DLOG_DEBUG(logDepth + 1, "Step 2");
  detail::OverlayWedgesPerHalfedge overlayWedgesPerHalfedge = {
    HalfedgeData<std::map<double, OverlayWedge>>(*mdataA.mesh),
    HalfedgeData<std::map<double, OverlayWedge>>(*mdataB.mesh),
    HalfedgeData<std::map<double, OverlayWedge>>(*cointriA.signpostTri->intrinsicMesh)
  };
  for (const auto& p : config.overlayWedgesPerOverlayVertex) {
    for (const OverlayWedge& wedge : p.second) {
      SurfaceMesh* mesh = wedge.halfedgeIn.getMesh();
      if (mesh == mdataA.mesh.get())
        overlayWedgesPerHalfedge.A[wedge.halfedgeIn][wedge.tHalfedgeIn] = wedge;
      else if (mesh == mdataB.mesh.get())
        overlayWedgesPerHalfedge.B[wedge.halfedgeIn][wedge.tHalfedgeIn] = wedge;
      else
        overlayWedgesPerHalfedge.intrinsic[wedge.halfedgeIn][wedge.tHalfedgeIn] = wedge;
    }
  }

  // Step 3: Generate polygons by finding next wedge using the sorted result
  DLOG_DEBUG(logDepth + 1, "Step 3");
  std::unordered_set<OverlayWedge> visited;
  for (const auto& p : config.overlayWedgesPerOverlayVertex) {
    for (const OverlayWedge& wedgeStart : p.second) {
      if (visited.count(wedgeStart))
        continue;

      config.overlayPolygons.push_back({});
      OverlayPolygon& polygon = config.overlayPolygons.back();

      OverlayWedge wedgeCurr = wedgeStart;
      std::unordered_set<OverlayVertex> duplicateCheck;
      do {
        CIT_ASSERT(polygon.wedges.size() < 100);

        // Check duplicates
        CIT_ASSERT(!duplicateCheck.count(wedgeCurr.overlayVertex));
        duplicateCheck.insert(wedgeCurr.overlayVertex);

        wedgeCurr.overlayPolygon_ptr = &polygon;
        polygon.wedges.push_back(wedgeCurr);

        visited.insert(wedgeCurr);

        wedgeCurr = detail::getNextWedge(overlayWedgesPerHalfedge, wedgeCurr);

      } while (wedgeCurr != wedgeStart);

    }
  }

  // Step 4: Set fIntrinsic for polygons adjacent to intrinsic edges
  DLOG_DEBUG(logDepth + 1, "Step 4");
  std::set<OverlayPolygon*> unassignedPolygons;
  for (OverlayPolygon& polygon : config.overlayPolygons) {
    for (const OverlayWedge& wedge : polygon.wedges) {
      if (wedge.halfedgeOut.getMesh() == cointriA.signpostTri->intrinsicMesh.get()) {
        polygon.fIntrinsic = wedge.halfedgeOut.face();
        break;
      }
    }
    if (!polygon.fIntrinsic.getMesh())
      unassignedPolygons.insert(&polygon);
  }

  // Step 5: Set fIntrinsic for remaining polygons (those surrounded by input edges only) via floodfill
  DLOG_DEBUG(logDepth + 1, "Step 5");
  // Assign overlayPolygon_ptr for each wedge in overlayWedgesPerHalfedge.*
  for (OverlayPolygon& polygon : config.overlayPolygons) {
    for (const OverlayWedge& wedge : polygon.wedges) {
      SurfaceMesh* mesh = wedge.halfedgeIn.getMesh();

      std::map<double, OverlayWedge>& wedges = mesh == mdataA.mesh.get()
        ? overlayWedgesPerHalfedge.A[wedge.halfedgeIn]
        : mesh == mdataB.mesh.get()
          ? overlayWedgesPerHalfedge.B[wedge.halfedgeIn]
          : overlayWedgesPerHalfedge.intrinsic[wedge.halfedgeIn];

      wedges.at(wedge.tHalfedgeIn).overlayPolygon_ptr = &polygon;
    }
  }

  while (!unassignedPolygons.empty()) {
    for (auto i = unassignedPolygons.begin(); i != unassignedPolygons.end(); ) {
      OverlayPolygon& polygon = *(*i);
      bool found = false;
      for (const OverlayWedge& wedge : polygon.wedges) {
        OverlayWedge wedgeOpp = detail::getOppositeWedge(overlayWedgesPerHalfedge, wedge);
        if (wedgeOpp.overlayPolygon_ptr->fIntrinsic.getMesh()) {
          polygon.fIntrinsic = wedgeOpp.overlayPolygon_ptr->fIntrinsic;
          found = true;
          break;
        }
      }
      if (found)
        i = unassignedPolygons.erase(i);
      else
        ++i;
    }
  }
}

double cit::computeEnergy(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN computeEnergy");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   computeEnergy"); });

  if (!config.topologyValid)
    return (config.energy = std::numeric_limits<double>::infinity());

  config.energy = 0.;
  config.energyDensity = FaceData<double>(*config.cointriA.signpostTri->intrinsicMesh);
  kt84::MaxMinAverage mma;
  for (Face A_f : config.cointriA.signpostTri->intrinsicMesh->faces()) {
    Halfedge A_he = A_f.halfedge();
    Halfedge B_he = config.cointriA.correspondingHalfedge[A_he];
    Face B_f = B_he.face();

    std::array<Vector2d, 3> A_triangle;
    std::array<Vector2d, 3> B_triangle;
    A_triangle[0] = {0,0};
    B_triangle[0] = {0,0};
    A_triangle[1] = config.cointriA.signpostTri->halfedgeVectorsInFace[A_he].castTo<Vector2d>();
    B_triangle[1] = config.cointriB.signpostTri->halfedgeVectorsInFace[B_he].castTo<Vector2d>();
    A_triangle[2] = -config.cointriA.signpostTri->halfedgeVectorsInFace[A_he.next().next()].castTo<Vector2d>();
    B_triangle[2] = -config.cointriB.signpostTri->halfedgeVectorsInFace[B_he.next().next()].castTo<Vector2d>();

    double A_area;
    double B_area;
    config.energyDensity[A_f] = computeTriangleEnergy(A_triangle, B_triangle, &A_area, &B_area);
    config.energy += config.energyDensity[A_f];
    config.energyDensity[A_f] /= 0.5 * (A_area + B_area);
    mma.update(config.energyDensity[A_f]);
  }
  LOG_INFO(logger, logDepth, "  Energy value: {}", config.energy);
  LOG_INFO(logger, logDepth, "  Energy density stats: {}", mma);
  return config.energy;
}

void cit::computeDerivative(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeDerivative");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeDerivative"); });

  if (!config.topologyValid)
    return;

  if (!config.edgePathComputed)
    computeEdgePath(config, logDepth + 1);

  const int n = 2 * (mdataA.nV + mdataB.nV);
  double energy = 0.;
  config.gradient = VectorXd::Zero(n);
  std::vector<Eigen::Triplet<double>> triplets;

  for (Face A_f : config.cointriA.signpostTri->intrinsicMesh->faces()) {
    std::array<Vertex, 3> A_v;
    std::array<Vertex, 3> B_v;
    int k = 0;
    for (Vertex A_v_ : A_f.adjacentVertices()) {
      A_v[k] = A_v_;
      B_v[k] = config.cointriA.correspondingVertex[A_v[k]];
      ++k;
    }

    // compute per-face energy
    ADScalar_6v_2nd energy_f = computeTriangleEnergyAD(config, A_f);
    energy += energy_f.value().value();

    // copy gradient
    Vector6d gradient_f = energy_f.value().derivatives();
    std::set<int> ignoredRows;
    for (int i = 0; i < 3; ++i) {
      int idx;
      if (A_v[i].getIndex() < mdataA.nV && B_v[i].getIndex() < mdataB.nV) {
        // Both vertices on A & B are fixed
        CIT_ASSERT(mdataA.anchorVertexIDs.count(A_v[i].getIndex()));
        CIT_ASSERT(mdataB.anchorVertexIDs.count(B_v[i].getIndex()));
        idx = -1;
        ignoredRows.insert(i);

      } else if (A_v[i].getIndex() < mdataA.nV) {
        CIT_ASSERT(B_v[i].getIndex() >= mdataB.nV);
        // this vertex originates in A, thus fixed on modelA and mutable on modelB
        idx = A_v[i].getIndex();

      } else {
        // this vertex originates in B, thus fixed on modelB and mutable on modelA
        CIT_ASSERT(B_v[i].getIndex() < mdataB.nV);
        idx = mdataA.nV + B_v[i].getIndex();
      }

      for (int di = 0; di < 2; ++di) {
        if (idx == -1) {
          CIT_ASSERT(gradient_f[2*i + di] == 0.);
          continue;
        }
        config.gradient[2*idx + di] += gradient_f[2*i + di];
      }
    }

    // project hessian to positive definite
    Matrix6d hessian_f = getHessian(energy_f);

    if (ignoredRows.empty()) {
      hessian_f = projectHessian(hessian_f);

    } else if (ignoredRows.size() == 1) {
      int ignoredRow = *ignoredRows.begin();
      Matrix4d hessian_sub;
      if (ignoredRow == 0 || ignoredRow == 2) {
        hessian_sub = hessian_f.block<4, 4>(2 * ((ignoredRow + 1) % 3), 2 * ((ignoredRow + 1) % 3));
      } else {
        hessian_sub.block<2, 2>(0, 0) = hessian_f.block<2, 2>(0, 0);
        hessian_sub.block<2, 2>(2, 0) = hessian_f.block<2, 2>(4, 0);
        hessian_sub.block<2, 2>(0, 2) = hessian_f.block<2, 2>(0, 4);
        hessian_sub.block<2, 2>(2, 2) = hessian_f.block<2, 2>(4, 4);
      }
      hessian_sub = projectHessian(hessian_sub);
      if (ignoredRow == 0 || ignoredRow == 2) {
        hessian_f.block<4, 4>(2 * ((ignoredRow + 1) % 3), 2 * ((ignoredRow + 1) % 3)) = hessian_sub;
      } else {
        hessian_f.block<2, 2>(0, 0) = hessian_sub.block<2, 2>(0, 0);
        hessian_f.block<2, 2>(4, 0) = hessian_sub.block<2, 2>(2, 0);
        hessian_f.block<2, 2>(0, 4) = hessian_sub.block<2, 2>(0, 2);
        hessian_f.block<2, 2>(4, 4) = hessian_sub.block<2, 2>(2, 2);
      }

    } else if (ignoredRows.size() == 2) {
      int takenRow = !ignoredRows.count(0) ? 0 : !ignoredRows.count(1) ? 1 : 2;
      Matrix2d hessian_sub = hessian_f.block<2, 2>(2 * takenRow, 2 * takenRow);
      hessian_sub = projectHessian(hessian_sub);
      hessian_f.block<2, 2>(2 * takenRow, 2 * takenRow) = hessian_sub;

    } else {
      // All three vertices are anchors, no variable is relevant for this face's energy
      continue;
    }

    // copy hessian
    for (int i = 0; i < 3; ++i) {
      int idx_i;
      if (A_v[i].getIndex() < mdataA.nV && B_v[i].getIndex() < mdataB.nV) {
        idx_i = -1;
      } else if (A_v[i].getIndex() < mdataA.nV) {
        idx_i = A_v[i].getIndex();
      } else {
        idx_i = mdataA.nV + B_v[i].getIndex();
      }
      for (int j = 0; j < 3; ++j) {
        int idx_j;
        if (A_v[j].getIndex() < mdataA.nV && B_v[j].getIndex() < mdataB.nV) {
          idx_j = -1;
        } else if (A_v[j].getIndex() < mdataA.nV) {
          idx_j = A_v[j].getIndex();
        } else {
          idx_j = mdataA.nV + B_v[j].getIndex();
        }
        for (int di = 0; di < 2; ++di) {
          for (int dj = 0; dj < 2; ++dj) {
            int idx_ii = 2*idx_i + di;
            int idx_jj = 2*idx_j + dj;
            int ii = 2*i + di;
            int jj = 2*j + dj;
            if (idx_i == -1 || idx_j == -1) {
              CIT_ASSERT(hessian_f(ii, jj) == 0.);
              continue;
            }
            triplets.emplace_back(idx_ii, idx_jj, hessian_f(ii, jj));
          }
        }
      }
    }
  }

  config.hessian = SparseMatrixd{n, n};
  config.hessian.setFromTriplets(triplets.begin(), triplets.end());
  config.hessian.makeCompressed();

  DLOG_INFO(logDepth + 1, "energy: {} (sanity check: {}), gradient norm: {}", energy, std::fabs(energy - config.energy), config.gradient.norm());

  config.derivativeComputed = true;
}

void cit::computeLaplacianTransportOperator(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeLaplacianTransportOperator");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeLaplacianTransportOperator"); });

  if (!config.edgePathComputed)
    computeEdgePath(config, logDepth + 1);

  auto computeLaplacianTransportOperator_inner = [](CoInTri& cointriP, const CoInTri& cointriQ) {
    cointriP.laplacianTransportOperator = EdgeData<Matrix2d>(cointriP.signpostTri->inputMesh);
    cointriP.laplacianTransportOperator.fill(Matrix2d::Identity());
    for (Edge P_e : cointriP.signpostTri->inputMesh.edges()) {
      // Skip edges adjacent to anchor vertex
      if (cointriP.mdata->anchorVertexIDs.count(P_e.halfedge().vertex().getIndex())) continue;
      if (cointriP.mdata->anchorVertexIDs.count(P_e.halfedge().twin().vertex().getIndex())) continue;

      Matrix2d& t = cointriP.laplacianTransportOperator[P_e];
      const std::vector<SurfacePoint>& Q_edgePath = cointriP.inputEdgePathOnOtherInput[P_e];
      const size_t n = Q_edgePath.size();

      // If the first end is an edge point, and its halfedge().face() is not adjacent to the next point in the path, apply the transport operator
      SurfacePoint Q_spStart0 = Q_edgePath[0];
      if (Q_spStart0.type == SurfacePointType::Edge) {
        SurfacePoint Q_spStart1 = Q_edgePath[1];
        Face Q_f = sharedFace(Q_spStart0, Q_spStart1);
        CIT_ASSERT(Q_f != Face());

        if (Q_spStart0.edge.halfedge().face() != Q_f)
          t = cointriQ.mdata->transportOperator[Q_spStart0.edge];
      }

      // If any of the edge path is a vertex point, we'll take two paths (left side & right side) and average the result
      bool needsTwice = false;
      Matrix2d result      = Matrix2d::Identity();
      Matrix2d resultRight = Matrix2d::Identity();

      for (int rightSide = 0; rightSide < 2; ++rightSide) {
        for (size_t k = 1; k < n - 1; ++k) {
          SurfacePoint Q_spCurr = Q_edgePath[k];
          SurfacePoint Q_spPrev = Q_edgePath[k - 1];
          SurfacePoint Q_spNext = Q_edgePath[k + 1];

          if (Q_spCurr.type == SurfacePointType::Edge) {
            Edge Q_eCurr = Q_spCurr.edge;

            Face Q_f = sharedFace(Q_spCurr, Q_spPrev);
            CIT_ASSERT(Q_f != Face());

            bool reversed = Q_eCurr.halfedge().twin().face() == Q_f;

            Matrix2d tNext = cointriQ.mdata->transportOperator[Q_eCurr];
            if (reversed)
              tNext = Matrix2d(tNext).inverse();  // mat = mat.inverse() won't work

            if (rightSide)
              resultRight = tNext * resultRight;
            else
              result = tNext * result;

          } else if (Q_spCurr.type == SurfacePointType::Vertex) {
            needsTwice = true;
            Vertex Q_vCurr = Q_spCurr.vertex;

            // We walk over outgoing halfedges around Q_vCurr clockwise (left side) or counterclockwise (right side)
            // Determine the first halfedge of the walk
            Halfedge Q_heStart;
            if (Q_spPrev.type == SurfacePointType::Vertex) {
              Vertex Q_vPrev = Q_spPrev.vertex;

              Q_heStart = rightSide
                ? Q_vCurr.connectingHalfedge(Q_vPrev).next().next().twin()
                : Q_vPrev.connectingHalfedge(Q_vCurr).next();

            } else {
              Face Q_f = sharedFace(Q_spPrev, Q_spCurr);
              CIT_ASSERT(Q_f != Face());

              Q_heStart = Q_f.oppositeHalfedge(Q_vCurr);

              Q_heStart = rightSide
                ? Q_heStart.next().twin()
                : Q_heStart.next().next();
            }

            // Determine the last halfedge of the walk
            Halfedge Q_heEnd;
            if (Q_spNext.type == SurfacePointType::Vertex) {
              Vertex Q_vNext = Q_spNext.vertex;

              Q_heEnd = rightSide
                ? Q_vNext.connectingHalfedge(Q_vCurr).next()
                : Q_vCurr.connectingHalfedge(Q_vNext).next().next().twin();

            } else {
              Face Q_f = sharedFace(Q_spCurr, Q_spNext);
              CIT_ASSERT(Q_f != Face());

              Q_heEnd = Q_f.oppositeHalfedge(Q_vCurr);

              Q_heEnd = rightSide
                ? Q_heEnd.next().next()
                : Q_heEnd.next().twin();
            }

            // Do the walk
            Halfedge Q_he = Q_heStart;
            while (true) {
              Edge Q_e = Q_he.edge();

              bool reversed = rightSide
                ? Q_e.halfedge() == Q_he
                : Q_e.halfedge().twin() == Q_he;

              Matrix2d tNext = cointriQ.mdata->transportOperator[Q_e];
              if (reversed)
                tNext = Matrix2d(tNext).inverse();  // mat = mat.inverse() won't work

              if (rightSide)
                resultRight = tNext * resultRight;
              else
                result = tNext * result;

              if (Q_he == Q_heEnd)
                break;

              Q_he = rightSide
                ? Q_he.next().next().twin()   // Counterclockwise
                : Q_he.twin().next();         // Clockwise
            }

          } else {
            continue;
          }
        }

        if (!needsTwice)
          break;
      }

      if (needsTwice)
        result = 0.5 * (result + resultRight);

      t = result * t;

      // If the last end is an edge point, and its halfedge().face() is not adjacent to the previous point in the path, apply the transport operator
      SurfacePoint Q_spEnd0 = Q_edgePath[n - 1];
      if (Q_spEnd0.type == SurfacePointType::Edge) {
        SurfacePoint Q_spEnd1 = Q_edgePath[n - 2];
        Face Q_f = sharedFace(Q_spEnd0, Q_spEnd1);
        CIT_ASSERT(Q_f != Face());

        if (Q_spEnd0.edge.halfedge().face() != Q_f)
          t = cointriQ.mdata->transportOperator[Q_spEnd0.edge].inverse() * t;     // Note the need for inversion
      }
    }
  };
  computeLaplacianTransportOperator_inner(config.cointriA, config.cointriB);
  computeLaplacianTransportOperator_inner(config.cointriB, config.cointriA);

  config.laplacianTransportOperatorComputed = true;
}

void cit::computeTemporalTransportOperator(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeTemporalTransportOperator");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeTemporalTransportOperator"); });

  if (!config.topologyValid)
    return;

  auto computeTemporalTransportOperator_inner = [&logDepth](CoInTri& cointriP, const CoInTri& cointriQ) {
    const ModelData& mdataQ = *cointriQ.mdata;

    cointriP.temporalTransportOperator = VertexData<Matrix2d>(cointriP.signpostTri->inputMesh);
    cointriP.temporalTransportOperator.fill(Matrix2d::Identity());

    for (Vertex P_inputV : cointriP.signpostTri->inputMesh.vertices()) {
      // Skip anchor vertex
      if (cointriP.mdata->anchorVertexIDs.count(P_inputV.getIndex())) continue;

      CIT_ASSERT(!cointriP.vertexPath[P_inputV].empty());

      // If vertexPath[P_inputV].back() and vertexLocations[Q_intrinsicV] are different due to modifications
      // to the vertex positions (via merge/split, relocation), append additional faces in a strip connecting them
      Vertex P_intrinsicV = cointriP.signpostTri->intrinsicMesh->vertex(P_inputV.getIndex());
      Vertex Q_intrinsicV = cointriP.correspondingVertex[P_intrinsicV];
      SurfacePoint Q_spInput = cointriP.vertexPath[P_inputV].back();                  // resulting from displaceSolutionVector(), input to createConfiguration()
      SurfacePoint Q_spOutput = cointriQ.signpostTri->vertexLocations[Q_intrinsicV];  // resulting from createConfiguration()

      // If the output is an edge point, it must be the same as the input
      if (Q_spOutput.type == SurfacePointType::Edge) {
        // tEdge changes slightly due to equivalentPointOnIntrinsic() & isIntrinsicEdgePartiallyOriginal() called in insertVertices()
        CIT_ASSERT(Q_spInput.edge == Q_spOutput.edge);
        DLOG_DEBUG(logDepth + 1, "Edge point consistency check: {}", std::fabs(Q_spInput.tEdge - Q_spOutput.tEdge));
        CIT_ASSERT(std::fabs(Q_spInput.tEdge - Q_spOutput.tEdge) < 5.e-4);

      // The usual case: the output is a face point
      } else {
        CIT_ASSERT(Q_spOutput.type == SurfacePointType::Face);

        if (Q_spInput.type == SurfacePointType::Face) {
          if (Q_spInput.face != Q_spOutput.face) {
            // Find face strip between the two points
            DLOG_DEBUG(logDepth + 1, "Appending additional face strip to vertexPath (vertex ID: {})", cointriP.uniqueID_per_Vertex[P_intrinsicV]);
            std::vector<Halfedge> Q_inputFaceStrip = findFaceStrip(Q_spInput.face, Q_spOutput.face);

            // Add edge midpoints of the face strip to the vertex path, followed by the actual vertex location
            for (Halfedge Q_inputHE : Q_inputFaceStrip)
              cointriP.vertexPath[P_inputV].push_back(SurfacePoint(Q_inputHE.edge(), 0.5));
            cointriP.vertexPath[P_inputV].push_back(Q_spOutput);
          }

        } else {
          CIT_ASSERT(Q_spInput.type == SurfacePointType::Edge);

          CIT_ASSERT(cointriP.vertexPath[P_inputV].size() >= 2);
          SurfacePoint Q_spInputPrev = *----cointriP.vertexPath[P_inputV].end();

          Face Q_fEnd = Q_spInput.edge.halfedge().face();
          if (Q_fEnd == sharedFace(Q_spInput, Q_spInputPrev))
            Q_fEnd = Q_spInput.edge.halfedge().twin().face();

          // Add centroid of Q_fEnd to the vertex path so that the path becomes non-singular
          cointriP.vertexPath[P_inputV].push_back(SurfacePoint(Q_fEnd, Vector3::constant(1./3.)));

          if (Q_fEnd != Q_spOutput.face) {
            // Find face strip between Q_fEnd and Q_spOutput.face
            DLOG_DEBUG(logDepth + 1, "Appending additional face strip to vertexPath (vertex ID: {})", cointriP.uniqueID_per_Vertex[P_intrinsicV]);
            std::vector<Halfedge> Q_inputFaceStrip = findFaceStrip(Q_fEnd, Q_spOutput.face);

            // Add edge midpoints of the face strip to the vertex path, followed by the actual vertex location
            for (Halfedge Q_inputHE : Q_inputFaceStrip)
              cointriP.vertexPath[P_inputV].push_back(SurfacePoint(Q_inputHE.edge(), 0.5));
            cointriP.vertexPath[P_inputV].push_back(Q_spOutput);
          }
        }
      }

      Matrix2d& t = cointriP.temporalTransportOperator[P_inputV];
      const std::vector<SurfacePoint>& Q_vertexPath = cointriP.vertexPath[P_inputV];
      const size_t n = Q_vertexPath.size();

      // If the first end is an edge point, and its halfedge().face() is not adjacent to the next point in the path, apply the transport operator
      SurfacePoint Q_spStart0 = Q_vertexPath[0];
      if (Q_spStart0.type == SurfacePointType::Edge) {
        SurfacePoint Q_spStart1 = Q_vertexPath[1];
        Face Q_f = sharedFace(Q_spStart0, Q_spStart1);
        CIT_ASSERT(Q_f != Face());

        if (Q_spStart0.edge.halfedge().face() != Q_f)
          t = cointriQ.mdata->transportOperator[Q_spStart0.edge];
      }

      for (size_t k = 1; k < n - 1; ++k) {
        SurfacePoint Q_spCurr = Q_vertexPath[k];
        SurfacePoint Q_spPrev = Q_vertexPath[k - 1];

        if (Q_spCurr.type != SurfacePointType::Edge) {
          CIT_ASSERT(Q_spCurr.type == SurfacePointType::Face);
          continue;
        }
        Edge Q_eCurr = Q_spCurr.edge;

        Face Q_f = sharedFace(Q_spCurr, Q_spPrev);
        CIT_ASSERT(Q_f != Face());

        bool reversed = Q_eCurr.halfedge().twin().face() == Q_f;

        Matrix2d tNext = mdataQ.transportOperator[Q_eCurr];
        if (reversed)
          tNext = Matrix2d(tNext).inverse();  // mat = mat.inverse() won't work
        t = tNext * t;
      }

      // If the last end is an edge point, and its halfedge().face() is not adjacent to the previous point in the path, apply the transport operator
      SurfacePoint Q_spEnd0 = Q_vertexPath[n - 1];
      if (Q_spEnd0.type == SurfacePointType::Edge) {
        SurfacePoint Q_spEnd1 = Q_vertexPath[n - 2];
        Face Q_f = sharedFace(Q_spEnd0, Q_spEnd1);
        CIT_ASSERT(Q_f != Face());

        if (Q_spEnd0.edge.halfedge().face() != Q_f)
          t = mdataQ.transportOperator[Q_spEnd0.edge].inverse() * t;     // Note the need for inversion
      }
    }
  };
  computeTemporalTransportOperator_inner(config.cointriA, config.cointriB);
  computeTemporalTransportOperator_inner(config.cointriB, config.cointriA);

  DLOG_INFO(logDepth + 1, "assemblying sparse matrix ... ");

  const size_t N = 2 * (mdataA.nV + mdataB.nV);
  config.T = SparseMatrixd(N, N);
  config.Tinv = SparseMatrixd(N, N);
  {
    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<Eigen::Triplet<double>> tripletsInv;
    auto fillTriplets = [&triplets, &tripletsInv](const CoInTri& cointri, size_t offset) {
      for (Vertex v : cointri.signpostTri->inputMesh.vertices()) {
        const Matrix2d& temporalTransportOperator = cointri.temporalTransportOperator[v];
        Matrix2d temporalTransportOperatorInverse = temporalTransportOperator.inverse();

        size_t i = v.getIndex();

        for (size_t drow = 0; drow < 2; ++drow)
        for (size_t dcol = 0; dcol < 2; ++dcol)
        {
          size_t row = offset + 2 * i + drow;
          size_t col = offset + 2 * i + dcol;
          triplets.emplace_back(row, col, temporalTransportOperator(drow, dcol));
          tripletsInv.emplace_back(row, col, temporalTransportOperatorInverse(drow, dcol));
        }
      }
    };
    fillTriplets(config.cointriA, 0);
    fillTriplets(config.cointriB, 2 * mdataA.nV);
    config.T.setFromTriplets(triplets.begin(), triplets.end());
    config.Tinv.setFromTriplets(tripletsInv.begin(), tripletsInv.end());
  }
}

void cit::computePreconditioner(Configuration& config, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computePreconditioner");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computePreconditioner"); });

  if (!config.laplacianTransportOperatorComputed)
    computeLaplacianTransportOperator(config, logDepth + 1);

  const size_t A_nV = mdataA.nV;
  const size_t B_nV = mdataB.nV;
  const size_t N = 2 * (A_nV + B_nV);

  // for converting barycentric differential vector to actual 2D vector with comparable length
  config.S = SparseMatrixd(N, N);
  {
    std::vector<Eigen::Triplet<double>> triplets;
    auto fillTriplets = [&triplets](const CoInTri& cointriP, const CoInTri& cointriQ, size_t offset) {
      for (size_t i = 0; i < cointriP.mdata->nV; ++i) {
        Vertex P_intrinsicV = cointriP.signpostTri->intrinsicMesh->vertex(i);
        Vertex Q_intrinsicV = cointriP.correspondingVertex[P_intrinsicV];

        // Skip anchor vertex
        if (Q_intrinsicV.getIndex() < cointriQ.mdata->nV) {
          CIT_ASSERT(cointriP.mdata->anchorVertexIDs.count(i));
          CIT_ASSERT(cointriQ.mdata->anchorVertexIDs.count(Q_intrinsicV.getIndex()));
          continue;
        }

        SurfacePoint Q_inputSP = cointriQ.signpostTri->vertexLocations[Q_intrinsicV];

        // An edge point must be converted to a face point so that its barycentric coordinates can be defined
        if (Q_inputSP.type == SurfacePointType::Edge)
          Q_inputSP = convertEdgePointToFacePoint(Q_inputSP);

        CIT_ASSERT(Q_inputSP.type == SurfacePointType::Face);

        Face Q_inputF = Q_inputSP.face;
        std::array<Halfedge, 2> Q_inputHE = {
          Q_inputF.halfedge().next().next(),
          Q_inputF.halfedge().next()
        };
        std::array<Vector2, 2> heVec = {
          cointriQ.signpostTri->inputGeom.halfedgeVectorsInFace[Q_inputHE[0]],
          -cointriQ.signpostTri->inputGeom.halfedgeVectorsInFace[Q_inputHE[1]]
        };

        // careful with indexing!
        for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
        {
          int row = offset + 2 * i + j;
          int col = offset + 2 * i + k;
          triplets.emplace_back(row, col, heVec[k][j]);
        }
      }
    };
    fillTriplets(config.cointriA, config.cointriB, 0);
    fillTriplets(config.cointriB, config.cointriA, 2 * mdataA.nV);
    config.S.setFromTriplets(triplets.begin(), triplets.end());
  }

  // connection Laplacian
  config.L = SparseMatrixd(N, N);
  {
    std::vector<Eigen::Triplet<double>> triplets;
    auto fillTriplets = [&triplets](const CoInTri& cointri, size_t offset) {
      for (Edge e : cointri.signpostTri->inputMesh.edges()) {
        const Matrix2d& laplacianTransportOperator = cointri.laplacianTransportOperator[e];
        Matrix2d laplacianTransportOperatorInverse = laplacianTransportOperator.inverse();

        double cotanWeight = cointri.signpostTri->inputGeom.edgeCotanWeights[e];

        size_t i = e.halfedge().tailVertex().getIndex();
        size_t j = e.halfedge().tipVertex().getIndex();

        if (cointri.mdata->anchorVertexIDs.count(i) || cointri.mdata->anchorVertexIDs.count(j))
          continue;

        // (i,i) & (j,j) blocks
        for (size_t d = 0; d < 2; ++d) {
          size_t i_idx = offset + 2 * i + d;
          size_t j_idx = offset + 2 * j + d;
          triplets.emplace_back(i_idx, i_idx, -cotanWeight);
          triplets.emplace_back(j_idx, j_idx, -cotanWeight);
        }

        // (i,j) block
        for (size_t drow = 0; drow < 2; ++drow)
        for (size_t dcol = 0; dcol < 2; ++dcol)
        {
          size_t row = offset + 2 * i + drow;
          size_t col = offset + 2 * j + dcol;
          triplets.emplace_back(row, col, cotanWeight * laplacianTransportOperatorInverse(drow, dcol));    // careful about the orientation!
        }

        // (j,i) block
        for (size_t drow = 0; drow < 2; ++drow)
        for (size_t dcol = 0; dcol < 2; ++dcol)
        {
          size_t row = offset + 2 * j + drow;
          size_t col = offset + 2 * i + dcol;
          triplets.emplace_back(row, col, cotanWeight * laplacianTransportOperator(drow, dcol));    // careful about the orientation!
        }
      }
    };
    fillTriplets(config.cointriA, 0);
    fillTriplets(config.cointriB, 2 * mdataA.nV);
    config.L.setFromTriplets(triplets.begin(), triplets.end());
  }

  // Add ones to diagonal corresponding to anchor vertices, so that the system becomes full-rank
  SparseMatrixd ones(N, N);
  {
    std::vector<Eigen::Triplet<double>> triplets;
    auto A_i = mdataA.anchorVertexIDs.begin();
    auto B_i = mdataB.anchorVertexIDs.begin();
    for (; A_i != mdataA.anchorVertexIDs.end(); ++A_i, ++B_i) {
      CIT_ASSERT(B_i != mdataB.anchorVertexIDs.end());
      for (size_t di = 0; di < 2; ++di) {
        size_t rowA = 2 * (*A_i) + di;
        size_t rowB = 2 * (mdataA.nV + *B_i) + di;
        triplets.emplace_back(rowA, rowA, 1.);
        triplets.emplace_back(rowB, rowB, 1.);
      }
    }
    CIT_ASSERT(B_i == mdataB.anchorVertexIDs.end());
    ones.setFromTriplets(triplets.begin(), triplets.end());
  }

  SparseMatrixd SL = config.S * config.L;
  config.P = SparseMatrixd(SL.transpose()) * massMatrix * SL + ones;

  config.preconditionerComputed = true;
}

VectorXd dbg_descentDirection;
VectorXd cit::computeDescentDirection(Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, int l, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN computeDescentDirection");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   computeDescentDirection"); });

  if (!config.topologyValid)
    return {};

  if (!config.derivativeComputed)
    computeDerivative(config, logDepth + 1);

  if (!config.preconditionerComputed)
    computePreconditioner(config, logDepth + 1);

  double s = 0.5;
  config.smoothedGradient = s * config.gradient;
  config.smoothedHessian = s * config.hessian;
  for (const auto& p : history) {
    s *= 0.5;
    config.smoothedGradient += s * p.first;
    config.smoothedHessian += s * p.second;
  }

  VectorXd descentDirection;

  if (sysParam.firstOrderMode) {
    DLOG_INFO(logDepth + 1, "First-order mode: descent direction is the same as smoothed gradient negated");
    descentDirection = -config.smoothedGradient;

  } else {
    DLOG_INFO(logDepth + 1, "Computing descent direction (wL={}, l={}) ... ", wL, l);

    Eigen::SimplicialLDLT<SparseMatrixd> solver;
    if (wL) {
      if (l)
        solver.compute(config.smoothedHessian + wL * std::pow(sysParam.alpha, l) * config.P);
      else
        solver.compute(config.smoothedHessian + wL * config.P);
    }
    else {
      solver.compute(config.smoothedHessian);
    }
    if(solver.info() != Eigen::Success) throw std::runtime_error("decomposition failed!");

    descentDirection = solver.solve(-config.smoothedGradient);
    if(solver.info() != Eigen::Success) throw std::runtime_error("solving failed!");
  }

  // Sanity check that we have zero at anchor vertices
  auto A_i = mdataA.anchorVertexIDs.begin();
  auto B_i = mdataB.anchorVertexIDs.begin();
  for (; A_i != mdataA.anchorVertexIDs.end(); ++A_i, ++B_i) {
    CIT_ASSERT(B_i != mdataB.anchorVertexIDs.end());
    for (size_t di = 0; di < 2; ++di) {
      CIT_ASSERT(descentDirection[2 * (*A_i) + di] == 0.);
      CIT_ASSERT(descentDirection[2 * (mdataA.nV + *B_i) + di] == 0.);
    }
  }
  CIT_ASSERT(B_i == mdataB.anchorVertexIDs.end());

  DLOG_INFO(logDepth + 1, "  Done, norm = {}", descentDirection.norm());
  return (dbg_descentDirection = descentDirection);
}

namespace cit { namespace detail {
// Utility struct for computing energy for a diamond before/after flipping edge
struct SplitDiamond {
  std::array<std::array<Vector2d, 3>, 2> A_triangles;
  std::array<std::array<Vector2d, 3>, 2> B_triangles;
  SplitDiamond(const std::array<Vector2, 4>& A_diamond, const std::array<Vector2, 4>& B_diamond, bool flipped) {
    std::function<std::array<std::array<Vector2d, 3>, 2>(const std::array<Vector2, 4>&)> getTriangles;
    if (flipped) {
      getTriangles = [](const std::array<Vector2, 4>& diamond)->std::array<std::array<Vector2d, 3>, 2> {
        /*
             1-----0
            / \   /
           /   \ /
          2-----3
        */
        return { {
          {diamond[1].castTo<Vector2d>(), diamond[2].castTo<Vector2d>(), diamond[3].castTo<Vector2d>()},
          {diamond[0].castTo<Vector2d>(), diamond[1].castTo<Vector2d>(), diamond[3].castTo<Vector2d>()}
        } };
      };
    } else {
      getTriangles = [](const std::array<Vector2, 4>& diamond)->std::array<std::array<Vector2d, 3>, 2> {
        /*
          1-----0
           \   / \
            \ /   \
             2-----3
        */
        return { {
          {diamond[0].castTo<Vector2d>(), diamond[1].castTo<Vector2d>(), diamond[2].castTo<Vector2d>()},
          {diamond[0].castTo<Vector2d>(), diamond[2].castTo<Vector2d>(), diamond[3].castTo<Vector2d>()}
        } };
      };
    }
    A_triangles = getTriangles(A_diamond);
    B_triangles = getTriangles(B_diamond);
  }
};
double computeDiamondEnergy(const SplitDiamond& sd) {
  double res = 0;
  for (int i = 0; i < 2; ++i) {
    res += computeTriangleEnergy(sd.A_triangles[i], sd.B_triangles[i]);
  }
  return res;
}

} // namespace detail
} // namespace cit

void cit::flipToMinimizeEnergy(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipToMinimizeEnergy");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipToMinimizeEnergy"); });

  double initial_energy = computeEnergy(config, logger, logDepth + 1), energy = initial_energy;
  LOG_INFO(logger, logDepth + 1, "Initial energy: {}, flipping edges ...", initial_energy);

  // as per SignpostIntrinsicTriangulation::flipToDelaunay
  std::deque<Edge> edgesToCheck;
  EdgeData<char> inQueue(*config.cointriA.signpostTri->intrinsicMesh, true);
  for (Edge A_e : config.cointriA.signpostTri->intrinsicMesh->edges()) {
    edgesToCheck.push_back(A_e);
  }
  size_t nFlips = 0;
  while (!edgesToCheck.empty()) {
    Edge A_e = edgesToCheck.front();
    edgesToCheck.pop_front();
    if (!inQueue[A_e]) continue;
    inQueue[A_e] = false;

    // get corresponding edge on B and unique edge ID
    Edge B_e = config.cointriA.correspondingEdge[A_e];

    // get diamond shape
    Halfedge A_he = A_e.halfedge();
    Halfedge B_he = config.cointriA.correspondingHalfedge[A_he];
    std::array<Vector2, 4> A_diamond = getLayoutDiamond(*config.cointriA.signpostTri, A_he);
    std::array<Vector2, 4> B_diamond = getLayoutDiamond(*config.cointriB.signpostTri, B_he);

    // compute energy for the two configurations
    detail::SplitDiamond diamond_before(A_diamond, B_diamond, false);
    detail::SplitDiamond diamond_after (A_diamond, B_diamond, true );
    double diamond_energy_before = detail::computeDiamondEnergy(diamond_before);
    double diamond_energy_after  = detail::computeDiamondEnergy(diamond_after);

    if (diamond_energy_after >= diamond_energy_before) continue;

    // try flipping
    bool wasFlipped = false;
    if (isEdgeFlippable(config.cointriA, A_e) && isEdgeFlippable(config.cointriB, B_e)) {
      config.cointriA.signpostTri->flipEdgeIfPossible(A_e);
      config.cointriB.signpostTri->flipEdgeIfPossible(B_e);
      wasFlipped = true;

      // update energy
      double energy_delta = diamond_energy_after - diamond_energy_before;
      CIT_ASSERT(energy_delta < 0);
      energy += energy_delta;
    }

    if (!wasFlipped) continue;

    // Handle the aftermath of a flip
    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge A_heN = A_he.next();
    Halfedge A_heT = A_he.twin();
    Halfedge A_heTN = A_heT.next();
    std::vector<Edge> neighEdges = {A_heN.edge(), A_heN.next().edge(), A_heTN.edge(), A_heTN.next().edge()};
    for (Edge nE : neighEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }

  computeEnergy(config, logger, logDepth + 1);

  LOG_INFO(logger, logDepth + 1, "Flipped {} edges", nFlips);
  LOG_INFO(logger, logDepth + 1, "  final energy: {}", energy);
  LOG_INFO(logger, logDepth + 1, "  (sanity check: {}, error: {})", config.energy, std::fabs(energy - config.energy));
  LOG_INFO(logger, logDepth + 1, "  (delta: {})", energy - initial_energy);

  config.cointriA.signpostTri->refreshQuantities();
  config.cointriB.signpostTri->refreshQuantities();
  updateCorrespondence(config, logger, logDepth + 1);
}

namespace cit { namespace detail {
double getDiamondMinAngle(const std::array<Vector2, 4>& diamond, bool flipped) {
  if (flipped) {
    /*
       1-----0
      / \   /
     /   \ /
    2-----3
    */
    double angle231 = angle(diamond[2] - diamond[3], diamond[1] - diamond[3]);
    double angle312 = angle(diamond[3] - diamond[1], diamond[2] - diamond[1]);
    double angle130 = angle(diamond[1] - diamond[3], diamond[0] - diamond[3]);
    double angle013 = angle(diamond[0] - diamond[1], diamond[3] - diamond[1]);
    double angle123 = angle(diamond[1] - diamond[2], diamond[3] - diamond[2]);
    double angle301 = angle(diamond[3] - diamond[0], diamond[1] - diamond[0]);
    return std::min<double>({angle231, angle312, angle130, angle013, angle123, angle301});

  } else {
    /*
    1-----0
     \   / \
      \ /   \
       2-----3
    */
    double angle120 = angle(diamond[1] - diamond[2], diamond[0] - diamond[2]);
    double angle201 = angle(diamond[2] - diamond[0], diamond[1] - diamond[0]);
    double angle023 = angle(diamond[0] - diamond[2], diamond[3] - diamond[2]);
    double angle302 = angle(diamond[3] - diamond[0], diamond[2] - diamond[0]);
    double angle012 = angle(diamond[0] - diamond[1], diamond[2] - diamond[1]);
    double angle230 = angle(diamond[2] - diamond[3], diamond[0] - diamond[3]);
    return std::min<double>({angle120, angle201, angle023, angle302, angle012, angle230});
  }
};
} // namespace cit
} // namespace detail

void cit::flipToMaximizeMinAngle(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN flipToMaximizeMinAngle");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   flipToMaximizeMinAngle"); });

  double A_minAngle = getMinAngle(config.cointriA);
  double B_minAngle = getMinAngle(config.cointriB);

  double initial_minAngle = std::min<double>(A_minAngle, B_minAngle);
  double initial_energy = computeEnergy(config, logger, logDepth + 1);

  LOG_INFO(logger, logDepth + 1, "Flipping edges to maximize minimum angle ... ");
  LOG_INFO(logger, logDepth + 1, "  initial min angle: {} (A: {}, B: {}), flipping edges ...", initial_minAngle, A_minAngle, B_minAngle);
  LOG_INFO(logger, logDepth + 1, "  initial energy: {}", initial_energy);

  // as per SignpostIntrinsicTriangulation::flipToDelaunay
  std::deque<Edge> edgesToCheck;
  EdgeData<char> inQueue(*config.cointriA.signpostTri->intrinsicMesh, true);
  for (Edge A_e : config.cointriA.signpostTri->intrinsicMesh->edges()) {
    edgesToCheck.push_back(A_e);
  }
  size_t nFlips = 0;
  while (!edgesToCheck.empty()) {
    Edge A_e = edgesToCheck.front();
    edgesToCheck.pop_front();
    if (!inQueue[A_e]) continue;
    inQueue[A_e] = false;

    Edge B_e = config.cointriA.correspondingEdge[A_e];

    if (!isEdgeFlippable(config.cointriA, A_e)) continue;
    if (!isEdgeFlippable(config.cointriB, B_e)) continue;

    // get diamond shape
    Halfedge A_he = A_e.halfedge();
    Halfedge B_he = config.cointriA.correspondingHalfedge[A_he];
    std::array<Vector2, 4> A_diamond = getLayoutDiamond(*config.cointriA.signpostTri, A_he);
    std::array<Vector2, 4> B_diamond = getLayoutDiamond(*config.cointriB.signpostTri, B_he);

    double A_minAngle_before = detail::getDiamondMinAngle(A_diamond, false);
    double B_minAngle_before = detail::getDiamondMinAngle(B_diamond, false);
    double minAngle_before = std::min<double>(A_minAngle_before, B_minAngle_before);

    double A_minAngle_after = detail::getDiamondMinAngle(A_diamond, true);
    double B_minAngle_after = detail::getDiamondMinAngle(B_diamond, true);
    double minAngle_after = std::min<double>(A_minAngle_after, B_minAngle_after);

    // try flipping
    bool wasFlipped = false;
    if (minAngle_after > minAngle_before) {
      config.cointriA.signpostTri->flipEdgeIfPossible(A_e);
      config.cointriB.signpostTri->flipEdgeIfPossible(B_e);
      wasFlipped = true;
    }

    if (!wasFlipped) continue;

    // Handle the aftermath of a flip
    nFlips++;

    // Add neighbors to queue, as they may need flipping now
    Halfedge A_heN = A_he.next();
    Halfedge A_heT = A_he.twin();
    Halfedge A_heTN = A_heT.next();
    std::vector<Edge> neighEdges = {A_heN.edge(), A_heN.next().edge(), A_heTN.edge(), A_heTN.next().edge()};
    for (Edge nE : neighEdges) {
      if (!inQueue[nE]) {
        edgesToCheck.push_back(nE);
        inQueue[nE] = true;
      }
    }
  }

  config.cointriA.signpostTri->refreshQuantities();
  config.cointriB.signpostTri->refreshQuantities();
  updateCorrespondence(config, logger, logDepth + 1);

  A_minAngle = getMinAngle(config.cointriA);
  B_minAngle = getMinAngle(config.cointriB);

  double final_minAngle = std::min<double>(A_minAngle, B_minAngle);
  double final_energy = computeEnergy(config, logger, logDepth + 1);

  LOG_INFO(logger, logDepth + 1, "Flipped {} edges", nFlips);
  LOG_INFO(logger, logDepth + 1, "  final min angle: {} (A: {}, B: {})", final_minAngle, A_minAngle, B_minAngle);
  LOG_INFO(logger, logDepth + 1, "    delta: {}", final_minAngle - initial_minAngle);
  LOG_INFO(logger, logDepth + 1, "  final energy: {}", final_energy);
  LOG_INFO(logger, logDepth + 1, "    delta: {}", final_energy - initial_energy);
}

#if 0
#include <nlopt.hpp>

namespace {

//--------------------------------------------
// >>>> structs and functions for use in NLopt
struct CircleConstraint {
  Vector2 center;
  double radius;
};

struct LineConstraint {
  Vector2 point;
  Vector2 normal;
};

using ConstraintTuple = std::tuple<CircleConstraint, LineConstraint, LineConstraint>;

struct ObjectiveFuncData {
  std::vector<ConstraintTuple>* constraints = nullptr;
  std::vector<double>* x_out = nullptr;
  nlopt::opt* opt = nullptr;
};

double circleConstraintFunc(const std::vector<double> &x_, std::vector<double> &grad, void *dataPtr) {
  Vector2 x = Vector2::castFrom(x_);
  const CircleConstraint& data = *(CircleConstraint*)(dataPtr);
  // Equation:
  //    |x - c|^2 <= r^2
  //    (x - c).(x - c) - r^2 <= 0
  if (!grad.empty()) {
    grad = (2. * (x - data.center)).castTo<std::vector<double>>();
  }
  return norm2(x - data.center) - data.radius * data.radius;
};

double lineConstraintFunc(const std::vector<double> &x_, std::vector<double> &grad, void *dataPtr) {
  Vector2 x = Vector2::castFrom(x_);
  const LineConstraint& data = *(LineConstraint*)(dataPtr);
  // Equation:
  //    (x - p).n >= 0
  //    -x.n + p.n <= 0
  if (!grad.empty()) {
    grad = (-data.normal).castTo<std::vector<double>>();
  }
  return -dot(x, data.normal) + dot(data.point, data.normal);
};

bool isFeasible(const std::vector<double>& x, const std::vector<ConstraintTuple>& constraints) {
  std::vector<double> dummy;
  for (auto& t : constraints) {
    if (circleConstraintFunc(x, dummy, (void*)&std::get<0>(t)) >= 0) return false;
    if (lineConstraintFunc(x, dummy, (void*)&std::get<1>(t)) >= 0) return false;
    if (lineConstraintFunc(x, dummy, (void*)&std::get<2>(t)) >= 0) return false;
  }
  return true;
}

double objectiveFunc(const std::vector<double> &x_, std::vector<double> &grad, void *dataPtr) {
  Vector2 x = Vector2::castFrom(x_);
  // Function: x^2
  if (!grad.empty()) {
      grad = (2. * x).castTo<std::vector<double>>();
  }

  const ObjectiveFuncData& data = *(ObjectiveFuncData*)(dataPtr);
  if (isFeasible(x_, *data.constraints)) {
    *data.x_out = x_;
    data.opt->force_stop();
  }

  return x.norm2();
};
// <<<< structs and functions for use in NLopt
//--------------------------------------------

} // namespace

size_t cit::relocateVerticesToRemoveDegenerateFaces(Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  LOG_TRACE(logger, logDepth, "BEGIN relocateVerticesToRemoveDegenerateFaces");
  auto scopeExit = kt84::make_ScopeExit([&](){ LOG_TRACE(logger, logDepth, "END   relocateVerticesToRemoveDegenerateFaces"); });

  std::set<std::tuple<int, double, double>> res;

  auto relocateVerticesToRemoveDegenerateFaces_inner = [&logger, &logDepth, &res](CoInTri& cointri) {

    for (Vertex v : cointri.signpostTri->intrinsicMesh->vertices()) {
      // Skip all fixed vertices
      if (v.getIndex() < cointri.mdata->nV)
        continue;

      // Skip vertices inserted on input edges
      if (cointri.signpostTri->vertexLocations[v].type == SurfacePointType::Edge)
        continue;

      kt84::MaxMinAverage angleBefore;
      for (Halfedge he : v.outgoingHalfedges()) {
        for (Halfedge he2 : {he.next().next(), he.twin().next().next()}) {
          double angle = cointri.signpostTri->cornerAngle(he2.corner());
          angleBefore.update(angle);
        }
      }
      if (angleBefore.min() > sysParam.angleThreshold)
        continue;

      LOG_DEBUG(logger, logDepth + 1, "Bad angle {} detected around vertex {}", angleBefore.min(), cointri.uniqueID_per_Vertex[v]);

      // Figure out a new position closest to the current one while satisfying all the angle constraints, using NLopt
      const double theta = 1.05 * sysParam.angleThreshold;

      std::vector<ConstraintTuple> constraints(v.degree());
      std::vector<double> x_out;
      nlopt::opt opt(nlopt::LD_MMA, 2);

      ObjectiveFuncData data = {
        &constraints,
        &x_out,
        &opt
      };

      opt.set_min_objective(objectiveFunc, &data);
      opt.set_maxeval(1000);

      Halfedge he = v.halfedge();
      for (size_t i = 0; i < v.degree(); ++i, he = he.twin().next()) {
        Vector2 pA = cointri.signpostTri->halfedgeVectorsInVertex[he];
        Vector2 pB = cointri.signpostTri->halfedgeVectorsInVertex[he.next().next().twin()];
        Vector2 pM = 0.5 * (pA + pB);

        double L = norm(pB - pA);
        double R = 0.5 * L / std::sin(theta);
        Vector2 dAB = (pB - pA) / L;
        Vector2 dMC = {-dAB.y, dAB.x};

        Vector2 pC = pM + R * std::cos(theta) * dMC;

        CircleConstraint& circleConstraint = std::get<0>(constraints[i]);
        circleConstraint.center = pC;
        circleConstraint.radius = R;
        opt.add_inequality_constraint(circleConstraintFunc, &circleConstraint);

        LineConstraint& lineConstraint1 = std::get<1>(constraints[i]);
        lineConstraint1.point = pA;
        lineConstraint1.normal = dMC * Vector2::fromAngle(theta);
        opt.add_inequality_constraint(lineConstraintFunc, &lineConstraint1);

        LineConstraint& lineConstraint2 = std::get<2>(constraints[i]);
        lineConstraint2.point = pB;
        lineConstraint2.normal = dMC * Vector2::fromAngle(-theta);
        opt.add_inequality_constraint(lineConstraintFunc, &lineConstraint2);
      }

      try {
        std::vector<double> x = {0, 0};
        double minf;
        nlopt::result opt_res = opt.optimize(x, minf);
        LOG_DEBUG(logger, logDepth + 1, "nlopt return code: {}, numevals: {}", nlopt_result_to_string((nlopt_result)opt_res), opt.get_numevals());

      } catch (const nlopt::forced_stop& e) {
        LOG_DEBUG(logger, logDepth + 1, "nlopt forced stop, numevals: {}", opt.get_numevals());

      } catch(const std::exception &e) {
        LOG_WARN(logger, logDepth + 1, "nlopt failed: {}, numevals", e.what(), opt.get_numevals());
      }

      if (x_out.empty()) {
        LOG_DEBUG(logger, logDepth + 1, "Couldn't find a feasible solution, skipping");
        continue;
      }

      if (!isFeasible(x_out, constraints))
        throw std::logic_error("This shouldn't happen...");

      // Trace geodesic
      Vector2 traceVec = Vector2::castFrom(x_out);
      if (traceVec == Vector2::zero())
        throw std::logic_error("This shouldn't happen...");

      TraceGeodesicResult traceResult = traceGeodesic(*cointri.signpostTri, {v}, traceVec);

      // Try to relocate the vertex to the new position
      int u = cointri.uniqueID_per_Vertex[v];
      if (!cointri.signpostTri->relocateInsertedVertex(v, traceResult.endPoint)) {
        LOG_WARN(logger, logDepth + 1, "Failed to relocate vertex (id: {}) with bad corner", u);
        continue;
      }

      kt84::MaxMinAverage angleAfter;
      for (Halfedge he : v.outgoingHalfedges()) {
        for (Halfedge he2 : {he.next().next(), he.twin().next().next()}) {
          double angle = cointri.signpostTri->cornerAngle(he2.corner());
          angleAfter.update(angle);
        }
      }

      LOG_DEBUG(logger, logDepth + 1, "Relocated to remove the bad angle: {}", angleAfter.min());

      res.emplace(u, traceVec.norm(), angleAfter.min() - angleBefore.min());
    }
  };

  relocateVerticesToRemoveDegenerateFaces_inner(config.cointriA);
  relocateVerticesToRemoveDegenerateFaces_inner(config.cointriB);

  if (!res.empty()) {
    LOG_INFO(logger, logDepth + 1, "===== Relocated {} vertices to remove degenerate corners =====", res.size());
    for (auto& t : res) {
      int u;
      double dist, delta;
      std::tie(u, dist, delta) = t;
      LOG_INFO(logger, logDepth + 1, "======= processed vertex: {} (dist: {}, delta: {}) =======", u, dist, delta);
    }
  }
  return res.size();
}
#endif

VectorXsp cit::displaceSolutionVector(const VectorXsp& x, const VectorXd& delta, VertexData<std::vector<SurfacePoint>>& A_vertexPath, VertexData<std::vector<SurfacePoint>>& B_vertexPath, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN displaceSolutionVector");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   displaceSolutionVector"); });

  CIT_ASSERT((size_t)x.size() == mdataA.nV + mdataB.nV);
  CIT_ASSERT(delta.size() == 2 * x.size());
  VectorXsp res(x.size());

  A_vertexPath = VertexData<std::vector<SurfacePoint>>(*mdataA.mesh);
  B_vertexPath = VertexData<std::vector<SurfacePoint>>(*mdataB.mesh);

  auto displaceSolutionVector_inner = [&x, &delta, &res, &logDepth](const ModelData& mdataP, const ModelData& mdataQ, VertexData<std::vector<SurfacePoint>>& P_vertexPath, size_t offset, size_t& nSnapped) {
    for (size_t i = offset; i < offset + mdataP.nV; ++i) {
      // Skip anchor vertex
      if (x[i].type == SurfacePointType::Vertex) {
        res[i] = x[i];
        continue;
      }

      Face f = x[i].face;

      // For an edge point, convert it to a face point
      if (f == Face())
        f = convertEdgePointToFacePoint(x[i]).face;

      std::array<Vector2, 3> vertCoords = {{
        {0., 0.},
        mdataQ.geometry->halfedgeVectorsInFace[f.halfedge()],
        -mdataQ.geometry->halfedgeVectorsInFace[f.halfedge().next().next()]
      }};
      Vector3 baryCoords = {
        delta[2*i],
        delta[2*i + 1],
        -(delta[2*i] + delta[2*i + 1])
      };
      Vector2 traceVec = {0,0};
      for (int j = 0; j < 3; ++j)
        traceVec += baryCoords[j] * vertCoords[j];

      // Something is wrong when traceVec is infinitesimal
      CIT_ASSERT(norm(traceVec) > 0.);

      // Transform traceVec to the edge's local coordinate system
      if (x[i].type == SurfacePointType::Edge)
        traceVec /= mdataQ.geometry->halfedgeVectorsInFace[x[i].edge.halfedge()].normalize();

      TraceOptions options;
      options.includePath = true;
      TraceGeodesicResult traceResult = traceGeodesic(*mdataQ.geometry, x[i], traceVec, options);

      // Snap the endPoint if close to an edge
      CIT_ASSERT(traceResult.endPoint.type == SurfacePointType::Face);
      if (snapFacePointToEdge(traceResult.endPoint)) {
        traceResult.pathPoints.push_back(traceResult.endPoint);
        DLOG_DEBUG(logDepth + 1, "Snapped {}'s vertex {} onto {}'s edge {}", mdataP.name, i - offset, mdataQ.name, traceResult.endPoint.edge);
        ++nSnapped;
      }

      res[i] = traceResult.endPoint;
      CIT_ASSERT(res[i].type != SurfacePointType::Vertex);

      Vertex P_v = mdataP.mesh->vertex(i - offset);
      P_vertexPath[P_v] = traceResult.pathPoints;
    }
  };
  size_t A_nSnapped = 0;
  size_t B_nSnapped = 0;
  displaceSolutionVector_inner(mdataA, mdataB, A_vertexPath, 0        , A_nSnapped);
  displaceSolutionVector_inner(mdataB, mdataA, B_vertexPath, mdataA.nV, B_nSnapped);

  DLOG_INFO(logDepth + 1, "Snapped {} vertices to edges ({} from A, {} from B)", A_nSnapped + B_nSnapped, A_nSnapped, B_nSnapped);

  return res;
}

bool cit::createConfiguration(uint64_t id, VectorXsp x, Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  DLOG_INFO(logDepth, "BEGIN createConfiguration (insertion)");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_INFO(logDepth, "END   createConfiguration (insertion)"); });

  if (config.id) {
    CIT_ASSERT(config.id == id);
    LOG_INFO(logger, logDepth, "");
    LOG_INFO(logger, logDepth, "Falling back to vertex insertion mode ----------------------------------------");

  } else {
    config.id = id;
    config.cointriA.mdata = &mdataA;
    config.cointriB.mdata = &mdataB;
    LOG_INFO(logger, logDepth, "");
    LOG_INFO(logger, logDepth, "createConfiguration BEGIN ----------------------------------------");
  }
  config.topologyValid = false;
  // config.angleValid = false;
  kt84::Timer timer;

  // Initialize intrinsic triangulation
  config.cointriA.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataA.mesh, *mdataA.geometry));
  config.cointriB.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataB.mesh, *mdataB.geometry));

  // Insert verticies to each other
  CIT_ASSERT((size_t)x.size() == mdataA.nV + mdataB.nV);
  insertVertices(x.head(mdataA.nV), config.cointriA, config.cointriB, logger, logDepth + 1);
  insertVertices(x.tail(mdataB.nV), config.cointriB, config.cointriA, logger, logDepth + 1);
  CIT_ASSERT(config.cointriA.signpostTri->intrinsicMesh->nVertices() == config.cointriB.signpostTri->intrinsicMesh->nVertices());
  CIT_ASSERT(config.cointriA.signpostTri->intrinsicMesh->nEdges() == config.cointriB.signpostTri->intrinsicMesh->nEdges());
  CIT_ASSERT(config.cointriA.signpostTri->intrinsicMesh->nFaces() == config.cointriB.signpostTri->intrinsicMesh->nFaces());

  // Get initial roughly compatible connectivity via Delaunay flipping
  config.cointriA.signpostTri->flipToDelaunay([&config](Edge A_e) { return isEdgeFlippable(config.cointriA, A_e); });
  config.cointriB.signpostTri->flipToDelaunay([&config](Edge B_e) { return isEdgeFlippable(config.cointriB, B_e); });

  updateCorrespondence(config, logger, logDepth + 1);
  mergeNearbyVertexPairs(config, logger, logDepth + 1);

  // Flip edges as much as possible
  flipToCompatible(config, logger, logDepth + 1);

  flipPatchInteriorCompatibleEdgesToExploreVariations(config, logger, logDepth + 1);

  mergeInconsistentVertexPairs(config, false, logger, logDepth + 1);
  flipToCompatible(config, logger, logDepth + 1);

  mergeInconsistentVertexPairs(config, true, logger, logDepth + 1);
  flipToCompatible(config, logger, logDepth + 1);

  convexifyConcaveIncompatiblePatch(config, logger, logDepth + 1);
  flipToCompatible(config, logger, logDepth + 1);

  bool success = config.nCompatibleEdges == config.cointriA.signpostTri->intrinsicMesh->nEdges();

  if (success) {
    CIT_ASSERT(config.nCompatibleFaces == config.cointriA.signpostTri->intrinsicMesh->nFaces());

    // // Remove degeneracy also before splitMergedVertices, to increase the chance of success
    // flipToRemoveDegenerateFaces(config, logger, logDepth + 1);
    // for (int i = 0; i < 100; ++i) {
    //   if (!relocateVerticesToRemoveDegenerateFaces(config, logger, logDepth + 1))
    //     break;
    // }

    success = splitMergedVertices(config, logger, logDepth + 1);
  }

  if (success) {
    CIT_ASSERT(config.nCompatibleEdges == config.cointriA.signpostTri->intrinsicMesh->nEdges());
    CIT_ASSERT(config.nCompatibleFaces == config.cointriA.signpostTri->intrinsicMesh->nFaces());

    config.topologyValid = true;

    // flipToMinimizeEnergy(config, logger, logDepth + 1);
    flipToMaximizeMinAngle(config, logger, logDepth + 1);

    // for (int i = 0; i < 100; ++i) {
    //   if (!relocateVerticesToRemoveDegenerateFaces(config, logger, logDepth + 1))
    //     break;
    // }

    double A_minAngle = getMinAngle(config.cointriA);
    double B_minAngle = getMinAngle(config.cointriB);
    double minAngle = std::min<double>(A_minAngle, B_minAngle);

    LOG_INFO(logger, logDepth, "");
    LOG_INFO(logger, logDepth, "Minimum angle: {}", minAngle);
    LOG_INFO(logger, logDepth, "  A: {}", A_minAngle);
    LOG_INFO(logger, logDepth, "  B: {}", B_minAngle);

    // config.angleValid = minAngle > sysParam.angleThreshold;

    makeIntrinsicEdgesConsistentlyOriented(config.cointriA, config.cointriB);
    makeIntrinsicEdgesConsistentlyOriented(config.cointriB, config.cointriA);
    updateCorrespondence(config, logger, logDepth + 1);
    config.cointriA.signpostTri->refreshQuantities();
    config.cointriB.signpostTri->refreshQuantities();

    computeEnergy(config, logger, logDepth + 1);

    LOG_INFO(logger, logDepth, "Successfully generated a CIT: ");
    // if (config.angleValid) {
    //   LOG_INFO(logger, logDepth, "Successfully generated a CIT: ");
    // } else {
    //   LOG_INFO(logger, logDepth, "Generated a topologycally valid CIT, but with bad corner angle: {} ", minAngle);
    // }
    LOG_INFO(logger, logDepth, "  # vertices: {}", config.cointriA.signpostTri->intrinsicMesh->nVertices());
    LOG_INFO(logger, logDepth, "  # edges: {}", config.cointriA.signpostTri->intrinsicMesh->nEdges());
    LOG_INFO(logger, logDepth, "  # faces: {}", config.cointriA.signpostTri->intrinsicMesh->nFaces());
    LOG_INFO(logger, logDepth, "  Energy: {}", config.energy);

  } else {
    LOG_INFO(logger, logDepth, "Failed to generate a CIT! this solution is infeasible.");
  }

  LOG_INFO(logger, logDepth, "  Elapsed time (ms): {}", timer.milliseconds_str());

  setSolutionVector(config);

  return config.topologyValid;// && config.angleValid;
}

bool cit::createConfiguration(uint64_t id, const Configuration& prevConfig, const VectorXd& delta, Configuration& config, spdlog::logger_ptr logger, size_t logDepth) {
  DLOG_INFO(logDepth, "BEGIN createConfiguration (relocation)");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_INFO(logDepth, "END   createConfiguration (relocation)"); });

  kt84::Timer timer;

  config = {};

  config.id = id;
  config.previous_id = prevConfig.id;
  config.cointriA.mdata = &mdataA;
  config.cointriB.mdata = &mdataB;

  // We start with a topologically valid configuration
  config.topologyValid = true;
  // config.angleValid = false;

  copyConfiguration(prevConfig, config, logDepth + 1);

  VectorXsp x = displaceSolutionVector(prevConfig.x, delta, config.cointriA.vertexPath, config.cointriB.vertexPath, logDepth + 1);

// Temporarily disable vertex relocation mode
#if 1
config.id = 0;
return createConfiguration(id, x, config, logger, logDepth + 1);
#endif

  // Try to relocate vertices to x
  auto createConfiguration_inner = [&x, &prevConfig](const CoInTri& cointriP, CoInTri& cointriQ, size_t offset) -> bool {
    for (size_t i = 0; i < cointriP.mdata->nV; ++i) {
      // Skip anchor vertex
      if (cointriP.mdata->anchorVertexIDs.count(i))
        continue;
      CIT_ASSERT(x[offset + i].type != SurfacePointType::Vertex);

      Vertex P_vIntrinsic = cointriP.signpostTri->intrinsicMesh->vertex(i);
      Vertex Q_vIntrinsic = cointriP.correspondingVertex[P_vIntrinsic];
      CIT_ASSERT(cointriQ.signpostTri->vertexLocations[Q_vIntrinsic] == prevConfig.x[offset + i]);

      SurfacePoint Q_pointOnInput = x[offset + i];
      SurfacePoint Q_pointOnIntrinsic = cointriQ.signpostTri->equivalentPointOnIntrinsic(Q_pointOnInput);
      if (cointriQ.signpostTri->relocateInsertedVertex(Q_vIntrinsic, Q_pointOnIntrinsic)) {
        // If the input was an edge point, the relocated position should also be converted to an edge point
        if (Q_pointOnInput.type == SurfacePointType::Edge)
          cointriQ.signpostTri->convertToEdgePoint(Q_vIntrinsic, 1.e-3);

      } else {
        return false;
      }
    }
    return true;
  };
  bool success =
    createConfiguration_inner(config.cointriA, config.cointriB, 0) &&
    createConfiguration_inner(config.cointriB, config.cointriA, mdataA.nV);

  if (success) {
    mergeNearbyVertexPairs(config, logger, logDepth + 1);

    success = splitMergedVertices(config, logger, logDepth + 1);
  }

  if (success) {
    // flipToMinimizeEnergy(config, logger, logDepth + 1);
    flipToMaximizeMinAngle(config, logger, logDepth + 1);

    // for (int i = 0; i < 100; ++i) {
    //   if (!relocateVerticesToRemoveDegenerateFaces(config, logger, logDepth + 1))
    //     break;
    // }

    double A_minAngle = getMinAngle(config.cointriA);
    double B_minAngle = getMinAngle(config.cointriB);
    double minAngle = std::min<double>(A_minAngle, B_minAngle);

    LOG_INFO(logger, logDepth, "");
    LOG_INFO(logger, logDepth, "Minimum angle: {}", minAngle);
    LOG_INFO(logger, logDepth, "  A: {}", A_minAngle);
    LOG_INFO(logger, logDepth, "  B: {}", B_minAngle);

    // success = minAngle > sysParam.angleThreshold;
  }

  if (success) {
    // config.angleValid = true;

    makeIntrinsicEdgesConsistentlyOriented(config.cointriA, config.cointriB);
    makeIntrinsicEdgesConsistentlyOriented(config.cointriB, config.cointriA);
    updateCorrespondence(config, logger, logDepth + 1);
    config.cointriA.signpostTri->refreshQuantities();
    config.cointriB.signpostTri->refreshQuantities();

    computeEnergy(config, logger, logDepth + 1);

    LOG_INFO(logger, logDepth, "Successfully generated a CIT by relocating vertices");
    LOG_INFO(logger, logDepth, "  Energy: {}", config.energy);
    LOG_INFO(logger, logDepth, "createConfiguration END ----------------------------------------");
    LOG_INFO(logger, logDepth, "  elapsed time (ms): {}", timer.milliseconds_str());

    setSolutionVector(config);

    return true;

  } else {
    LOG_INFO(logger, logDepth, "Failed to generate a CIT by relocating vertices");
    LOG_INFO(logger, logDepth, "  elapsed time (ms): {}", timer.milliseconds_str());
    return createConfiguration(id, x, config, logger, logDepth + 1);
  }
}
