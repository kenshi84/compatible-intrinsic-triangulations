#include "cit.hpp"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/utilities/combining_hash_functions.h"
#include<Eigen/SparseCholesky>
#include <cereal/archives/json.hpp>

#include "kt84/graphics/graphics_util.hh"
using namespace kt84::graphics_util;

void cit::readInputMesh(ModelData& mdata, const std::string& inputFilename) {
  std::tie(mdata.mesh, mdata.geometry, mdata.texCoords) = readManifoldSurfaceMesh(inputFilename);

  mdata.nV = mdata.mesh->nVertices();
  DLOG_INFO(0, "Model {}:", mdata.shortName);
  DLOG_INFO(0, "  #vertices: {}", mdata.nV);
  DLOG_INFO(0, "  #faces: {}", mdata.mesh->nFaces());

  // required quantities
  mdata.geometry->requireEdgeLengths();
  mdata.geometry->requireVertexDualAreas();
  mdata.geometry->requireVertexNormals();
  mdata.geometry->requireEdgeCotanWeights();
  mdata.geometry->requireFaceNormals();

  // normalize surface area
  double area = 0;
  for (Face f : mdata.mesh->faces())
    area += mdata.geometry->faceAreas[f];
  DLOG_INFO(0, "  surface area: {} (before normalization)", area);

  for (Vertex v : mdata.mesh->vertices())
    mdata.geometry->inputVertexPositions[v] /= std::sqrt(area);
  mdata.geometry->refreshQuantities();

  // compute bounding box
  for (Vertex v : mdata.mesh->vertices())
    mdata.bbox.extend(mdata.geometry->inputVertexPositions[v].castTo<Vector3d>());

  // setup transport operator
  mdata.transportOperator = EdgeData<Matrix2d>(*mdata.mesh);
  for (Edge e : mdata.mesh->edges()) {
    std::array<Vector2, 4> diamondP = getLayoutDiamond(*mdata.geometry, e.halfedge());
    std::array<Halfedge, 4> diamondHE = e.diamondBoundary();
    std::array<Vertex, 4> diamondV;
    for (int i = 0; i < 4; ++i)
      diamondV[i] = diamondHE[i].vertex();

    Face fThis = e.halfedge().face();
    Face fOther = e.halfedge().twin().face();
    std::array<int, 3> viThis = {0, 1, 2};
    std::array<int, 3> viOther = {2, 3, 0};

    auto readInputMesh_getMatrix = [&diamondP, &diamondV](Face f, const std::array<int, 3>& vi) {
      int offset = 0;
      bool found = false;
      for (; offset < 3; ++offset) {
        if (f.halfedge().vertex() == diamondV[vi[offset]]) {
          found = true;
          break;
        }
      }
      CIT_ASSERT(found);

      Vector2d p0 = diamondP[vi[offset]].castTo<Vector2d>();
      Vector2d p1 = diamondP[vi[(offset + 1) % 3]].castTo<Vector2d>();
      Vector2d p2 = diamondP[vi[(offset + 2) % 3]].castTo<Vector2d>();

      Matrix2d mat;
      mat << p0 - p2, p1 - p2;
      return mat;
    };
    Matrix2d matThis = readInputMesh_getMatrix(fThis, viThis);
    Matrix2d matOther = readInputMesh_getMatrix(fOther, viOther);

    // be careful about the orientation: the transport operator maps barycentric coordinate in this face to the one in the other face
    mdata.transportOperator[e] = matOther.inverse() * matThis;
  }
}

VectorXsp cit::readBarycentricCoordinates(ModelData& mdataFrom, const ModelData& mdataTo, const std::string& filename) {
  DLOG_INFO(0, "Reading barycentric coordinates mapping {} vertices onto {}", mdataFrom.shortName, mdataTo.shortName);

  // Read sparse matrix data from file, store per row
  std::map<int, std::set<std::pair<int, double>>> rows;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open())
    throw std::runtime_error("Couldn't open file " + filename);

  int mode = 0;
  int v1, v2;
  double d;
  while (fin) {
    if (mode == 0) {
      fin >> v1;
      // v1--;
      mode = 1;
    } else if (mode == 1) {
      fin >> v2;
      // v2--;
      mode = 2;
    } else if (mode == 2) {
      fin >> d;
      rows[v1].insert({v2, d});
      mode = 0;
    } else {
      throw std::runtime_error("wrong mode");
    }
  }

  // convert each row data to SurfacePoint
  VectorXsp res(mdataFrom.nV);
  int i = 0;
  for (auto& row : rows) {
    CIT_ASSERT(row.first == i);
    std::array<Vertex, 3> v;
    Vector3 baryCoord;
    int j = 0;
    for (auto& p : row.second) {
      v[j] = mdataTo.mesh->vertex(p.first);
      baryCoord[j] = p.second;
      ++j;
    }

    switch (row.second.size()) {
    case 1:
      // New vertex is being created on an existing vertex --> skip
      DLOG_INFO(0, "vertex point {}, {}", i, row.second.begin()->first);
      CIT_ASSERT(baryCoord[0] == 1.);
      res[i] = SurfacePoint(v[0]);
      if (sysParam.fixAnchors) {
        mdataFrom.anchorVertexIDs.insert(i);
      }
      break;
    case 2: {
      // New vertex is being created on an existing edge
      DLOG_INFO(0, "edge point: {}, {}", v[0].getIndex(), v[1].getIndex());
      CIT_ASSERT(std::fabs(baryCoord[0] + baryCoord[1] - 1) < 0.001);

      Halfedge he = v[0].connectingHalfedge(v[1]);
      if (he == Halfedge())
        throw std::runtime_error("Specified edge could not be found");

      Edge e = he.edge();
      if (e.halfedge() == he)
        res[i] = SurfacePoint(e, baryCoord[1]);
      else
        res[i] = SurfacePoint(e, baryCoord[0]);

      break;
    }
    case 3: {
      // New vertex is being created on an existing face
      // std::DLOG_INFO(0, "face point");
      CIT_ASSERT(std::fabs(baryCoord[0] + baryCoord[1] + baryCoord[2] - 1) < 0.001);

      Face f = v[0].spanningFace(v[1], v[2]);
      if (f == Face())
        throw std::runtime_error("Specified face could not be found");

      res[i] = SurfacePoint(f, {});
      for (int j = 0; j < 3; ++j)
        res[i].faceCoords[vertexIndexInTriangle(v[j], f)] = baryCoord[j];

      break;
    }
    default:
      throw std::runtime_error("Found a row with more than 3 nonzeros");
    }
    ++i;
  }

  return res;
}

std::pair<int, int> cit::getUniqueEdgeID(const CoInTri& cointri, const std::array<Vertex, 2>& vp) {
  int u0 = cointri.uniqueID_per_Vertex[vp[0]];
  int u1 = cointri.uniqueID_per_Vertex[vp[1]];
  if (u1 < u0)
    std::swap(u0, u1);
  return {u0, u1};
}

std::pair<int, int> cit::getUniqueEdgeID(const CoInTri& cointri, Edge e) {
  return getUniqueEdgeID(cointri, e.adjacentVertices());
}

Edge cit::getEdgeByUniqueID(const CoInTri& cointri, const std::pair<int, int>& ue) {
  if (!cointri.uniqueID_to_Vertex.count(ue.first)) return Edge();
  if (!cointri.uniqueID_to_Vertex.count(ue.second)) return Edge();
  Vertex v0 = cointri.uniqueID_to_Vertex.at(ue.first);
  Vertex v1 = cointri.uniqueID_to_Vertex.at(ue.second);
  return v0.connectingEdge(v1);
}

void cit::setCameraMatrix(const ModelData& mdata) {
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

std::array<Vector2, 4> cit::getLayoutDiamond(const IntrinsicGeometryInterface& geometry, Halfedge iHe) {
  // Blatant copy of SignpostIntrinsicTriangulation::layoutDiamond (due to being protected)

  // Conventions:
  //  - iHe points from vertex 2 to vertex 0, other vertices are numbered ccw
  //  - iHe is incident on face A, other is face B
  //  - halfedges within face are numbered CCW as A0, A1, A2 (etc),
  //    starting with iHe and twin(iHe)
  //  - When we lay out the triangle, p3 is at the origin and
  //    edge 3-0 is along the X-axis
  //  - flips is always ccw, so iHe points from vertex 3 --> 1 after

  // Gather index values
  Halfedge iHeA0 = iHe;
  Halfedge iHeA1 = iHeA0.next();
  Halfedge iHeA2 = iHeA1.next();
  Halfedge iHeB0 = iHe.twin();
  Halfedge iHeB1 = iHeB0.next();
  Halfedge iHeB2 = iHeB1.next();

  // Gather length values
  double l01 = geometry.edgeLengths[iHeA1.edge()];
  double l12 = geometry.edgeLengths[iHeA2.edge()];
  double l23 = geometry.edgeLengths[iHeB1.edge()];
  double l30 = geometry.edgeLengths[iHeB2.edge()];
  double l02 = geometry.edgeLengths[iHeA0.edge()];

  // Lay out the vertices of the diamond
  Vector2 p3{0., 0.};
  Vector2 p0{l30, 0.};
  Vector2 p2 = layoutTriangleVertexFromLength(p3, p0, l02, l23); // involves more arithmetic than strictly necessary
  Vector2 p1 = layoutTriangleVertexFromLength(p2, p0, l01, l12);

  return {p0, p1, p2, p3};
}

std::unordered_map<Vertex, Vector2> cit::layoutBaseFaces(const CoInTri& cointri, Face intrinsic_f) {
  std::unordered_map<Vertex, Vector2> layoutVertexPosition;

  // collect all input faces to be laid out
  std::set<Face> input_faces;
  std::set<Edge> input_internal_edges;

  std::array<SurfacePoint, 3> sp = {
    cointri.signpostTri->vertexLocations[intrinsic_f.halfedge().vertex()],
    cointri.signpostTri->vertexLocations[intrinsic_f.halfedge().next().vertex()],
    cointri.signpostTri->vertexLocations[intrinsic_f.halfedge().next().next().vertex()]
  };

  for (int i = 0; i < 3; ++i) {
    // For each intrinsic vertex, if its location is a face point, add that face to input_faces
    if (sp[i].type == SurfacePointType::Face) {
      input_faces.insert(sp[i].face);

    // If edge point, add its two adjacent faces to input_faces, and its edge to input_internal_edges
    } else if (sp[i].type == SurfacePointType::Edge) {
      input_faces.insert(sp[i].edge.halfedge().face());
      input_faces.insert(sp[i].edge.halfedge().twin().face());
      input_internal_edges.insert(sp[i].edge);
    }
  }

  // Add input faces touched by intrinsic halfedges
  for (Halfedge intrinsic_he : intrinsic_f.adjacentHalfedges()) {
    // If this intrinsic edge is partially original, add the relevant face to input_faces
    Edge intrinsic_e = intrinsic_he.edge();
    Edge input_e;
    double tEdgeMin, tEdgeMax;
    bool reversed;
    if (cointri.signpostTri->isIntrinsicEdgePartiallyOriginal(intrinsic_e, &input_e, &tEdgeMin, &tEdgeMax, &reversed)) {
      Halfedge input_he = input_e.halfedge();
      if (intrinsic_e.halfedge() == intrinsic_he) {
        if (reversed)
          input_he = input_he.twin();
      } else {
        if (!reversed)
          input_he = input_he.twin();
      }
      input_faces.insert(input_he.face());

    // Otherwise, visit each edge point in this intrinsic edge's path, and add two faces adjacent to input edge
    } else {
      for (size_t i = 1; i < cointri.intrinsicEdgePath[intrinsic_e].size() - 1; ++i) {
        const SurfacePoint& p = cointri.intrinsicEdgePath[intrinsic_e][i];
        CIT_ASSERT(p.type == SurfacePointType::Edge);
        if (std::fabs(p.tEdge) < 1.0e-5 || std::fabs(1. - p.tEdge) < 1.0e-5) {
          // ignore garbage erroneous intersections at endpoints
          continue;
        }
        Edge input_e = p.edge;
        input_faces.insert(input_e.halfedge().face());
        input_faces.insert(input_e.halfedge().twin().face());
        input_internal_edges.insert(input_e);
      }
    }
  }

  CIT_ASSERT(!input_faces.empty());

  // ensure input_faces has no internal vertex
  auto isVertexOnBoundary = [&input_faces, &input_internal_edges](Vertex v) {
    Halfedge he = v.halfedge();
    while (true) {
      if (!input_faces.count(he.face()))
        return true;
      if (!input_internal_edges.count(he.edge()))
        return true;
      he = he.twin().next();
      if (he == v.halfedge())
        break;
    }
    return false;
  };
  for (Face f : input_faces) {
    for (Vertex v : f.adjacentVertices())
      CIT_ASSERT(isVertexOnBoundary(v));
  }

  // initial layout for the first face in input_faces
  Halfedge he0 = input_faces.begin()->halfedge();
  Halfedge he1 = he0.next();
  Halfedge he2 = he1.next();
  layoutVertexPosition[he0.vertex()] = {0,0};
  layoutVertexPosition[he1.vertex()] = cointri.signpostTri->inputGeom.halfedgeVectorsInFace[he0];
  layoutVertexPosition[he2.vertex()] = -cointri.signpostTri->inputGeom.halfedgeVectorsInFace[he2];

  // lay out the rest by floodfill
  size_t floodFillCount = 1;
  std::deque<Halfedge> floodFillFront;
  if (input_internal_edges.count(he0.edge())) floodFillFront.push_back(he0);
  if (input_internal_edges.count(he1.edge())) floodFillFront.push_back(he1);
  if (input_internal_edges.count(he2.edge())) floodFillFront.push_back(he2);

  while (!floodFillFront.empty()) {
    Halfedge he = floodFillFront.front();
    floodFillFront.pop_front();
    Halfedge heAB = he.twin();
    CIT_ASSERT(input_faces.count(heAB.face()));
    ++floodFillCount;
    Halfedge heBC = heAB.next();
    Halfedge heCA = heBC.next();
    Vertex vA = heAB.vertex();
    Vertex vB = heBC.vertex();
    Vertex vC = heCA.vertex();
    layoutVertexPosition[vC] = layoutTriangleVertexFromLength(layoutVertexPosition.at(vA),
                                                              layoutVertexPosition.at(vB),
                                                              cointri.signpostTri->inputGeom.edgeLengths[heBC.edge()],
                                                              cointri.signpostTri->inputGeom.edgeLengths[heCA.edge()]);

    if (input_internal_edges.count(heBC.edge())) floodFillFront.push_back(heBC);
    if (input_internal_edges.count(heCA.edge())) floodFillFront.push_back(heCA);
  }
  CIT_ASSERT(floodFillCount == input_faces.size());

  return layoutVertexPosition;
}

namespace cit { namespace detail {

ADVector2d_6v_2nd getTriangleCornerAD(const std::unordered_map<Vertex, Vector2>& layoutVertexPosition, SurfacePoint sp, size_t varIndex = (size_t)-1) {
  if (sp.type == SurfacePointType::Vertex) {
    CIT_ASSERT(varIndex == (size_t)-1);

    return layoutVertexPosition.at(sp.vertex).castTo<ADVector2d_6v_2nd>();

  } else {
    CIT_ASSERT(varIndex != (size_t)-1);

    // If edge point, reinterpret it as a face point
    if (sp.type == SurfacePointType::Edge)
      sp = convertEdgePointToFacePoint(sp);

    // Get the 2x3 matrix holding triangle's corner positions
    Eigen::Matrix<double, 2, 3> M;
    M <<
      layoutVertexPosition.at(sp.face.halfedge().vertex()).castTo<Vector2d>(),
      layoutVertexPosition.at(sp.face.halfedge().next().vertex()).castTo<Vector2d>(),
      layoutVertexPosition.at(sp.face.halfedge().next().next().vertex()).castTo<Vector2d>();

    // Get the 3D AutoDiff variable vector for this face coords
    ADVector3d_6v_2nd bc = sp.faceCoords.castTo<ADVector3d_6v_2nd>();
    for (size_t j = 0; j < 2; ++j) {
      bc[j].value().derivatives() = Vector6d::Unit(6, 2 * varIndex + j);
      bc[j].derivatives() = Vector6d::Unit(6, 2 * varIndex + j);
    }
    bc[2] = 1. - bc[0] - bc[1];

    // Get the output 2D AutoDiff variable vector by matrix-vector product
    return M.cast<ADScalar_6v_2nd>() * bc;
  }
}

} }

ADScalar_6v_2nd cit::computeTriangleEnergyAD(Configuration& config, Face A_f) {
  Face B_f = config.cointriA.correspondingFace[A_f];

  const std::unordered_map<Vertex, Vector2> A_layoutVertexPosition = layoutBaseFaces(config.cointriA, A_f);
  const std::unordered_map<Vertex, Vector2> B_layoutVertexPosition = layoutBaseFaces(config.cointriB, B_f);

  Halfedge A_he = A_f.halfedge();
  Halfedge B_he = config.cointriA.correspondingHalfedge[A_he];

  std::array<ADVector2d_6v_2nd, 3> A_triangle;
  std::array<ADVector2d_6v_2nd, 3> B_triangle;

  for (int i = 0; i < 3; ++i) {
    Vertex A_v = A_he.vertex();
    Vertex B_v = B_he.vertex();

    A_he = A_he.next();
    B_he = B_he.next();

    SurfacePoint A_sp = config.cointriA.signpostTri->vertexLocations[A_v];
    SurfacePoint B_sp = config.cointriB.signpostTri->vertexLocations[B_v];

    if (A_sp.type == SurfacePointType::Vertex && B_sp.type == SurfacePointType::Vertex) {
      // Both A and B are fixed (anchor vertex)
      CIT_ASSERT(mdataA.anchorVertexIDs.count(A_v.getIndex()));
      CIT_ASSERT(mdataB.anchorVertexIDs.count(B_v.getIndex()));

      A_triangle[i] = detail::getTriangleCornerAD(A_layoutVertexPosition, A_sp);
      B_triangle[i] = detail::getTriangleCornerAD(B_layoutVertexPosition, B_sp);

    } else if (A_sp.type == SurfacePointType::Vertex) {
      // A is fixed, B is variable
      A_triangle[i] = detail::getTriangleCornerAD(A_layoutVertexPosition, A_sp);
      B_triangle[i] = detail::getTriangleCornerAD(B_layoutVertexPosition, B_sp, i);

    } else {
      // A is variable, B is fixed
      CIT_ASSERT(B_sp.type == SurfacePointType::Vertex);
      A_triangle[i] = detail::getTriangleCornerAD(A_layoutVertexPosition, A_sp, i);
      B_triangle[i] = detail::getTriangleCornerAD(B_layoutVertexPosition, B_sp);
    }
  }

  return computeTriangleEnergy(A_triangle, B_triangle);
}

void cit::setSolutionVector(Configuration& config) {
  auto setSolutionVector_inner = [&] (const CoInTri& cointriP, const CoInTri& cointriQ) -> VectorXsp {
    VectorXsp res(cointriP.mdata->nV);
    for (size_t i = 0; i < cointriP.mdata->nV; ++i) {
      Vertex P_v = cointriP.signpostTri->intrinsicMesh->vertex(i);
      Vertex Q_v = cointriP.correspondingVertex[P_v];
      res[i] = cointriQ.signpostTri->vertexLocations[Q_v];
    }
    return res;
  };
  config.x.resize(mdataA.nV + mdataB.nV);
  config.x.head(mdataA.nV) = setSolutionVector_inner(config.cointriA, config.cointriB);
  config.x.tail(mdataB.nV) = setSolutionVector_inner(config.cointriB, config.cointriA);
}

bool cit::isEdgeFlippable(const CoInTri& cointri, Edge e) {
  if (!cointri.signpostTri->flipEdgeIfPossible(e, 1e-6, true)) return false;
  if (e.halfedge().vertex().degree() < 4) return false;
  if (e.halfedge().twin().vertex().degree() < 4) return false;
  // Check if there arleady exists another edge connecting the same pair, and reject if yes
  auto v = getOppositeVertices(e);
  if (v[0].connectingEdge(v[1]) != Edge())
    return false;
  return true;
}

void cit::checkVertexDegree(const CoInTri& cointri) {
  for (Vertex v : cointri.signpostTri->intrinsicMesh->vertices())
    CIT_ASSERT(v.degree() >= 3);
}

void cit::writeBlobToFile(const std::string& filename, const std::string& blob) {
  std::ofstream ofs(filename.c_str(), std::ios::binary);
  if (!ofs.is_open())
    throw std::runtime_error("Couldn't open file: " + filename);

  const uint64_t sz = blob.size();
  ofs.write((const char*)&sz, sizeof(uint64_t));
  ofs.write(&blob[0], sz);
}

std::string cit::readBlobFromFile(const std::string& filename) {
  if (!filesystem::exists(filename))
    throw std::runtime_error("File doesn't exist: " + filename);

  if (!filesystem::is_regular_file(filename))
    throw std::runtime_error("Not a regular file: " + filename);

  std::ifstream ifs(filename.c_str(), std::ios::binary);
  if (!ifs.is_open())
    throw std::runtime_error("Couldn't open file: " + filename);

  uint64_t sz;
  ifs.read((char*)&sz, sizeof(uint64_t));
  if (ifs.gcount() != sizeof(uint64_t))
    throw std::runtime_error("Failed to read the blob size; read bytes: " + std::to_string(ifs.gcount()));

  std::string blob(sz, '\0');
  ifs.read(&blob[0], sz);
  if ((uint64_t)ifs.gcount() != sz)
    throw std::runtime_error("Failed to read expected number of bytes for a blob: expected=" + std::to_string(sz) + ", read=" + std::to_string(ifs.gcount()));

  return blob;
}

std::vector<cit::PatchPair> cit::extractInconsistentPatches(const Configuration& config, bool useVertexAdjacency) {
  std::vector<PatchPair> patchPairs;

  const CoInTri& cointriA = config.cointriA;
  const CoInTri& cointriB = config.cointriB;

  FaceData<char> A_visited{*cointriA.signpostTri->intrinsicMesh, false};
  FaceData<char> B_visited{*cointriB.signpostTri->intrinsicMesh, false};
  for (Edge A_e : cointriA.signpostTri->intrinsicMesh->edges()) {
    if (!isEdgeCompatible(cointriA, A_e))
      continue;

    Face A_f0 = A_e.halfedge().face();
    Face A_f1 = A_e.halfedge().twin().face();
    if (isFaceCompatible(cointriA, A_f0) + isFaceCompatible(cointriA, A_f1) != 1)
      continue;

    if (isFaceCompatible(cointriA, A_f0))
      std::swap(A_f0, A_f1);

    if (A_visited[A_f0])
      continue;

    Edge B_e = cointriA.correspondingEdge[A_e];
    Face B_f0 = B_e.halfedge().face();
    Face B_f1 = B_e.halfedge().twin().face();
    if (isFaceCompatible(cointriB, B_f0))
      std::swap(B_f0, B_f1);
    CIT_ASSERT(!isFaceCompatible(cointriB, B_f0));
    // CIT_ASSERT(!B_visited[B_f0]);
    if (B_visited[B_f0]) continue; // TODO: figure out why this assertion is wrong

    auto extractInconsistentPatches_inner = [useVertexAdjacency](const CoInTri& cointri, Face f0, FaceData<char>& visited) -> std::set<Face> {
      std::set<Face> faces;
      std::deque<Face> facesToProcess;
      facesToProcess.push_back(f0);
      while (!facesToProcess.empty()) {
        Face f1 = facesToProcess.front();
        facesToProcess.pop_front();
        CIT_ASSERT(!visited[f1]);
        visited[f1] = true;
        faces.insert(f1);
        if (useVertexAdjacency) {
          for (Vertex fv : f1.adjacentVertices()) {
            for (Face f2 : fv.adjacentFaces()) {
              if (std::find(facesToProcess.begin(), facesToProcess.end(), f2) != facesToProcess.end())
                continue;
              if (!isFaceCompatible(cointri, f2) && !visited[f2])
                facesToProcess.push_back(f2);
            }
          }
        } else {
          for (Face f2 : f1.adjacentFaces()) {
            if (std::find(facesToProcess.begin(), facesToProcess.end(), f2) != facesToProcess.end())
              continue;
            if (!isFaceCompatible(cointri, f2) && !visited[f2])
              facesToProcess.push_back(f2);
          }
        }
      }
      return faces;
    };
    Patch A_patch, B_patch;
    A_patch.faces = extractInconsistentPatches_inner(cointriA, A_f0, A_visited);
    B_patch.faces = extractInconsistentPatches_inner(cointriB, B_f0, B_visited);
    patchPairs.push_back({A_patch, B_patch});
  }
  return patchPairs;
}

void cit::analyzeInconsistentPatch(const CoInTri& cointri, Patch& patch) {
  patch.edges.clear();
  patch.vertices.clear();
  patch.interiorEdges.clear();
  patch.interiorVertices.clear();
  patch.boundaryHalfedges.clear();
  patch.vertexAngleSums.clear();

  // Identify set of edges & vertices
  for (Face f : patch.faces) {
    for (Halfedge he : f.adjacentHalfedges()) {
      patch.edges.insert(he.edge());
      patch.vertices.insert(he.vertex());
    }
  }

  // Identify set of boundary/interior edges
  for (Edge e : patch.edges) {
    bool isBoundary = false;
    for (Face f : e.adjacentFaces()) {
      if (!patch.faces.count(f)) {
        isBoundary = true;
        break;
      }
    }
    if (isBoundary)
      patch.boundaryEdges.insert(e);
    else
      patch.interiorEdges.insert(e);
  }

  // Identify set of boundary/interior vertices
  for (Vertex v : patch.vertices) {
    bool isBoundary = false;
    for (Face f : v.adjacentFaces()) {
      if (!patch.faces.count(f)) {
        isBoundary = true;
        break;
      }
    }
    if (isBoundary)
      patch.boundaryVertices.insert(v);
    else
      patch.interiorVertices.insert(v);
  }

  // Inequality holds when any of the boundary vertex's topology is non-half-disk
  CIT_ASSERT(patch.boundaryEdges.size() >= patch.boundaryVertices.size());

  // Identify ordered list of boundary halfedges
  std::set<Edge> visited;
  size_t count = 0;
  while (visited.size() < patch.boundaryEdges.size()) {
    Halfedge heStart;
    for (Edge e : patch.boundaryEdges) {
      if (!visited.count(e)) {
        heStart = e.halfedge();
        if (!patch.faces.count(heStart.face()))
          heStart = heStart.twin();
        CIT_ASSERT(patch.faces.count(heStart.face()));
        break;
      }
    }
    CIT_ASSERT(heStart.getMesh());
    patch.boundaryHalfedges.push_back({});
    for (Halfedge he = heStart; ; ) {
      visited.insert(he.edge());
      patch.boundaryHalfedges.back().push_back(he);

      he = he.twin();
      while (!patch.faces.count(he.face()))
        he = he.next().next().twin();

      if (he == heStart)
        break;
    }
    count += patch.boundaryHalfedges.back().size();
  }
  CIT_ASSERT(count == patch.boundaryEdges.size());

  // Set normalized vertex angle sums
  for (Vertex v : patch.vertices) {
    patch.vertexAngleSums[v] = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (patch.faces.count(he.face()))
        patch.vertexAngleSums[v] += cointri.signpostTri->cornerAngles[he.corner()];
    }
    patch.vertexAngleSums[v] *= 2. * M_PI / cointri.signpostTri->intrinsicVertexAngleSums[v];
  }
}

double cit::getMinAngle(const CoInTri& cointri) {
  kt84::MinSelector<Halfedge> minAngleHE;
  for (Halfedge he : cointri.signpostTri->intrinsicMesh->halfedges())
    minAngleHE.update(cointri.signpostTri->cornerAngle(he.corner()), he);
  return minAngleHE.score;
}

void cit::fillResultInfo(Configuration& newConfig, double currentConfigEnergy, const VectorXd& previousGradient, const VectorXd& previousDescentDirection, int l, double stepSize) {
  ResultInfo& ri = resultInfoTable[newConfig.id];

  ri.nIncompatibleEdges = newConfig.cointriA.signpostTri->intrinsicMesh->nEdges() - newConfig.nCompatibleEdges;
  ri.nIncompatibleFaces = newConfig.cointriA.signpostTri->intrinsicMesh->nFaces() - newConfig.nCompatibleFaces;

  ri.energy = computeEnergy(newConfig);
  ri.energyDelta = newConfig.energy - currentConfigEnergy;

  ri.A_minAngle = getMinAngle(newConfig.cointriA);
  ri.B_minAngle = getMinAngle(newConfig.cointriB);

  ri.previous_wL = current_wL;
  ri.previous_l = l;
  if (previousGradient.size()) {
    ri.previous_gradientNorm = previousGradient.norm();
    ri.previous_descentDirectionNorm = previousDescentDirection.norm();
    ri.previous_slope = previousGradient.dot(previousDescentDirection);
  }
  ri.stepSize = stepSize;
}

void cit::writeResult(const Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, double smax, size_t nSteps, ResultContext context) {
  const ResultInfo& ri = resultInfoTable[config.id];

  std::string filename_base = getResultFileName(context, config.id, config.previous_id, ri.previous_l);

  // Serialize
  if (context == ResultContext::UpdateConfiguration) {
    DLOG_INFO(0, "");
    DLOG_INFO(0, "");
    DLOG_INFO(0, "++++ Serializing to {}.dat", filename_base);
    writeBlobToFile(filename_base + ".dat", toSerializedBlob(config, history, wL, smax, sysParam.firstOrderMode, nSteps));
  }

  // Write meta info
  std::ofstream fout_json((filename_base + ".json").c_str());
  cereal::JSONOutputArchive ar(fout_json);
  ar(cereal::make_nvp("nIncompatibleEdges", ri.nIncompatibleEdges));
  ar(cereal::make_nvp("nIncompatibleFaces", ri.nIncompatibleFaces));
  ar(cereal::make_nvp("energy", ri.energy));
  ar(cereal::make_nvp("energyDelta", ri.energyDelta));
  ar(cereal::make_nvp("A_minAngle", ri.A_minAngle));
  ar(cereal::make_nvp("B_minAngle", ri.B_minAngle));
  ar(cereal::make_nvp("A_maxSignpostError", ri.A_maxSignpostError));
  ar(cereal::make_nvp("B_maxSignpostError", ri.B_maxSignpostError));
  ar(cereal::make_nvp("previous_id", config.previous_id));
  ar(cereal::make_nvp("previous_wL", ri.previous_wL));
  ar(cereal::make_nvp("previous_l", ri.previous_l));
  ar(cereal::make_nvp("previous_gradientNorm", ri.previous_gradientNorm));
  ar(cereal::make_nvp("previous_descentDirectionNorm", ri.previous_descentDirectionNorm));
  ar(cereal::make_nvp("previous_slope", ri.previous_slope));
  ar(cereal::make_nvp("stepSize", ri.stepSize));
  ar(cereal::make_nvp("firstOrderMode", sysParam.firstOrderMode));
  ar(cereal::make_nvp("nSteps", nSteps));

  DLOG_INFO(0, "++++ Meta info written to {}.json", filename_base);
}

uint64_t cit::makeConfigurationID() {
  static std::mutex mtx;
  static std::unordered_set<uint64_t> existingIDs;
  std::lock_guard<std::mutex> lock(mtx);
  uint64_t id = (uint64_t)(std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1));
  while (existingIDs.count(id))
    ++id;
  existingIDs.insert(id);
  return id;
}

namespace cit { namespace detail {
bool isSegmentOnEdge(SurfacePoint sp0, SurfacePoint sp1) {
  CIT_ASSERT(sp0.getMesh() == sp1.getMesh());

  if (sp0.type == SurfacePointType::Face) return false;
  if (sp1.type == SurfacePointType::Face) return false;

  if (sp0.type == SurfacePointType::Vertex && sp1.type == SurfacePointType::Vertex)
    return sp0.vertex.connectingHalfedge(sp1.vertex) != Halfedge();

  if (sp0.type == SurfacePointType::Edge && sp1.type == SurfacePointType::Edge)
    return sp0.edge == sp1.edge;

  if (sp0.type != SurfacePointType::Edge)
    std::swap(sp0, sp1);

  return sp0.edge.halfedge().vertex() == sp1.vertex || sp0.edge.halfedge().tipVertex() == sp1.vertex;
}
}}

void cit::writeOverlayMesh(const Configuration& config) {
  DLOG_TRACE(0, "BEGIN writeOverlayMesh");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(0, "END   writeOverlayMesh"); });

  const size_t nV = config.overlayWedgesPerOverlayVertex.size();
  const size_t nF = config.overlayPolygons.size();

  std::unordered_map<OverlayVertex, size_t> vid_per_overlayVertex;
  std::unordered_map<size_t, OverlayVertex> vid_to_overlayVertex;
  for (const auto& p : config.overlayWedgesPerOverlayVertex) {
    OverlayVertex v = p.first;
    size_t vid = vid_per_overlayVertex.size();
    vid_per_overlayVertex[v] = vid;
    vid_to_overlayVertex[vid] = v;
  }
  CIT_ASSERT(vid_per_overlayVertex.size() == nV && vid_to_overlayVertex.size() == nV);

#if 0
  std::vector<std::vector<size_t>> polygons(nF);
  size_t fid = 0;
  for (const OverlayPolygon& op : config.overlayPolygons) {
    for (const OverlayWedge& wedge : op.wedges)
      polygons[fid].push_back(vid_per_overlayVertex.at(wedge.overlayVertex));;
    ++fid;
  }
  CIT_ASSERT(fid == nF);

  ManifoldSurfaceMesh overlayMesh(polygons);
  VertexPositionGeometry A_omGeometry(overlayMesh);
  VertexPositionGeometry B_omGeometry(overlayMesh);
  CornerData<Vector2> A_texCoords(overlayMesh);
  CornerData<Vector2> B_texCoords(overlayMesh);

  // Fill position data
  for (const auto& p : config.overlayWedgesPerOverlayVertex) {
    size_t vid = vid_per_overlayVertex.at(p.first);
    Vertex v = overlayMesh.vertex(vid);

    const OverlayWedge& wedge = p.second.front();
    Vector3 A_pos = wedge.facePointA.interpolate(mdataA.geometry->inputVertexPositions);
    Vector3 B_pos = wedge.facePointB.interpolate(mdataB.geometry->inputVertexPositions);

    A_omGeometry.inputVertexPositions[v] = A_pos;
    B_omGeometry.inputVertexPositions[v] = B_pos;
  }

  // Fill texCoords
  for (const OverlayPolygon& op : config.overlayPolygons) {
    const size_t nW = op.wedges.size();
    for (size_t j0 = 0; j0 < nW; ++j0) {
      size_t j1 = (j0 + 1) % nW;

      // Identify corner
      const OverlayWedge& wedge0 = op.wedges[j0];
      const OverlayWedge& wedge1 = op.wedges[j1];

      size_t vid0 = vid_per_overlayVertex.at(wedge0.overlayVertex);
      size_t vid1 = vid_per_overlayVertex.at(wedge1.overlayVertex);

      Vertex v0 = overlayMesh.vertex(vid0);
      Vertex v1 = overlayMesh.vertex(vid1);

      Halfedge h = v0.connectingHalfedge(v1);
      CIT_ASSERT(h != Halfedge());

      Corner c = h.corner();

      // Get interpolated texture coordinate at the wedge
      A_texCoords[c] = {0., 0.};
      B_texCoords[c] = {0., 0.};

      SurfacePoint facePointA = wedge0.facePointA;
      SurfacePoint facePointB = wedge0.facePointB;
      CIT_ASSERT(facePointA.type == SurfacePointType::Face);
      CIT_ASSERT(facePointB.type == SurfacePointType::Face);

      Halfedge A_he = facePointA.face.halfedge();
      Halfedge B_he = facePointB.face.halfedge();
      for (int i = 0; i < 3; ++i) {
        if (mdataA.texCoords) A_texCoords[c] += facePointA.faceCoords[i] * (*mdataA.texCoords)[A_he.corner()];
        if (mdataB.texCoords) B_texCoords[c] += facePointB.faceCoords[i] * (*mdataB.texCoords)[B_he.corner()];
        A_he = A_he.next();
        B_he = B_he.next();
      }
    }
  }
#endif

  filesystem::create_directory(getResultFileNameBase());
  filesystem::create_directory(getResultFileNameBase() + "/meshes");

  std::string filename_AonB = getResultFileNameBase() + "/meshes/" + std::to_string(config.id) + "_AonB";
  std::string filename_BonA = getResultFileNameBase() + "/meshes/" + std::to_string(config.id) + "_BonA";

#if 0
  if (mdataB.texCoords) writeSurfaceMesh(overlayMesh, A_omGeometry, B_texCoords, filename_BonA + ".obj"); else writeSurfaceMesh(overlayMesh, A_omGeometry, filename_BonA + ".obj");
  if (mdataA.texCoords) writeSurfaceMesh(overlayMesh, B_omGeometry, A_texCoords, filename_AonB + ".obj"); else writeSurfaceMesh(overlayMesh, B_omGeometry, filename_AonB + ".obj");

  DLOG_INFO(0, "+++++ Overlay mesh written to {}", filename_AonB + ".obj");
  DLOG_INFO(0, "+++++ Overlay mesh written to {}", filename_BonA + ".obj");
#endif

  //--------------------------------+
  // Write edge images as ply files |
  //--------------------------------+
  // Collect edges
  std::vector<std::array<size_t, 2>> A_edges;
  std::vector<std::array<size_t, 2>> B_edges;
  for (const OverlayPolygon& op : config.overlayPolygons) {
    const size_t nW = op.wedges.size();
    for (size_t j0 = 0; j0 < nW; ++j0) {
      size_t j1 = (j0 + 1) % nW;

      const OverlayWedge& wedge0 = op.wedges[j0];
      const OverlayWedge& wedge1 = op.wedges[j1];

      const OverlayVertex& v0 = wedge0.overlayVertex;
      const OverlayVertex& v1 = wedge1.overlayVertex;

      size_t vid0 = vid_per_overlayVertex.at(v0);
      size_t vid1 = vid_per_overlayVertex.at(v1);
      if (vid0 > vid1)
        continue;

      SurfacePoint A_sp0 = wedge0.facePointA.reduced();
      SurfacePoint A_sp1 = wedge1.facePointA.reduced();
      SurfacePoint B_sp0 = wedge0.facePointB.reduced();
      SurfacePoint B_sp1 = wedge1.facePointB.reduced();

      if (detail::isSegmentOnEdge(A_sp0, A_sp1)) A_edges.push_back({vid0, vid1});
      if (detail::isSegmentOnEdge(B_sp0, B_sp1)) B_edges.push_back({vid0, vid1});
    }
  }

  // How to draw lines in .ply files using Meshlab?
  // https://stackoverflow.com/a/58648176
  /*
    ply
    format ascii 1.0
    element vertex 2
    property float x
    property float y
    property float z
    element edge 1                        
    property int vertex1                  
    property int vertex2                  
    end_header
    0 0 0 
    0 0 1 
    0 1
  */

  std::ostringstream oss_header_pre;
  oss_header_pre << "ply" << std::endl;
  oss_header_pre << "format ascii 1.0" << std::endl;
  oss_header_pre << "element vertex " << nV << std::endl;
  oss_header_pre << "property float x" << std::endl;
  oss_header_pre << "property float y" << std::endl;
  oss_header_pre << "property float z" << std::endl;
  oss_header_pre << "element edge ";

  std::ostringstream oss_header_post;
  oss_header_post << std::endl;
  oss_header_post << "property int vertex1" << std::endl;
  oss_header_post << "property int vertex2" << std::endl;
  oss_header_post << "end_header" << std::endl;

  std::ofstream fout_AonB((filename_AonB + ".ply").c_str());
  std::ofstream fout_BonA((filename_BonA + ".ply").c_str());

  fout_AonB << oss_header_pre.str() << A_edges.size() << oss_header_post.str();
  fout_BonA << oss_header_pre.str() << B_edges.size() << oss_header_post.str();

  // Write vertices
  for (size_t vid = 0; vid < nV; ++vid) {
    OverlayVertex v = vid_to_overlayVertex.at(vid);
    OverlayWedge wedge = config.overlayWedgesPerOverlayVertex.at(v).front();

    Vector3 A_pos = wedge.facePointA.interpolate(mdataA.geometry->inputVertexPositions);
    Vector3 B_pos = wedge.facePointB.interpolate(mdataB.geometry->inputVertexPositions);

    fout_BonA << A_pos.x << " " << A_pos.y << " " << A_pos.z << std::endl;
    fout_AonB << B_pos.x << " " << B_pos.y << " " << B_pos.z << std::endl;
  }

  // Write edges
  for (std::array<size_t, 2> e : A_edges) fout_AonB << e[0] << " " << e[1] << std::endl;
  for (std::array<size_t, 2> e : B_edges) fout_BonA << e[0] << " " << e[1] << std::endl;

  DLOG_INFO(0, "+++++ Overlay mesh edges written to {}", filename_AonB + ".ply");
  DLOG_INFO(0, "+++++ Overlay mesh edges written to {}", filename_BonA + ".ply");
}

spdlog::logger_ptr cit::makeLogger(const std::string& filename, bool alongWithDefault, std::ostringstream* additional_oss) {
  std::vector<spdlog::sink_ptr> sinks = { std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, true) };

  if (alongWithDefault) {
    for (auto& s : spdlog::default_logger()->sinks())
      sinks.push_back(s);
  }

  if (additional_oss) {
    sinks.push_back(std::make_shared<spdlog::sinks::ostream_sink_mt>(*additional_oss));
  }

  spdlog::logger_ptr logger = std::make_shared<spdlog::logger>("", sinks.begin(), sinks.end());
  spdlog::initialize_logger(logger);
  return logger;
}

std::vector<Halfedge> cit::findFaceStrip(Face fFrom, Face fTo) {
  SurfaceMesh* mesh = fFrom.getMesh();
  CIT_ASSERT(mesh && mesh == fTo.getMesh());

  if (fFrom == fTo)
    return {};

  std::deque<Halfedge> halfedgesToProcess = {
    fFrom.halfedge(),
    fFrom.halfedge().next(),
    fFrom.halfedge().next().next()
  };
  FaceData<char> visited(*mesh, false);
  visited[fFrom] = true;
  std::map<Halfedge, Halfedge> prev;
  Halfedge heFound;

  while (true) {
    Halfedge hePrev = halfedgesToProcess.front();
    halfedgesToProcess.pop_front();
    Face fCurr = hePrev.twin().face();
    if (visited[fCurr])
      continue;
    visited[fCurr] = true;

    if (fCurr == fTo) {
      heFound = hePrev;
      break;
    }

    for (Halfedge heNext : fCurr.adjacentHalfedges()) {
      if (!visited[heNext.twin().face()]) {
        halfedgesToProcess.push_back(heNext);
        prev[heNext] = hePrev;
      }
    }
  }

  std::vector<Halfedge> res;
  Halfedge he = heFound;
  while (true) {
    res.push_back(he);
    if (he.face() == fFrom)
      break;
    he = prev[he];
  }

  std::reverse(res.begin(), res.end());
  return res;
}

bool cit::isPatchSizeConsistent(const Patch& p1, const Patch& p2) {
  if (p1.boundaryHalfedges.size() != p2.boundaryHalfedges.size()) return false;
  for (size_t i = 0; i < p1.boundaryHalfedges.size(); ++i) {
    if (p1.boundaryHalfedges[i].size() != p2.boundaryHalfedges[i].size()) return false;
  }
  return
    p1.faces.size() == p2.faces.size() &&
    p1.edges.size() == p2.edges.size() &&
    p1.vertices.size() == p2.vertices.size() &&
    p1.interiorEdges.size() == p2.interiorEdges.size() &&
    p1.interiorVertices.size() == p2.interiorVertices.size();
}

void cit::makeIntrinsicEdgesConsistentlyOriented(CoInTri& cointriP, const CoInTri& cointriQ) {
  for (Edge P_eIntrinsic : cointriP.signpostTri->intrinsicMesh->edges()) {
    Edge Q_eIntrinsic = cointriP.correspondingEdge[P_eIntrinsic];

    // Figure out if P_eIntrinsic and Q_eIntrinsic are in the same orientation or not
    std::array<Vertex, 2> P_vIntrinsic = P_eIntrinsic.adjacentVertices();
    std::array<Vertex, 2> Q_vIntrinsic = Q_eIntrinsic.adjacentVertices();
    for (int i = 0; i < 2; ++i) {
      Q_vIntrinsic[i] = cointriQ.correspondingVertex[Q_vIntrinsic[i]];   // Slight abuse of prefix Q_
    }
    if (P_vIntrinsic != Q_vIntrinsic) {
      CIT_ASSERT(P_vIntrinsic[0] == Q_vIntrinsic[1]);
      CIT_ASSERT(P_vIntrinsic[1] == Q_vIntrinsic[0]);
    }

    // Dummy variables for isIntrinsicEdgePartiallyOriginal()
    Edge eInput;
    double tEdgeMin, tEdgeMax;

    // Check if this intrinsic edge is original on P or Q
    bool P_reversed;
    bool Q_reversed;
    bool P_original = cointriP.signpostTri->isIntrinsicEdgePartiallyOriginal(P_eIntrinsic, &eInput, &tEdgeMin, &tEdgeMax, &P_reversed);
    bool Q_original = cointriQ.signpostTri->isIntrinsicEdgePartiallyOriginal(Q_eIntrinsic, &eInput, &tEdgeMin, &tEdgeMax, &Q_reversed);

    // The sense of 'reversed' in Q is dependent on the relative orientation of P_eIntrinsic & Q_eIntrinsic
    if (Q_original)
      Q_reversed = (Q_reversed && P_vIntrinsic == Q_vIntrinsic) || (!Q_reversed && P_vIntrinsic != Q_vIntrinsic);

    bool doSwitch = false;

    // This intrinsic edge is original on both P & Q simultaneously; the orientation judgement may be in conflict between P & Q,
    // so we adopt the judgement in A
    if (P_original && Q_original) {
      if (cointriP.mdata->shortName == "A") {
        if (P_reversed)
          doSwitch = true;
      } else {
        if (Q_reversed)
          doSwitch = true;
      }

    // If P_eIntrinsic is original on P, orient it the same as P_eInput
    } else if (P_original) {
      if (P_reversed)
        doSwitch = true;

    // If Q_eIntrinsic is original on Q, orient P_eIntrinsic the same as Q_eInput
    } else if (Q_original) {
      if (Q_reversed)
        doSwitch = true;

    // For non-original edge, orient it according to the unique ID
    } else {
      if (cointriP.uniqueID_per_Vertex[P_vIntrinsic[0]] >= cointriP.uniqueID_per_Vertex[P_vIntrinsic[1]])
        doSwitch = true;
    }

    if (doSwitch)
      cointriP.signpostTri->switchHalfedgeSides(P_eIntrinsic);
  }
}

void cit::copyConfiguration(const Configuration& origConfig, Configuration& copiedConfig, size_t logDepth) {
  DLOG_INFO(logDepth, "BEGIN copyConfiguration");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_INFO(logDepth, "END   copyConfiguration"); });

  // Copy meta info
  copiedConfig.topologyValid = origConfig.topologyValid;
  // copiedConfig.angleValid = origConfig.angleValid;
  copiedConfig.x = origConfig.x;
  copiedConfig.energy = origConfig.energy;
  copiedConfig.nCompatibleEdges = origConfig.nCompatibleEdges;
  copiedConfig.nCompatibleFaces = origConfig.nCompatibleFaces;

  // Copy signpostTri & uniqueID_per_Vertex through serialization
  copiedConfig.cointriA.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataA.mesh, *mdataA.geometry, ::geometrycentral::toSerializedBlob(*origConfig.cointriA.signpostTri)));
  copiedConfig.cointriB.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataB.mesh, *mdataB.geometry, ::geometrycentral::toSerializedBlob(*origConfig.cointriB.signpostTri)));
  copiedConfig.cointriA.uniqueID_per_Vertex = VertexData<int>(*copiedConfig.cointriA.signpostTri->intrinsicMesh, ::geometrycentral::toSerializedBlob(origConfig.cointriA.uniqueID_per_Vertex));
  copiedConfig.cointriB.uniqueID_per_Vertex = VertexData<int>(*copiedConfig.cointriB.signpostTri->intrinsicMesh, ::geometrycentral::toSerializedBlob(origConfig.cointriB.uniqueID_per_Vertex));

  updateCorrespondence(copiedConfig, spdlog::default_logger(), logDepth + 1);
}

bool cit::snapFacePointToEdge(SurfacePoint& sp) {
  CIT_ASSERT(sp.type == SurfacePointType::Face);
  size_t minIndex;
  if (min(sp.faceCoords, &minIndex) < sysParam.snapThreshold) {
    sp.faceCoords[minIndex] = 0.;
    sp.faceCoords /= sum(sp.faceCoords);
    sp = sp.reduced();
    CIT_ASSERT(sp.type == SurfacePointType::Edge);
    return true;
  }
  return false;
}

double cit::estimateMaxStep(Configuration& config, double wL, size_t logDepth) {
  DLOG_TRACE(logDepth, "BEGIN estimateMaxStep");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(logDepth, "END   estimateMaxStep"); });

  VectorXd deltaBC = computeDescentDirection(config, {}, wL, 0, logDepth + 1);
  VectorXd delta = config.S * deltaBC;     // Convert barycentric coordinate differential to a 2D vector whose length is in object space

  // Get the largest displacement
  double maxDisplacement = 0.;
  CIT_ASSERT((size_t)delta.size() == 2 * (mdataA.nV + mdataB.nV));
  for (size_t i = 0; i < mdataA.nV + mdataB.nV; ++i) {
    Vector2d d = delta.segment(2 * i, 2);
    maxDisplacement = std::max<double>(maxDisplacement, d.norm());
  }

  double smax = 0.01 / maxDisplacement;
  DLOG_INFO(logDepth + 1, "The largest displacement magnitude is {}, so we estimate smax to {}", maxDisplacement, smax);
  return smax;
}

// Useful for debugging
std::set<Edge> cit::getCompatibleEdges(const CoInTri& cointri) {
  std::set<Edge> res;
  for (Edge e : cointri.signpostTri->intrinsicMesh->edges()) {
    if (isEdgeCompatible(cointri, e))
      res.insert(e);
  }
  return res;
}
std::set<Edge> cit::getIncompatibleEdges(const CoInTri& cointri) {
  std::set<Edge> res;
  for (Edge e : cointri.signpostTri->intrinsicMesh->edges()) {
    if (!isEdgeCompatible(cointri, e))
      res.insert(e);
  }
  return res;
}
std::set<Face> cit::getCompatibleFaces(const CoInTri& cointri) {
  std::set<Face> res;
  for (Face e : cointri.signpostTri->intrinsicMesh->faces()) {
    if (isFaceCompatible(cointri, e))
      res.insert(e);
  }
  return res;
}
std::set<Face> cit::getIncompatibleFaces(const CoInTri& cointri) {
  std::set<Face> res;
  for (Face e : cointri.signpostTri->intrinsicMesh->faces()) {
    if (!isFaceCompatible(cointri, e))
      res.insert(e);
  }
  return res;
}
std::set<Face> cit::getReversedFaces(const CoInTri& cointri) {
  std::set<Face> res;
  for (Face e : cointri.signpostTri->intrinsicMesh->faces()) {
    if (isFaceReversed(cointri, e))
      res.insert(e);
  }
  return res;
}
std::vector<Vertex> cit::adjacentVertices(Vertex v) {
  std::vector<Vertex> res;
  for (Vertex v2 : v.adjacentVertices())
    res.push_back(v2);
  return res;
}
std::vector<Halfedge> cit::outgoingHalfedges(Vertex v) {
  std::vector<Halfedge> res;
  for (Halfedge he : v.outgoingHalfedges())
    res.push_back(he);
  return res;
}
std::vector<Face> cit::adjacentFaces(Vertex v) {
  std::vector<Face> res;
  for (Face f: v.adjacentFaces())
    res.push_back(f);
  return res;
}
std::array<Face, 2> cit::adjacentFaces(Edge e) {
  return {
    e.halfedge().face(),
    e.halfedge().twin().face()
  };
}
std::array<Vertex, 3> cit::adjacentVertices(Face f) {
  return {
    f.halfedge().vertex(),
    f.halfedge().next().vertex(),
    f.halfedge().next().next().vertex()
  };
}
std::array<Halfedge, 3> cit::adjacentHalfedges(Face f) {
  return {
    f.halfedge(),
    f.halfedge().next(),
    f.halfedge().next().next()
  };
}
std::array<Face, 3> cit::adjacentFaces(Face f) {
  return {
    f.halfedge().twin().face(),
    f.halfedge().next().twin().face(),
    f.halfedge().next().next().twin().face()
  };
}
