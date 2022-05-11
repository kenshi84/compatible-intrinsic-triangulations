#include "cit.hpp"

#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>

std::string cit::toSerializedBlob(const Configuration& config, const std::list<std::pair<VectorXd, SparseMatrixd>>& history, double wL, double smax, bool firstOrderMode, size_t nSteps) {
  std::array<std::string, 11> splitBlobs = {
    ::geometrycentral::toSerializedBlob(*config.cointriA.signpostTri),
    ::geometrycentral::toSerializedBlob(config.cointriA.uniqueID_per_Vertex),
    ::geometrycentral::toSerializedBlob(*config.cointriB.signpostTri),
    ::geometrycentral::toSerializedBlob(config.cointriB.uniqueID_per_Vertex),
    ::geometrycentral::toSerializedBlob(config.nCompatibleEdges),
    ::geometrycentral::toSerializedBlob(config.nCompatibleFaces),
    ::geometrycentral::toSerializedBlob(history),
    ::geometrycentral::toSerializedBlob(wL),
    ::geometrycentral::toSerializedBlob(smax),
    ::geometrycentral::toSerializedBlob(firstOrderMode),
    ::geometrycentral::toSerializedBlob(nSteps)
  };

  std::pair<unsigned char, std::string> versionedBlob{ 2, ::geometrycentral::toSerializedBlob(splitBlobs) };
  return ::geometrycentral::toSerializedBlob(versionedBlob);
}

void cit::fromSerializedBlob(const std::string& serializedBlob, Configuration& config, std::list<std::pair<VectorXd, SparseMatrixd>>& history, double& wL, double& smax, bool& firstOrderMode, size_t& nSteps) {
  std::array<std::string, 11> splitBlobs;

  // Set default values to newly introduced items
  splitBlobs[9 ] = ::geometrycentral::toSerializedBlob<bool>(false);  // firstOrderMode
  splitBlobs[10] = ::geometrycentral::toSerializedBlob<size_t>(0);    // nSteps

  // Try to deserialize assuming the blob is versioned
  std::pair<unsigned char, std::string> versionedBlob;
  try {
    ::geometrycentral::fromSerializedBlob(serializedBlob, versionedBlob);

    if (versionedBlob.first == 2) {
      ::geometrycentral::fromSerializedBlob(versionedBlob.second, splitBlobs);

    } else if (versionedBlob.first == 1) {
      std::array<std::string, 10> splitBlobs10;
      ::geometrycentral::fromSerializedBlob(versionedBlob.second, splitBlobs10);
      std::copy(splitBlobs10.begin(), splitBlobs10.end(), splitBlobs.begin());

    } else {
      DLOG_ERROR(0, "Unknown blob version");
      CIT_ASSERT(false);
    }

  // Retry assuming the blob is in legacy unversioned format with 9 items
  } catch (const std::exception& e) {
    DLOG_INFO(0, "Error while trying to deserialize assuming versioned blob: {}", e.what());
    DLOG_INFO(0, "Retrying assuming legacy unversioned blob ...");

    try {
      std::array<std::string, 9> splitBlobs9;
      ::geometrycentral::fromSerializedBlob(serializedBlob, splitBlobs9);
      std::copy(splitBlobs9.begin(), splitBlobs9.end(), splitBlobs.begin());

      DLOG_INFO(0, "Succeeded in deserializing legacy unversioned blob!");

    } catch (const std::exception& e) {
      DLOG_ERROR(0, "Error while deserializing: {}", e.what());
      throw;
    }
  }

  config = {};
  CoInTri& cointriA = config.cointriA;
  CoInTri& cointriB = config.cointriB;

  cointriA.mdata = &mdataA;
  cointriB.mdata = &mdataB;

  cointriA.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataA.mesh, *mdataA.geometry, splitBlobs[0]));
  cointriB.signpostTri.reset(new SignpostIntrinsicTriangulation(*mdataB.mesh, *mdataB.geometry, splitBlobs[2]));

  cointriA.uniqueID_per_Vertex = VertexData<int>(*cointriA.signpostTri->intrinsicMesh, splitBlobs[1]);
  cointriB.uniqueID_per_Vertex = VertexData<int>(*cointriB.signpostTri->intrinsicMesh, splitBlobs[3]);

  // This suppresses cout printing
  ::geometrycentral::fromSerializedBlob(splitBlobs[4], config.nCompatibleEdges);
  ::geometrycentral::fromSerializedBlob(splitBlobs[5], config.nCompatibleFaces);

  config.topologyValid = config.nCompatibleEdges == cointriA.signpostTri->intrinsicMesh->nEdges();

  // double A_minAngle = getMinAngle(config.cointriA);
  // double B_minAngle = getMinAngle(config.cointriB);
  // double minAngle = std::min<double>(A_minAngle, B_minAngle);
  // config.angleValid = minAngle > sysParam.angleThreshold;

  updateCorrespondence(config);

  computeEnergy(config);

  setSolutionVector(config);

  ::geometrycentral::fromSerializedBlob(splitBlobs[6], history);
  ::geometrycentral::fromSerializedBlob(splitBlobs[7], wL);
  ::geometrycentral::fromSerializedBlob(splitBlobs[8], smax);
  ::geometrycentral::fromSerializedBlob(splitBlobs[9], firstOrderMode);
  ::geometrycentral::fromSerializedBlob(splitBlobs[10], nSteps);
}
