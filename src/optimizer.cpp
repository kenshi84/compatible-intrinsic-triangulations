#include "cit.hpp"
#include "imgui.h"
#include "imfilebrowser.h"

// #include <kt84/tbb_util.hh>

namespace {
int initTweaks = [](){
  cit::tweaks.d["stepSize"] = 1.e-3;
  cit::tweaks.b["updateConfiguration"] = false;
  cit::tweaks.b["dontUpdateConfigurationIfWorse"] = true;
  cit::tweaks.b["dontUpdateConfigurationIfInvalid"] = true;
  cit::tweaks.b["doTemporalSmoothing"] = true;
  return 0;
}();
}

extern VectorXd dbg_descentDirection;

std::deque<VectorXsp> dbg_previousX;
std::deque<VectorXd> dbg_previousGradient;
std::deque<VertexData<std::vector<SurfacePoint>>> dbg_vertexPathA;
std::deque<VertexData<std::vector<SurfacePoint>>> dbg_vertexPathB;

namespace cit {

void updateDescentDirection() {
  computeDescentDirection(currentConfig, tweaks.b["doTemporalSmoothing"] ? currentHistory : std::list<std::pair<VectorXd, SparseMatrixd>>{}, current_wL, 0);
  mdataA.dispList.invalidate();
  mdataB.dispList.invalidate();
}

}

bool cit::findMaxStep(const VectorXd& descentDirection, int l) {
  DLOG_INFO(0, "");
  DLOG_INFO(0, "+-----------------+");
  DLOG_INFO(0, "| findMaxStep ... |");
  DLOG_INFO(0, "+-----------------+");

  uint64_t id;
  Configuration newConfig;
  bool upward;

  while (true) {
    DLOG_INFO(0, "Examining step size current_smax = {} ... ", current_smax);
    id = makeConfigurationID();
    upward = createConfiguration(id, currentConfig, current_smax * descentDirection, newConfig, makeLogger(getResultFileName(ResultContext::FindMaxStep, id, currentConfig.id, l) + ".log", true));
    fillResultInfo(newConfig, currentConfig.energy, currentConfig.gradient, descentDirection, l, current_smax);
    writeResult(newConfig, {}, 0., 0., 0, ResultContext::FindMaxStep);

    if (upward) {
      DLOG_INFO(0, "current_smax is rather small; increase it as long as feasible");

    } else {
      size_t nIncompatibleEdges = newConfig.cointriA.signpostTri->intrinsicMesh->nEdges() - newConfig.nCompatibleEdges;
      if (nIncompatibleEdges > 100) {
        DLOG_INFO(0, "current_smax is probably too large! scaling it down substantially");
        current_smax *= 1e-10;
        continue;
      }
      DLOG_INFO(0, "current_smax is rather large; decrease it until reaching feasible");
    }
    break;
  }

  // preparing for parallel execution
  std::vector<double> parallel_smax(sysParam.nThreads);
  parallel_smax[0] = (upward ? sysParam.smax_scaleFactor : (1. / sysParam.smax_scaleFactor)) * current_smax;
  for (unsigned int i = 1; i < sysParam.nThreads; ++i) {
    parallel_smax[i] = (upward ? sysParam.smax_scaleFactor : (1. / sysParam.smax_scaleFactor)) * parallel_smax[i - 1];
  }
  for (size_t nIter = 0; nIter < sysParam.nMaxIter; ++nIter) {
    DLOG_INFO(0, "Running iteration #{}...", nIter);
    // parallel execution
    std::array<uint64_t          , MAX_NTHREADS> parallel_id;
    std::array<char              , MAX_NTHREADS> parallel_feasible;     parallel_feasible.fill(0);
    std::array<Configuration     , MAX_NTHREADS> parallel_newConfig;    // use std::array instead of std::vector which seems buggy (when combind with std::shared_ptr?) in Ubuntu
    std::array<std::ostringstream, MAX_NTHREADS> parallel_oss;
    std::array<spdlog::logger_ptr, MAX_NTHREADS> parallel_logger;

    for (unsigned int i = 0; i < sysParam.nThreads; ++i)
      parallel_id[i] = makeConfigurationID();

    auto findMaxStep_threadFunc = [&](size_t i) {
      parallel_logger[i]->info("Examining step size {} ... ", parallel_smax[i]);
      parallel_feasible[i] = createConfiguration(parallel_id[i], currentConfig, parallel_smax[i] * descentDirection, parallel_newConfig[i], parallel_logger[i]);
    };

    if (sysParam.nThreads == 1) {
      // no parallel execution
      parallel_logger[0] = makeLogger(getResultFileName(ResultContext::FindMaxStep, parallel_id[0], currentConfig.id, l) + ".log", true);
      findMaxStep_threadFunc(0);

    } else {
      // prepare logger for each thread
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        parallel_logger[i] = makeLogger(getResultFileName(ResultContext::FindMaxStep, parallel_id[i], currentConfig.id, l) + ".log", false, &parallel_oss[i]);

      // // parallel execution
      // kt84::tbb_util::parallel_for(sysParam.nThreads, findMaxStep_threadFunc, 1);

      // PARALLEL EXECUTION DOESN'T WORK!
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        findMaxStep_threadFunc(i);
    }

    // print output at once
    if (sysParam.nThreads > 1) {
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        DLOG_INFO(0, parallel_oss[i].str());
    }

    // Write results to file
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      fillResultInfo(parallel_newConfig[i], currentConfig.energy, currentConfig.gradient, descentDirection, l, parallel_smax[i]);
      writeResult(parallel_newConfig[i], {}, 0., 0., 0, ResultContext::FindMaxStep);
    }

    // check parallel execution results
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      if ((upward && !parallel_feasible[i]) || (!upward && parallel_feasible[i])) {
        if (!upward) {
          current_smax = parallel_smax[i];
        }
        DLOG_INFO(0, "+-----------------------------------");
        DLOG_INFO(0, "| findMaxStep done! ");
        DLOG_INFO(0, "|   current_smax = {}", current_smax);
        DLOG_INFO(0, "+-----------------------------------");
        DLOG_INFO(0, "");
        return true;
      }
      current_smax = parallel_smax[i];
    }
    if (current_smax < 1.e-30) {
      DLOG_WARN(0, "+---------------------------------------");
      DLOG_WARN(0, "| findMaxStep failed:");
      DLOG_WARN(0, "|   Couldn't find any feasible step size");
      DLOG_WARN(0, "+---------------------------------------");
      return false;
    }
    // prepare next batch of parallel_smax
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      parallel_smax[i] *= std::pow(upward ? sysParam.smax_scaleFactor : (1. / sysParam.smax_scaleFactor), sysParam.nThreads);
    }
  }
  DLOG_ERROR(0, "Couldn't find smax after the maximum number of iterations!");
  return false;
}

bool cit::findOptimalStep(const VectorXd& descentDirection, int l, double& s, Configuration& newConfig) {
  DLOG_INFO(0, "");
  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "| findOptimalStep ... |");
  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "");

  double slope = currentConfig.gradient.dot(descentDirection);

  std::vector<double> parallel_prevEnergy(sysParam.nThreads, 0.);     // To check that energy is still changing

  // -- preparing for parallel execution
  std::vector<double> parallel_s(sysParam.nThreads);
  parallel_s[0] = current_smax;
  for (unsigned int i = 1; i < sysParam.nThreads; ++i) {
    parallel_s[i] = sysParam.s_scaleFactor * parallel_s[i - 1];
  }
  for (size_t nIter = 0; nIter < sysParam.nMaxIter; ++nIter) {
    DLOG_INFO(0, "Running iteration #{}...", nIter);
    // -- parallel execution
    std::array<uint64_t          , MAX_NTHREADS> parallel_id;
    std::array<char              , MAX_NTHREADS> parallel_satisfied;      parallel_satisfied.fill(0);
    std::array<Configuration     , MAX_NTHREADS> parallel_newConfig;      // use std::array instead of std::vector which seems buggy (when combind with std::shared_ptr?) in Ubuntu
    std::array<std::ostringstream, MAX_NTHREADS> parallel_oss;
    std::array<spdlog::logger_ptr, MAX_NTHREADS> parallel_logger;

    for (unsigned int i = 0; i < sysParam.nThreads; ++i)
      parallel_id[i] = makeConfigurationID();

    auto findOptimalStep_threadFunc = [&](size_t i) {
      parallel_logger[i]->info("Examining step size {} ... ", parallel_s[i]);
      if (createConfiguration(parallel_id[i], currentConfig, parallel_s[i] * descentDirection, parallel_newConfig[i], parallel_logger[i])) {
        if (parallel_newConfig[i].energy < currentConfig.energy + sysParam.armijo * parallel_s[i] * slope) {
          parallel_satisfied[i] = true;
        }
      }
    };

    if (sysParam.nThreads == 1) {
      // no parallel execution
      parallel_logger[0] = makeLogger(getResultFileName(ResultContext::FindOptimalStep, parallel_id[0], currentConfig.id, l) + ".log", true);
      findOptimalStep_threadFunc(0);

    } else {
      // prepare logger for each thread
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        parallel_logger[i] = makeLogger(getResultFileName(ResultContext::FindOptimalStep, parallel_id[i], currentConfig.id, l) + ".log", false, &parallel_oss[i]);

      // // parallel execution
      // kt84::tbb_util::parallel_for(sysParam.nThreads, findOptimalStep_threadFunc, 1);

      // PARALLEL EXECUTION DOESN'T WORK!
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        findOptimalStep_threadFunc(i);
    }

    // print output at once
    if (sysParam.nThreads > 1) {
      for (unsigned int i = 0; i < sysParam.nThreads; ++i)
        DLOG_INFO(0, parallel_oss[i].str());
    }

    // Write results to file
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      fillResultInfo(parallel_newConfig[i], currentConfig.energy, currentConfig.gradient, descentDirection, l, parallel_s[i]);
      writeResult(parallel_newConfig[i], {}, 0., 0., 0, ResultContext::FindOptimalStep);
    }

    // check parallel execution results
    bool changed = false;
    bool found = false;
    kt84::MinSelector<int> chosen;
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      if (parallel_satisfied[i]) {
        found = true;
        chosen.update(parallel_newConfig[i].energy, i);
      }
      if (parallel_newConfig[i].energy == std::numeric_limits<double>::infinity() || parallel_prevEnergy[i] != parallel_newConfig[i].energy)
        changed = true;
    }
    if (found) {
      s = parallel_s[chosen.value];
      newConfig = parallel_newConfig[chosen.value];
      DLOG_INFO(0, "+---------------------------------------");
      DLOG_INFO(0, "| findOptimalStep done!");
      DLOG_INFO(0, "|   s = {}", s);
      DLOG_INFO(0, "|   energy delta = {}", newConfig.energy - currentConfig.energy);
      DLOG_INFO(0, "+---------------------------------------");
      DLOG_INFO(0, "");
      return true;
    }
    if (!changed) {
      DLOG_WARN(0, "+---------------------------------------");
      DLOG_WARN(0, "| findOptimalStep failed:");
      DLOG_WARN(0, "|   No change in energy");
      DLOG_WARN(0, "+---------------------------------------");
      return false;
    }
    // prepare next batch of parallel_s
    for (unsigned int i = 0; i < sysParam.nThreads; ++i) {
      parallel_prevEnergy[i] = parallel_newConfig[i].energy;
      parallel_s[i] *= std::pow(sysParam.s_scaleFactor, sysParam.nThreads);
    }
    if (descentDirection.norm() * parallel_s[0] < 1.e-30) {
      DLOG_WARN(0, "+---------------------------------------");
      DLOG_WARN(0, "| findOptimalStep failed:");
      DLOG_WARN(0, "|   Step size became too small");
      DLOG_WARN(0, "+---------------------------------------");
      return false;
    }
  }
  DLOG_ERROR(0, "Couldn't find optimal step size after the maximum number of iterations!");
  return false;
}

void cit::updateConfiguration(const Configuration& newConfig, int l_chosen) {
  DLOG_TRACE(0, "BEGIN updateConfiguration");
  auto scopeExit = kt84::make_ScopeExit([&](){ DLOG_TRACE(0, "END   updateConfiguration"); });
  DLOG_INFO(0, "  Energy delta: {}", newConfig.energy - currentConfig.energy);

  // append currentConfig data to history
  currentHistory.push_front({currentConfig.gradient, currentConfig.hessian});
  dbg_previousX.push_front(currentConfig.x);
  dbg_previousGradient.push_front(currentConfig.gradient);
  dbg_vertexPathA.push_front(currentConfig.cointriA.vertexPath);
  dbg_vertexPathB.push_front(currentConfig.cointriB.vertexPath);

  while (currentHistory.size() > sysParam.historyMaxSize) {
    currentHistory.pop_back();
    dbg_previousX.pop_back();
    dbg_previousGradient.pop_back();
    dbg_vertexPathA.pop_back();
    dbg_vertexPathB.pop_back();
  }

  currentConfig = newConfig;

  // Needed for rendering, even when the new configuration is invalid
  computeEdgePath(currentConfig, 1);

  const ResultInfo& ri = resultInfoTable.at(currentConfig.id);
  if (ri.A_maxSignpostError > sysParam.maxSignpostError || ri.B_maxSignpostError > sysParam.maxSignpostError) {
    DLOG_INFO(0, "Accumulated signpost error became too high: {} on A, {} on B", ri.A_maxSignpostError, ri.B_maxSignpostError);
    DLOG_INFO(0, "  Sanitizing ...");
    std::thread tA([](){currentConfig.cointriA.signpostTri->sanitizeSignpost(*mdataA.geometry, currentConfig.cointriA.intrinsicEdgePath);});
    std::thread tB([](){currentConfig.cointriB.signpostTri->sanitizeSignpost(*mdataB.geometry, currentConfig.cointriB.intrinsicEdgePath);});
    tA.join();
    tB.join();
    DLOG_INFO(0, "  Computing edge path again with sanitized signpost data");
    computeEdgePath(currentConfig, 1);
    DLOG_INFO(0, "  Signpost error after sanitization: {} on A, {} on B", ri.A_maxSignpostError, ri.B_maxSignpostError);
  }

  computeOverlayPolygons(currentConfig, 1);

  if (currentConfig.topologyValid) {
    computeTemporalTransportOperator(currentConfig, 1);

    // apply temporal transport operator to each entry in the derivative history
    for (auto& p : currentHistory) {
      p.first = currentConfig.T * p.first;
      p.second = SparseMatrixd(currentConfig.Tinv.transpose()) * p.second * currentConfig.Tinv;
    }

    computeDescentDirection(currentConfig, tweaks.b["doTemporalSmoothing"] ? currentHistory : std::list<std::pair<VectorXd, SparseMatrixd>>{}, current_wL, 0, 1);
  }

  // also update wL
  if (current_wL && !sysParam.firstOrderMode)
    current_wL *= std::pow(sysParam.alpha, l_chosen);

  ++current_nSteps;

  writeResult(currentConfig, currentHistory, current_wL, current_smax, current_nSteps, ResultContext::UpdateConfiguration);

  if (currentConfig.topologyValid)
    writeOverlayMesh(currentConfig);

  mdataA.dispList.invalidate();
  mdataB.dispList.invalidate();
}

bool cit::optimizeOneStep(bool doTemporalSmoothing) {
  DLOG_INFO(0, "");
  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "| optimizeOneStep ...");
  DLOG_INFO(0, "|   current_wL = {}", current_wL);
  DLOG_INFO(0, "|   current energy = {}", currentConfig.energy);
  DLOG_INFO(0, "|   nSteps = {}", current_nSteps);
  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "");
  kt84::Timer timer;

  double prevEnergy = currentConfig.energy;

  std::array<Configuration, 3> newConfig;
  std::array<double, 3> s;
  for (int l = -1; l <= 1; ++l) {
    DLOG_INFO(0, "Examining l = {} ... ", l);
    VectorXd descentDirection;
    double slope;
    while (true) {
      descentDirection = computeDescentDirection(currentConfig, doTemporalSmoothing ? currentHistory : decltype(currentHistory){}, current_wL, l, 1);

      slope = currentConfig.gradient.dot(descentDirection);

      // reset history if smoothing leads to non-descent direction or the smoothed gradient contradicts the current gradient
      if (slope > 0 || currentConfig.gradient.dot(currentConfig.smoothedGradient) < 0) {
        CIT_ASSERT(doTemporalSmoothing && !currentHistory.empty());
        DLOG_WARN(0, "");
        DLOG_WARN(0, "");
        DLOG_WARN(0, "### Reset history as smoothing generated bad descent direction ###");
        DLOG_WARN(0, "");
        currentHistory.clear();
        continue;
      }
      break;
    }

    if (!findMaxStep(descentDirection, l)) {
      DLOG_WARN(0, "findMaxStep failed!");
      continue;
    }

    if (!findOptimalStep(descentDirection, l, s[l + 1], newConfig[l + 1]))
      DLOG_WARN(0, "findOptimalStep failed!");

    if (!current_wL || sysParam.firstOrderMode)
      break;  // don't repeat when Laplacian preconditioner is disabled, or we're in first-order mode
  }

  // replace currentConfig with best newConfig
  kt84::MinSelector<int> l_chosen;
  for (int l = -1; l <= 1; ++l) {
    if (newConfig[l + 1].energy < std::numeric_limits<double>::infinity())
      l_chosen.update(newConfig[l + 1].energy, l);
  }
  if (!l_chosen.count) {
    DLOG_INFO(0, "+---------------------+");
    DLOG_INFO(0, "| optimizeOneStep failed:");
    DLOG_INFO(0, "|   None of l in {-1,0,1} resulted in valid next step");
    DLOG_INFO(0, "+---------------------+");
    return false;
  }
  updateConfiguration(newConfig[l_chosen.value + 1], l_chosen.value);
  fs::copy_file(
    getResultFileName(ResultContext::FindOptimalStep, currentConfig.id, currentConfig.previous_id, l_chosen.value) + ".log",
    getResultFileName(ResultContext::UpdateConfiguration, currentConfig.id) + ".log"
  );
  cit::tweaks.d["stepSize"] = s[l_chosen.value + 1];

  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "| optimizeOneStep done!");
  DLOG_INFO(0, "|   chosen l = {}", l_chosen.value);
  DLOG_INFO(0, "|   (new) current_wL = {}", current_wL);
  DLOG_INFO(0, "|   (new) Current energy = {}", currentConfig.energy);
  DLOG_INFO(0, "|   Energy delta = {}", currentConfig.energy - prevEnergy);
  DLOG_INFO(0, "|   elapsed time (ms) = {}", timer.milliseconds_str());
  DLOG_INFO(0, "+---------------------+");
  DLOG_INFO(0, "");
  return true;
}

void cit::optimize(int N) {
  if (N == 0) return;

  size_t nSuccessPrev = 1;
  size_t nSuccess = 0;

  for (int i = 0; N == -1 || i < N; ++i) {
    double energyBefore = currentConfig.energy;

    bool success = optimizeOneStep(tweaks.b["doTemporalSmoothing"]);

    if (success) {
      double energyAfter = currentConfig.energy;
      success = sysParam.energyDeltaThreshold < energyBefore - energyAfter;

      if (!success)
        DLOG_ERROR(0, "Energy decrease is very small!");

    } else {
      DLOG_ERROR(0, "Some error occurred in optimizeOneStep!");
    }

    if (success) {
      DLOG_INFO(0, "Step {} was successful! Moving on to the next step...", current_nSteps);
      ++nSuccess;
      continue;
    }

    // If we haven't had any success in this mode and in the other mode, give up
    if (!nSuccess && !nSuccessPrev) {
      DLOG_INFO(0, "");
      DLOG_INFO(0, "");
      DLOG_INFO(0, "!!!! Optimization seems to have converged !!!!");
      break;
    }

    // Otherwise, switch to the other mode
    DLOG_INFO(0, "We're stuck in this {}-order mode! Switching to the other {}-order mode...",
      sysParam.firstOrderMode ? "first" : "second",
      sysParam.firstOrderMode ? "second" : "first");

    sysParam.firstOrderMode = !sysParam.firstOrderMode;
    nSuccessPrev = nSuccess;
    nSuccess = 0;

    current_smax = estimateMaxStep(currentConfig, current_wL, 1);
  }
}

void cit::optimizerMenu() {

  ImGui::Begin("Optimizer");

  // ImGui::SliderInt("nThreads", (int*)&sysParam.nThreads, 1, std::thread::hardware_concurrency());
  // ImGui::Separator();

  ImGui::Text("#Steps: %lu", current_nSteps);
  ImGui::Text("Energy: %.10f", currentConfig.energy);
  ImGui::Text("History size: %lu", currentHistory.size());
  ImGui::Checkbox("Do temporal smoothing", &tweaks.b["doTemporalSmoothing"]);

  ImGui::Separator();

  if (ImGui::Checkbox("First-order mode", &sysParam.firstOrderMode)) {
    updateDescentDirection();
  }

  ImGui::Separator();
  if (ImGui::Button("Export overlay mesh") && currentConfig.topologyValid) {
    writeOverlayMesh(currentConfig);
  }

  ImGui::Separator();

  static int N = 1;
  ImGui::DragInt("N", &N, 5.f, -1, 1000000);
  ImGui::SliderScalar("armijo", ImGuiDataType_Double, &sysParam.armijo, staticValuePtr(0.), staticValuePtr(1.), "%.10f", 5.f);
  if (ImGui::Button("Optimize N steps")) {
    optimize(N);
  }

  ImGui::Separator();
  if (ImGui::SliderScalar("Current wL", ImGuiDataType_Double, &current_wL, staticValuePtr(0.), staticValuePtr(1.e+10), "%.10f", 5.f)) {
    updateDescentDirection();
  }
  if (ImGui::Button("Increase wL")) {
    current_wL *= sysParam.alpha;
    updateDescentDirection();
  }
  if (ImGui::Button("Decrease wL")) {
    current_wL /= sysParam.alpha;
    updateDescentDirection();
  }

  ImGui::Separator();
  if (ImGui::Button("Find max step size")) {
    findMaxStep(dbg_descentDirection, 0);
  }
  ImGui::SliderScalar("Current smax", ImGuiDataType_Double, &current_smax, staticValuePtr(0.), staticValuePtr(100.), "%.30f", 100.f);

  ImGui::Separator();
  ImGui::Checkbox("Update configuration", &tweaks.b["updateConfiguration"]);

  ImGui::Separator();
  if (ImGui::Button("Find optimal step size")) {
    Configuration newConfig;
    findOptimalStep(dbg_descentDirection, 0, tweaks.d["stepSize"], newConfig);
    if (tweaks.b["updateConfiguration"]) {
      updateConfiguration(newConfig, 0);
    }
  }
  ImGui::SliderScalar("Step size", ImGuiDataType_Double, &tweaks.d["stepSize"], staticValuePtr(1.e-10), staticValuePtr(10.), "%.30f", 100.f);
  ImGui::Checkbox("Don't update config if worse", &tweaks.b["dontUpdateConfigurationIfWorse"]);
  ImGui::Checkbox("Don't update config if invalid", &tweaks.b["dontUpdateConfigurationIfInvalid"]);
  if (ImGui::Button("Manual stepping")) {
    if (!currentConfig.topologyValid)
      return;

    DLOG_INFO(0, "");
    DLOG_INFO(0, "Manual stepping with s = {}", tweaks.d["stepSize"]);

    Configuration newConfig;
    uint64_t id = makeConfigurationID();
    bool success = createConfiguration(id, currentConfig, tweaks.d["stepSize"] * dbg_descentDirection, newConfig, makeLogger(getResultFileName(ResultContext::UpdateConfiguration, id) + ".log", true));

    if (tweaks.b["updateConfiguration"] && ((success && (newConfig.energy < currentConfig.energy || !tweaks.b["dontUpdateConfigurationIfWorse"])) || !tweaks.b["dontUpdateConfigurationIfInvalid"])) {

      VectorXsp previousX = currentConfig.x;
      VectorXd previousGradient = currentConfig.gradient;

      fillResultInfo(newConfig, currentConfig.energy, currentConfig.gradient, dbg_descentDirection, 0, tweaks.d["stepSize"]);
      updateConfiguration(newConfig, 0);

#if 1
      if (currentConfig.topologyValid) {
        kt84::MaxMinAverage lengthError, angleChange;
        auto testTemporalTransportOperator = [&previousX, &previousGradient, &lengthError, &angleChange](const ModelData& mdataP, const ModelData& mdataQ, size_t offset) {
          for (size_t i = offset; i < offset + mdataP.nV; ++i) {
            // previous gradient vector in previous position, in 3D
            if (previousX[i].type == SurfacePointType::Vertex) {
              CIT_ASSERT(currentConfig.x[i].type == SurfacePointType::Vertex);
              continue;
            }

            if (previousX[i].type == SurfacePointType::Edge)
              previousX[i] = convertEdgePointToFacePoint(previousX[i]);

            SurfacePoint previousGradientSP(previousX[i].face, {});
            for (size_t di = 0; di < 2; ++di)
              previousGradientSP.faceCoords[di] = previousGradient[2 * i + di];
            previousGradientSP.faceCoords[2] = -(previousGradientSP.faceCoords[0] + previousGradientSP.faceCoords[1]);
            Vector3 previousGradient_3D = previousGradientSP.interpolate(mdataQ.geometry->inputVertexPositions);

            // transported previous gradient vector in current position, in 3D
            SurfacePoint transportedGradientSP;

            if (currentConfig.x[i].type == SurfacePointType::Face)
              transportedGradientSP = currentConfig.x[i];
            else
              transportedGradientSP = convertEdgePointToFacePoint(currentConfig.x[i]);

            for (size_t di = 0; di < 2; ++di)
              transportedGradientSP.faceCoords[di] = currentHistory.front().first[2 * i + di];
            transportedGradientSP.faceCoords[2] = -(transportedGradientSP.faceCoords[0] + transportedGradientSP.faceCoords[1]);
            Vector3 transportedGradient_3D = transportedGradientSP.interpolate(mdataQ.geometry->inputVertexPositions);

            lengthError.update(std::abs(previousGradient_3D.norm() - transportedGradient_3D.norm()));
            angleChange.update(angle(previousGradient_3D, transportedGradient_3D));
          }
        };
        testTemporalTransportOperator(mdataA, mdataB, 0);
        testTemporalTransportOperator(mdataB, mdataA, mdataA.nV);

        DLOG_INFO(0, "");
        DLOG_INFO(0, "");
        DLOG_INFO(0, "!!! Testing temporal transport operator... ");
        DLOG_INFO(0, "!!! Length error: max={}, min={}, avg={}", lengthError.max(), lengthError.min(), lengthError.avg());
        DLOG_INFO(0, "!!! Angle change: max={}, min={}, avg={}", angleChange.max(), angleChange.min(), angleChange.avg());
        DLOG_INFO(0, "");
      }
#endif
    }
  }

  ImGui::Separator();
  if (ImGui::Button("Reset history")) {
    currentHistory.clear();
  }
  if (ImGui::Button("Restore init")) {
    currentConfig = initialConfig;
    currentHistory.clear();
    current_wL = sysParam.wL_initial;
    computeDescentDirection(currentConfig, currentHistory, current_wL, 0);
    mdataA.dispList.invalidate();
    mdataB.dispList.invalidate();
  }
  static ImGui::FileBrowser fileDialog;
  fileDialog.SetTitle("Deserialize");
  fileDialog.SetTypeFilters({ ".dat" });
  if (ImGui::Button("Deserialize")) {
    fileDialog.Open();
  }
  ImGui::End();
  fileDialog.Display();
  if(fileDialog.HasSelected()) {
    std::string filename = fileDialog.GetSelected().string();
    DLOG_INFO(0, "");
    DLOG_INFO(0, "");
    DLOG_INFO(0, "++++ Deserializing from {}", filename);
    DLOG_INFO(0, "");
    try {
      fromSerializedBlob(readBlobFromFile(filename), currentConfig, currentHistory, current_wL, current_smax, sysParam.firstOrderMode, current_nSteps);
      currentHistory.clear();
      computeEdgePath(currentConfig);
      computeOverlayPolygons(currentConfig);
      computeDescentDirection(currentConfig, currentHistory, current_wL, 0);
      mdataA.dispList.invalidate();
      mdataB.dispList.invalidate();
    } catch (const std::exception& e) {
      DLOG_ERROR(0, "Failed to deserialize: {}", e.what());
    }
    fileDialog.ClearSelected();
  }
}
