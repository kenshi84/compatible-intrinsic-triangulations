#include "cit.hpp"

namespace cit {

Tweaks tweaks;

std::array<LightParam, 8> lightParam;

MaterialParam materialParam;

SystemParameter sysParam;

ModelData mdataA, mdataB;

Configuration currentConfig, initialConfig;

std::map<uint64_t, ResultInfo> resultInfoTable;

double current_wL = 0.;
double current_smax = 1.;
std::list<std::pair<VectorXd, SparseMatrixd>> currentHistory;
size_t current_nSteps = 0;

SparseMatrixd massMatrix;

Camera* camera_active = nullptr;
int window_width;
int window_height;
ImGuiIO* io;
uint64_t sessionID;
kt84::TextureObject colormapTexture;
std::string caseName;

int highlightedVertexID = 0;
std::pair<int, int> highlightedEdgeID = {0, 0};
int highlightedFaceID = -1;

}
