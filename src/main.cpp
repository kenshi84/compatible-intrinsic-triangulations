#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
// #include "tbb/tbb.h"
#pragma clang diagnostic pop

#include "cit.hpp"
using namespace cit;

#include "args/args.hxx"
#include "imgui/backends/imgui_impl_glfw.h"
#include "imgui/backends/imgui_impl_opengl2.h"
#include "kt84/glfw_util.hh"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#pragma clang diagnostic pop

#include "cit_resources.h"

int main(int argc, char** argv) {
  sessionID = (uint64_t)(std::chrono::system_clock::now().time_since_epoch() / std::chrono::seconds(1));

  // Configure the argument parser
  // clang-format off
  args::ArgumentParser parser("Compatible Intrinsic Triangulations Demo (version: " + std::string(VERSIONTAG) + ")");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::ValueFlag<std::string> argMapAtoB(parser, "file", "The barycentric coordinates of A's vertices on B", {"mapAtoB"});
  args::ValueFlag<std::string> argMapBtoA(parser, "file", "The barycentric coordinates of B's vertices on A", {"mapBtoA"});
  args::ValueFlag<std::string> argDeserialize(parser, "file", "Deserialize from file", {"deserialize"});
  args::ValueFlag<std::string> argTextureA(parser, "file", "Texture image file for A", {"textureA"});
  args::ValueFlag<std::string> argTextureB(parser, "file", "Texture image file for B", {"textureB"});
  // args::ValueFlag<unsigned int> argNThreads(parser, "number", "Number of threads for parallel execution", {"nThreads"});
  // args::ValueFlag<double> argAngleThreshold(parser, "number", "Minimum angle threshold", {"angleThreshold"});
  args::Flag argFixAnchors(parser, "", "Fix anchor vertices", {"fixAnchors"});
  args::ValueFlag<std::string> argCaseName(parser, "name", "Name for the experiment case", {"caseName"});
  args::ValueFlag<double> argSnapThreshold(parser, "number", "Threshold for snapping face points to edges", {"snapThreshold"});
  args::ValueFlag<double> argInitialWeightLap(parser, "number", "Initial value for Laplacian weight", {"initialWeightLap"});
  args::ValueFlag<double> argArgmijo(parser, "number", "Armijo constant", {"armijo"});
  args::ValueFlag<double> argEnergyDeltaThreshold(parser, "number", "Energy delta threshold", {"energyDeltaThreshold"});
  args::Flag argInteractive(parser, "", "Interavtive mode with OpenGL", {"interactive"});

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
  args::ValueFlag<std::string> argLogLevel(parser, "level", logLevelHelp, {"log-level"});

  args::Positional<std::string> argMeshA(parser, "meshA", "The mesh of model A");
  args::Positional<std::string> argMeshB(parser, "meshB", "The mesh of model B");

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
  if (!argMeshA || !argMeshB) {
    std::cout << parser;
    return EXIT_FAILURE;
  }

  // Some constants per model
  mdataA.shortName = "A";
  mdataB.shortName = "B";
  mdataA.fullName = guessNiceNameFromPath(args::get(argMeshA));
  mdataB.fullName = guessNiceNameFromPath(args::get(argMeshB));
  mdataA.inputEdgeColor = Vector3f{91.f, 155.f, 213.f} / 255.f;
  mdataB.inputEdgeColor = Vector3f{237.f, 125.f, 49.f} / 255.f;

  if (argCaseName)
    caseName = args::get(argCaseName);
  else
    caseName = mdataA.fullName + "-" + mdataB.fullName;

  // Configure logging
  spdlog::set_automatic_registration(false);
  spdlog::set_pattern("%D %T %L %s:%4# %v");
  spdlog::set_level(spdlog::level::trace);
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
      DLOG_ERROR(0, "Wrong log level was specified");
      std::cout << parser;
      return EXIT_FAILURE;
    }
  }
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(getResultFileNameBase() + ".log", true);
  file_sink->set_level(spdlog::level::trace);
  file_sink->set_pattern("%D %T %L %s:%4# %v");
  spdlog::default_logger()->sinks().push_back(file_sink);
  spdlog::flush_on(spdlog::level::info);

  DLOG_INFO(0, "Compatible Intrinsic Triangulations Demo (version: {})", VERSIONTAG);
  DLOG_INFO(0, "Session ID: {}", sessionID);

  // Print out all the arguments for easier debugging
  std::string argStr;
  {
    std::ostringstream oss;
    for (int i = 0; i < argc; ++i)
      oss << (i ? " " : "") << argv[i];
    argStr = oss.str();
  }
  DLOG_INFO(0, "Given arguments: {}", argStr);

  // Configure concurrency
  // if (argNThreads)
  //   sysParam.nThreads = std::max<unsigned int>(1, args::get(argNThreads));
  sysParam.nThreads = 1;
  if (sysParam.nThreads > MAX_NTHREADS) {
    DLOG_ERROR(0, "Number of threads must be no larger than {}", MAX_NTHREADS);
    return EXIT_FAILURE;
  }
  DLOG_INFO(0, "Number of threads: {}", sysParam.nThreads);
  // tbb::task_scheduler_init init(sysParam.nThreads);

  if (argFixAnchors)
    sysParam.fixAnchors = true;

#if 0
  if (argAngleThreshold) {
    sysParam.angleThreshold = args::get(argAngleThreshold);
    if (sysParam.angleThreshold <= 0) {
      DLOG_ERROR(0, "Non-positive angle threshold was specified");
      return EXIT_FAILURE;
    }
    if (sysParam.angleThreshold > 0.1) {
      DLOG_ERROR(0, "Angle threshold is probably too large");
      return EXIT_FAILURE;
    }
  }
#endif
  // sysParam.angleThreshold = 0.;

  if (argSnapThreshold) {
    sysParam.snapThreshold = args::get(argSnapThreshold);
    if (sysParam.snapThreshold < 0) {
      DLOG_ERROR(0, "Negative threshold for snapping was specified");
      return EXIT_FAILURE;
    }
    if (sysParam.snapThreshold > 0.5) {
      DLOG_ERROR(0, "Snapping threshold is probably too large");
      return EXIT_FAILURE;
    }
  }

  if (argInitialWeightLap) {
    sysParam.wL_initial = args::get(argInitialWeightLap);
    if (sysParam.wL_initial <= 0.) {
      DLOG_ERROR(0, "Non-positive value for initial wL was specified");
      return EXIT_FAILURE;
    }
  }

  if (argArgmijo) {
    sysParam.armijo = args::get(argArgmijo);
    if (sysParam.armijo < 0.) {
      DLOG_ERROR(0, "Negative value for Armijo constant was specified");
      return EXIT_FAILURE;
    }
  }

  if (argEnergyDeltaThreshold) {
    sysParam.energyDeltaThreshold = args::get(argEnergyDeltaThreshold);
    if (sysParam.energyDeltaThreshold <= 0.) {
      DLOG_ERROR(0, "Non-positive value for energy delta threshold was specified");
      return EXIT_FAILURE;
    }
  }

  // Read mesh from file and do some preprocessing
  try {
    readInputMesh(mdataA, args::get(argMeshA));
    readInputMesh(mdataB, args::get(argMeshB));
  } catch(const std::exception& e) {
    DLOG_ERROR(0, "Failed to read input mesh: {}", e.what());
    return EXIT_FAILURE;
  }

  if (mdataA.nV == mdataB.nV) {
    DLOG_ERROR(0, "The two models must have different number of vertices!");
    return EXIT_FAILURE;
  }

  // setup combined mass matrix
  massMatrix = SparseMatrixd(2 * (mdataA.nV + mdataB.nV), 2 * (mdataA.nV + mdataB.nV));
  {
    std::vector<Eigen::Triplet<double>> triplets;
    auto fillTriplets = [&triplets](const ModelData& mdata, size_t offset) {
      for (size_t i = 0; i < mdata.nV; ++i) {
        Vertex v = mdata.mesh->vertex(i);
        double vertexArea = mdata.geometry->vertexDualAreas[v];
        for (size_t di = 0; di < 2; ++di) {
          size_t idx = offset + 2 * i + di;
          triplets.emplace_back(idx, idx, vertexArea);
        }
      }
    };
    fillTriplets(mdataA, 0);
    fillTriplets(mdataB, 2 * mdataA.nV);
    massMatrix.setFromTriplets(triplets.begin(), triplets.end());
  }

  if (argDeserialize) {
    if (argMapAtoB || argMapBtoA) {
      DLOG_ERROR(0, "--deserialize and {--mapAtoB, --mapBtoA} are mutually exclusive.");
      std::cout << parser;
      return EXIT_FAILURE;
    }

    try {
      fromSerializedBlob(readBlobFromFile(args::get(argDeserialize)), currentConfig, currentHistory, current_wL, current_smax, sysParam.firstOrderMode, current_nSteps);
    } catch (const std::exception& e) {
      DLOG_ERROR(0, "Failed to deserialize: {}", e.what());
      return EXIT_FAILURE;
    }

    try {
      currentConfig.id = std::stoll(guessNiceNameFromPath(args::get(argDeserialize)));
    } catch (...) {}

    computeEdgePath(currentConfig);

    if (currentConfig.topologyValid)
      computeDescentDirection(currentConfig, currentHistory, current_wL, 0);

  } else {
    if (!argMapAtoB || !argMapBtoA) {
      DLOG_ERROR(0, "Both of --mapAtoB & --mapBtoA must be specified when not using --deserialize.");
      std::cout << parser;
      return EXIT_FAILURE;
    }

    // read solution vector from file
    VectorXsp x(mdataA.nV + mdataB.nV);
    try {
      x.head(mdataA.nV) = readBarycentricCoordinates(mdataA, mdataB, args::get(argMapAtoB)); // passing mdataB because the specified points live on B
      x.tail(mdataB.nV) = readBarycentricCoordinates(mdataB, mdataA, args::get(argMapBtoA));
    } catch (const std::exception& e) {
      DLOG_ERROR(0, "Failed to read barycentric coordinates: {}", e.what());
      return EXIT_FAILURE;
    }

    // For each face point, if its smallest barycentric coordinate is below threshold, snap it to its closest edge
    size_t A_nSnapped = 0;
    size_t B_nSnapped = 0;
    for (size_t i = 0; i < (size_t)x.size(); ++i) {
      if (x[i].type == SurfacePointType::Face && snapFacePointToEdge(x[i])) {
        DLOG_DEBUG(0, "Snapped {}'s vertex {} onto {}'s edge {}",
          i < mdataA.nV ? "A" : "B",
          i < mdataA.nV ? i : (i - mdataA.nV),
          i < mdataA.nV ? "B" : "A",
          x[i].edge);
        ++(i < mdataA.nV ? A_nSnapped : B_nSnapped);
      }
    }
    DLOG_INFO(0, "Snapped {} vertices to edges ({} from A, {} from B)", A_nSnapped + B_nSnapped, A_nSnapped, B_nSnapped);

    uint64_t id = makeConfigurationID();
    createConfiguration(id, x, currentConfig, makeLogger(getResultFileName(ResultContext::UpdateConfiguration, id) + ".log", true));

    current_wL = sysParam.wL_initial;

    computeEdgePath(currentConfig);

    if (currentConfig.topologyValid)
      current_smax = estimateMaxStep(currentConfig, current_wL, 1);  // This internally calls computeDescentDirection

    fillResultInfo(currentConfig, 0., {}, {}, 0, 0.);
    writeResult(currentConfig, currentHistory, current_wL, current_smax, current_nSteps, ResultContext::UpdateConfiguration);
  }

  if (currentConfig.topologyValid) {
    computeOverlayPolygons(currentConfig);
    writeOverlayMesh(currentConfig);
  }

  initialConfig = currentConfig;

  if (!args::get(argInteractive)) {
    // Non-interactive, headless mode

    if (!currentConfig.topologyValid)
      return EXIT_FAILURE;

    optimize(-1);
    return EXIT_SUCCESS;
  }

  // glewInit();

  // Initialize GLFW
  kt84::glfw_util::InitConfig initConfig;
  initConfig.window.title = "Compatible Intrinsic Triangulations (Session ID: " + std::to_string(sessionID) + ")";
  initConfig.callback.framebuffersize = callback_framebuffersize;
  initConfig.callback.mousebutton = callback_mousebutton;
  initConfig.callback.cursorpos = callback_cursorpos;
  initConfig.callback.key = callback_key;
  kt84::glfw_util::LoopControl loopControl = kt84::glfw_util::init(initConfig);
  // idle functions
  loopControl.idle_func[0] = []() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (io->WantCaptureMouse && io->MouseReleased[0]) {
      mdataA.dispList.invalidate();
      mdataB.dispList.invalidate();
    }
  };
  loopControl.idle_func[1] = draw;
  loopControl.idle_func[2] = optimizerMenu;
  loopControl.idle_func[3] = []() {
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
  };

  // OpenGL initial settings
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1,1);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glClearColor(1,1,1,1);

  lightParam[0].enabled = true;
  lightParam[0].ambient = 0.2;
  lightParam[0].diffuse = 1.0;
  lightParam[0].specular = 0.0;
  lightParam[0].position = Vector4f{ -0.40464258, 0.59512496, 0.6943284, 0.0 };
  lightParam[1].enabled = true;
  lightParam[1].ambient = 0.0;
  lightParam[1].diffuse = 0.5;
  lightParam[1].specular = 0.0;
  lightParam[1].position = Vector4f{ 0.8579609, -0.21099381, 0.46838602, 0.0 };
  lightParam[2].enabled = true;
  lightParam[2].ambient = 0.0;
  lightParam[2].diffuse = 0.25;
  lightParam[2].specular = 0.0;
  lightParam[2].position = Vector4f{ -0.5633799, -0.8006998, 0.2036816, 0.0 };
  materialParam.ambient = 0.1;
  materialParam.specular = 0.0;
  materialParam.shininess = 0.0;

  Vector3d A_center = 0.5 * (mdataA.bbox.min() + mdataA.bbox.max());
  Vector3d B_center = 0.5 * (mdataB.bbox.min() + mdataB.bbox.max());
  mdataA.camera.init(A_center + Vector3d(0,0,mdataA.bbox.diagonal().norm() * 1.2), A_center, Vector3d::UnitY());
  mdataB.camera.init(B_center + Vector3d(0,0,mdataB.bbox.diagonal().norm() * 1.2), B_center, Vector3d::UnitY());

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
      DLOG_ERROR(0, "Failed to load texture: {}", args::get(*argTexture));
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
      DLOG_ERROR(0, "Failed to load colormap texture");
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
  io = &ImGui::GetIO();
  io->IniFilename = nullptr;  // Suppress generating imgui.ini
  ImGui_ImplGlfw_InitForOpenGL(loopControl.window, true);
  ImGui_ImplOpenGL2_Init();

  kt84::glfw_util::start_loop(loopControl);

  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  return EXIT_SUCCESS;
}
