cmake_minimum_required(VERSION 3.10.0)

project(compatible_intrinsic_triangulations_demo)

set(CMAKE_CXX_STANDARD 20)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

SET(EXTRA_LIBS "")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-Wall -Wextra -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field -Wno-logical-op-parentheses -Wno-unknown-pragmas -Wno-deprecated -Wno-reorder-ctor -Wno-overloaded-virtual")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=0 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")

  if (NOT APPLE)
    SET(EXTRA_LIBS GL;GLU)
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  # add_definitions(/D "__TBB_NO_IMPLICIT_LINKAGE")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)

  SET(EXTRA_LIBS opengl32;glu32)

else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

add_definitions(-DIMGUIZMO_USE_ONLY_ROT)

# == Deps

# Needed to prevent the use of system-installed TBB via geometry-central/cmake/FindSuiteSparse.cmake
set(SUITESPARSE OFF CACHE BOOL " " FORCE)
set(GC_EIGEN_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen CACHE PATH "Eigen path")
add_subdirectory(deps/geometry-central)

# GLEW
add_subdirectory(deps/glew-2.2.0/build/cmake)
include_directories(deps/glew-2.2.0/include)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/spdlog/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/imGuIZMO.quat/imGuIZMO.quat/")
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
include_directories("/usr/local/include/")
include_directories("/opt/homebrew/include")

# Boost
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

# bin2c borrowed from https://github.com/wjakob/nanogui
file(GLOB resources "${CMAKE_CURRENT_SOURCE_DIR}/resources/*")
string (REGEX REPLACE "([^\\]|^);" "\\1," resources_string "${resources}")
string (REGEX REPLACE "[\\](.)" "\\1" resources_string "${resources_string}")
set(bin2c_cmdline
  -DOUTPUT_C=cit_resources.cpp
  -DOUTPUT_H=cit_resources.h
  "-DINPUT_FILES=${resources_string}"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/bin2c.cmake")
add_custom_command(
  OUTPUT cit_resources.cpp cit_resources.h
  COMMAND ${CMAKE_COMMAND} ARGS ${bin2c_cmdline}
  DEPENDS ${resources}
  COMMENT "Running bin2c"
  PRE_BUILD VERBATIM)

# GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/deps/glfw")

# ImGUI (CMakeLists.txt not provided)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/backends")
add_library(imgui
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/imgui.cpp
  # ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/imgui_demo.cpp       # Useful for testing via ImGui::ShowDemoWindow
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/imgui_draw.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/imgui_tables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/imgui_widgets.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_glfw.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imgui/backends/imgui_impl_opengl2.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/imGuIZMO.quat/imGuIZMO.quat/imGuIZMOquat.cpp)
target_link_libraries(imgui glfw)

# == Build our project stuff

find_package(Git QUIET)
set(VERSIONTAG "unknown")
if (GIT_FOUND OR Git_FOUND)
  execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD RESULT_VARIABLE RET OUTPUT_VARIABLE COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (NOT RET)
    set(VERSIONTAG "${COMMIT}")
  endif()
endif()
message(STATUS "Version tag: ${VERSIONTAG}")
configure_file("src/version.cpp.in" "version.cpp")

set(CIT_MAIN_SRCS
  cit_resources.cpp
  version.cpp
  src/core.cpp
  src/event.cpp
  src/global.cpp
  src/main.cpp
  src/optimizer.cpp
  src/render.cpp
  src/serialization.cpp
  src/utility.cpp
)

add_executable(cit_demo "${CIT_MAIN_SRCS}")
target_link_libraries(cit_demo geometry-central glew imgui ${EXTRA_LIBS})

add_executable(bc_extractor src/bc_extractor.cpp)
target_link_libraries(bc_extractor geometry-central)

add_executable(om_viewer
  cit_resources.cpp
  src/om_viewer.cpp
)
target_link_libraries(om_viewer geometry-central glew imgui ${EXTRA_LIBS})

option(BUILD_TESTS "Build tests" OFF)
if (BUILD_TESTS)
  enable_testing()

  ### Download and unpack googletest at configure time
  # (from the googletest docs)
  if(${CMAKE_VERSION} VERSION_LESS 3.14)
    include(FetchContentLocal)
  else()
    include(FetchContent)
  endif()

  FetchContent_Declare(
      googletest
      GIT_REPOSITORY    https://github.com/google/googletest.git
      GIT_TAG           release-1.10.0
      SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
      BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
  )

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  set_property(TARGET gtest PROPERTY CXX_STANDARD 20)

  set(CIT_TEST_SRCS
    test/test_main.cpp
    test/test_angle.cpp
    test/test_energy.cpp
  )

  add_executable(cit_test "${CIT_TEST_SRCS}")
  target_link_libraries(cit_test geometry-central gtest_main)
endif()
