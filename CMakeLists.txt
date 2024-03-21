
cmake_minimum_required(VERSION 3.2)
project(GDMesh)

set(CMAKE_CXX_STANDARD 17) # Set C++17 standard
set(CMAKE_CXX_STANDARD_REQUIRED True) # Make sure that the specified C++ standard is mandatory

# src path
set(ROOT_SRC ${CMAKE_SOURCE_DIR}/src)
# 3rdparty path deault vc14
set(ROOT_3RDPARTY ${CMAKE_SOURCE_DIR}/3rdparty/MSVC_vc14) 


# Platform
message("$$$$$$$$$ Platform: ${CMAKE_SYSTEM_NAME}")
if(WIN32)
    set(PLATFORM win32)
elseif(APPLE)
    set(PLATFORM macos)
elseif(UNIX)
    set(PLATFORM linux) 
endif()

# all cpp files
file(GLOB_RECURSE ALL_CPP_FILES ${ROOT_SRC}/*.cpp ${ROOT_SRC}/*.cxx ${ROOT_SRC}/*.cc ${ROOT_SRC}/*.c)
# all h files
file(GLOB_RECURSE ALL_H_FILES ${ROOT_SRC}/*.h ${ROOT_SRC}/*.hpp ${ROOT_SRC}/*.hxx)

set(global_include_dir ${ROOT_SRC})
# project include 
include_directories(${global_include_dir})

# 函数: 添加所有子目录
function(add_all_subdirectories dir)
    # 获取当前目录下的所有文件和目录
    file(GLOB children RELATIVE ${dir} ${dir}/*)

    # 遍历所有子项
    foreach(child ${children})
        # 如果子项是目录并且包含CMakeLists.txt，则添加它
        if(IS_DIRECTORY ${dir}/${child} AND EXISTS ${dir}/${child}/CMakeLists.txt)
            add_subdirectory(${child})
        endif()
    endforeach()
endfunction()

# 调用函数添加当前项目目录下的所有子目录
add_all_subdirectories(${CMAKE_CURRENT_SOURCE_DIR})


# Compiler
message("$$$$$$$$$ Compiler: ${CMAKE_CXX_COMPILER_ID}")


# occ include
set(OCC_INCLUDE_DIR ${ROOT_3RDPARTY}/occ-7.6.3/include)
# 3rdparty lib
file(GLOB_RECURSE OCC_LIB_FILES ${ROOT_3RDPARTY}/occ-7.6.3/lib/*.lib)
# Boost include
set(BOOST_INCLUDE_DIR ${ROOT_3RDPARTY}/boost-1.8.4/include/boost-1_84)
# Boost lib
file(GLOB_RECURSE BOOST_LIB_FILES ${ROOT_3RDPARTY}/boost-1.8.3/lib/*.lib)
# eigen include
set(EIGEN_INCLIUDE_DIR ${ROOT_3RDPARTY}/eigen)
# igl include
set(IGL_INCLUDE_DIR ${ROOT_3RDPARTY}/libgl)
# GLFW include
set(GLFW_INCLUDE_DIR ${ROOT_3RDPARTY}/glfw-3.4/include)
# GLFW lib
file(GLOB_RECURSE GLFW_LIB_FILES ${ROOT_3RDPARTY}/glfw-3.4/lib/*.lib)
# glad include
set(GLAD_INCLUDE_DIR ${ROOT_3RDPARTY}/glad/include)
# embree include
set(EMBREE_INCLUDE_DIR ${ROOT_3RDPARTY}/embree-4.3.0/include)
# embree lib
file(GLOB_RECURSE EMBREE_LIB_FILES ${ROOT_3RDPARTY}/embree-4.3.0/lib/*.lib)
# gmsh include
set(GMSH_INCLUDE_DIR ${ROOT_3RDPARTY}/gmsh-4.10.5/include)
# gmsh lib
file(GLOB_RECURSE GMSH_LIB_FILES ${ROOT_3RDPARTY}/gmsh-4.10.5/lib/*.lib)
# OSG include
set(OSG_INCLUDE_DIR ${ROOT_3RDPARTY}/osg-3.6.5/include)
# OSG lib
file(GLOB_RECURSE OSG_LIB_FILES ${ROOT_3RDPARTY}/osg-3.6.5/lib/*.lib)
# googletest include
set(GTEST_INCLUDE_DIR ${ROOT_3RDPARTY}/googletest-1.12.1/include)
# googletest lib
file(GLOB_RECURSE GTEST_LIB_FILES ${ROOT_3RDPARTY}/googletest-1.12.1/lib/*.lib)

message("$$$$$$$$$ CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
# group source files
file(GLOB SUB_DIRS RELATIVE "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/src/*")
foreach(subdir ${SUB_DIRS})
    file(GLOB_RECURSE subdir_cpp_sources "${CMAKE_SOURCE_DIR}/src/${subdir}/*.cpp" )
    file(GLOB_RECURSE subdir_h_sources "${CMAKE_SOURCE_DIR}/src/${subdir}/*.h")
    source_group("${subdir}" FILES ${subdir_cpp_sources} ${subdir_h_sources})
endforeach()


if(OCC_INCLUDE_DIR AND OCC_LIB_FILES)
    include_directories(${OCC_INCLUDE_DIR})
    link_libraries(${OCC_LIB_FILES})
    message("#####  OCC is OK ! ########")
else()
    message("#####  OCC was not found ! ########")
endif()


if(BOOST_INCLUDE_DIR AND BOOST_LIB_FILES)
    include_directories(${BOOST_INCLUDE_DIR})
    link_libraries(${BOOST_LIB_FILES})
    message("#####  BOOST is OK ! ########")
else()
    message("#####  BOOST was not found ! ########")
endif()

if(EIGEN_INCLIUDE_DIR)
    include_directories(${EIGEN_INCLIUDE_DIR})
    message("#####  EIGEN is OK ! ########")
else()
    message("#####  EIGEN was not found ! ########")
endif()

if(IGL_INCLUDE_DIR)
    include_directories(${IGL_INCLUDE_DIR})
    message("#####  IGL is OK ! ########")
else()
    message("#####  IGL was not found ! ########")
endif()

if(GLFW_INCLUDE_DIR AND GLFW_LIB_FILES)
    include_directories(${GLFW_INCLUDE_DIR})
    link_libraries(${GLFW_LIB_FILES})
    message("#####  GLFW is OK ! ########")
else()
    message("#####  GLFW was not found ! ########")
endif()


if(GLAD_INCLUDE_DIR)
    include_directories(${GLAD_INCLUDE_DIR})
    message("#####  GLAD is OK ! ########")
else()
    message("#####  GLAD was not found ! ########")
endif()

if(EMBREE_INCLUDE_DIR AND EMBREE_LIB_FILES)
    include_directories(${EMBREE_INCLUDE_DIR})
    link_libraries(${EMBREE_LIB_FILES})
    message("#####  EMBREE is OK ! ########")
else()
    message("#####  EMBREE was not found ! ########")
endif()


if(OSG_INCLUDE_DIR AND OSG_LIB_FILES)
    include_directories(${OSG_INCLUDE_DIR})
    link_libraries(${OSG_LIB_FILES})
    message("#####  OSG is OK ! ########")
else()
    message("#####  OSG was not found ! ########")
endif()


if(GMSH_INCLUDE_DIR AND GMSH_LIB_FILES)
    include_directories(${GMSH_INCLUDE_DIR})
    link_libraries(${GMSH_LIB_FILES})
    message("#####  GMSH is OK ! ########")
else()
    message("#####  GMSH was not found ! ########")
endif()

if(GTEST_INCLUDE_DIR AND GTEST_LIB_FILES)
    include_directories(${GTEST_INCLUDE_DIR})
    link_libraries(${GTEST_LIB_FILES})
    message("#####  GTEST is OK ! ########")
else()
    message("#####  GTEST was not found ! ########")
endif()



add_executable(${PROJECT_NAME} ${ALL_CPP_FILES} ${ALL_H_FILES})
target_compile_options(${PROJECT_NAME} PRIVATE "/bigobj" "/wd4430")
