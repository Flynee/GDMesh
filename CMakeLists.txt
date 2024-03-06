
cmake_minimum_required(VERSION 3.2)
set(CMAKE_CXX_STANDARD 17)
project(GDMesh)

set(CMAKE_CXX_STANDARD 17) # Set C++11 standard
set(CMAKE_CXX_STANDARD_REQUIRED True) # Make sure that the specified C++ standard is mandatory

# src path
set(ROOT_SRC ${CMAKE_SOURCE_DIR}/src)
# 3rdparty path
set(ROOT_3RDPARTY ${CMAKE_SOURCE_DIR}/3rdparty)


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

add_executable(GDMesh ${ALL_CPP_FILES} ${ALL_H_FILES})

# Compiler
message("$$$$$$$$$ Compiler: ${CMAKE_CXX_COMPILER_ID}")
if(MSVC)
    set(COMPILER vc14)
    target_compile_options(GDMesh PRIVATE "/bigobj" "/wd4430")
    # occ include
    set(OCC_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/occ/include)
    # 3rdparty lib
    file(GLOB_RECURSE OCC_LIB_FILES ${ROOT_3RDPARTY}/vc14/occ/lib/*.lib)
    # Boost include
    set(BOOST_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/Boost/include/boost-1_83)
    # Boost lib
    file(GLOB_RECURSE BOOST_LIB_FILES ${ROOT_3RDPARTY}/vc14/Boost/lib/*.lib)
    # eigen include
    set(EIGEN_INCLIUDE_DIR ${ROOT_3RDPARTY}/vc14/eigen)
    # igl include
    set(IGL_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/libgl)
    # GLFW include
    set(GLFW_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/GLFW/include)
    # GLFW lib
    file(GLOB_RECURSE GLFW_LIB_FILES ${ROOT_3RDPARTY}/vc14/GLFW/lib/*.lib)
    # glad include
    set(GLAD_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/glad/include)
    # embree include
    set(EMBREE_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/embree/include)
    # embree lib
    file(GLOB_RECURSE EMBREE_LIB_FILES ${ROOT_3RDPARTY}/vc14/embree/lib/*.lib)
    # OSG include
    set(OSG_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/osg/include)
    # OSG lib
    file(GLOB_RECURSE OSG_LIB_FILES ${ROOT_3RDPARTY}/vc14/osg/lib/*.lib)

elseif(CMAKE_COMPILER_IS_GNUCXX)
    set(COMPILER gcc)
endif()

message("$$$$$$$$$ CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
# group source files
file(GLOB SUB_DIRS RELATIVE "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/src/*")
foreach(subdir ${SUB_DIRS})
    file(GLOB_RECURSE subdir_cpp_sources "${CMAKE_SOURCE_DIR}/src/${subdir}/*.cpp" )
    file(GLOB_RECURSE subdir_h_sources "${CMAKE_SOURCE_DIR}/src/${subdir}/*.h")
    source_group("${subdir}" FILES ${subdir_cpp_sources} ${subdir_h_sources})
endforeach()


if(OCC_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${OCC_INCLUDE_DIR})
endif()

if(OCC_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${OCC_LIB_FILES})
endif()

if(BOOST_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${BOOST_INCLUDE_DIR})
endif()

if(BOOST_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${BOOST_LIB_FILES})
endif()

if(EIGEN_INCLIUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${EIGEN_INCLIUDE_DIR})
endif()

if(IGL_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${IGL_INCLUDE_DIR})
endif()

if(GLFW_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${GLFW_INCLUDE_DIR})
endif()

if(GLFW_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${GLFW_LIB_FILES})
endif()

if(GLAD_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${GLAD_INCLUDE_DIR})
endif()

if(EMBREE_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${EMBREE_INCLUDE_DIR})
endif()

if(EMBREE_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${EMBREE_LIB_FILES})
endif()

if(OSG_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${OSG_INCLUDE_DIR})
endif()

if(OSG_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${OSG_LIB_FILES})
endif()

# cmake debug mode
set(CMAKE_FIND_DEBUG_MODE 1)
# debug message
if(CMAKE_FIND_DEBUG_MODE)
    if(OCC_INCLUDE_DIR)
        message("$$$$$$$$$ find OCC_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ OCC_INCLUDE_DIR not found")
    endif()

    if(OCC_LIB_FILES)
        message("$$$$$$$$$ find OCC_LIB_FILES")
    else()
        message("$$$$$$$$$ OCC_LIB_FILES not found")
    endif()

    if(BOOST_INCLUDE_DIR)
        message("$$$$$$$$$ find BOOST_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ BOOST_INCLUDE_DIR not found")
    endif()

    if(BOOST_LIB_FILES)
        message("$$$$$$$$$ find BOOST_LIB_FILES")
    else()
        message("$$$$$$$$$ BOOST_LIB_FILES not found")
    endif()

    if(EIGEN_INCLIUDE_DIR)
        message("$$$$$$$$$ find EIGEN_INCLIUDE_DIR")
    else()
        message("$$$$$$$$$ EIGEN_INCLIUDE_DIR not found")
    endif()

    if(IGL_INCLUDE_DIR)
        message("$$$$$$$$$ find IGL_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ IGL_INCLUDE_DIR not found")
    endif()

    if(GLFW_INCLUDE_DIR)
        message("$$$$$$$$$ find GLFW_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ GLFW_INCLUDE_DIR not found")
    endif()

    if(GLFW_LIB_FILES)
        message("$$$$$$$$$ find GLFW_LIB_FILES")
    else()
        message("$$$$$$$$$ GLFW_LIB_FILES not found")
    endif()

    if(GLAD_INCLUDE_DIR)
        message("$$$$$$$$$ find GLAD_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ GLAD_INCLUDE_DIR not found")
    endif()

    if(EMBREE_INCLUDE_DIR)
        message("$$$$$$$$$ find EMBREE_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ EMBREE_INCLUDE_DIR not found")
    endif()

    if(EMBREE_LIB_FILES)
        message("$$$$$$$$$ find EMBREE_LIB_FILES")
    else()
        message("$$$$$$$$$ EMBREE_LIB_FILES not found")
    endif()

    if(OSG_INCLUDE_DIR)
        message("$$$$$$$$$ find OSG_INCLUDE_DIR")
    else()
        message("$$$$$$$$$ OSG_INCLUDE_DIR not found")
    endif()

    if(OSG_LIB_FILES)
        message("$$$$$$$$$ find OSG_LIB_FILES")
    else()
        message("$$$$$$$$$ OSG_LIB_FILES not found")
    endif()

endif()
