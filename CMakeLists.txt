
cmake_minimum_required(VERSION 3.2)
set(CMAKE_CXX_STANDARD 11)
project(GDMesh)


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

# Compiler
message("$$$$$$$$$ Compiler: ${CMAKE_CXX_COMPILER_ID}")
if(MSVC)
    set(COMPILER vc14)
    # occ include
    set(OCC_INCLUDE_DIR ${ROOT_3RDPARTY}/vc14/occ/include)
    # 3rdparty lib
    file(GLOB_RECURSE OCC_LIB_FILES ${ROOT_3RDPARTY}/vc14/occ/lib/*.lib)

elseif(CMAKE_COMPILER_IS_GNUCXX)
    set(COMPILER gcc)
endif()




add_executable(GDMesh ${ROOT_SRC}/main.cpp)

if(OCC_INCLUDE_DIR)
    target_include_directories(GDMesh PUBLIC ${OCC_INCLUDE_DIR})
endif()

if(OCC_LIB_FILES)
    target_link_libraries(GDMesh PUBLIC ${OCC_LIB_FILES})
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

endif()
