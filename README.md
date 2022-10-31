To make armadillo functions work, we need to

1 - include the library: 
    #include <armadillo>

2 - link blas and lapack in the CMake file:

    bdm_add_executable(${CMAKE_PROJECT_NAME}
    HEADERS ${PROJECT_HEADERS}
    SOURCES ${PROJECT_SOURCES}
    LIBRARIES ${BDM_REQUIRED_LIBRARIES} lapack blas)
