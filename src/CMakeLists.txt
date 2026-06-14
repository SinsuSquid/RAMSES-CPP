target_sources(${LIB_NAME} PRIVATE
    Initializer.cpp
)

add_subdirectory(core)
add_subdirectory(solvers)
add_subdirectory(particles)
add_subdirectory(io)
