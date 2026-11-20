cmake_minimum_required(VERSION 3.20)
# Builds all obdeect components; one ctest run covers all of them. Each component also
# builds on its own from src/<component> (wheels are built per component).
project(obdeect_all LANGUAGES C CXX)

enable_testing()

# The stream component first: the other components use its target.
add_subdirectory(src/stream)
add_subdirectory(src/propagation)
add_subdirectory(src/sources)
add_subdirectory(src/ray_tracing)
