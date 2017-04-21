cmake_minimum_required(VERSION 2.8.12)
project(supporter)

add_executable(${PROJECT_NAME}_bin main.cpp
        clipper.hpp
        clipper.cpp
        mesh_slicer.h hash_edge.h
        normalizing_model.h
        scene_organizer.h slice_overhang_detector.h wall_support_generator.h fermat_spirals.h)
        #testing_models_path.h
        #detecting_overhangs.h
        #rendering_tree_support.h
        #scene_organizer.h generating_support.h normalizing_model.h settings.h slice.h clipper.cpp clipper.hpp)
target_include_directories(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME}_bin PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(
        ${PROJECT_NAME}_bin
        ${LIBIGL_LIBRARIES}
        ${LIBIGL_VIEWER_EXTRA_LIBRARIES}
        ${LIBIGL_OPENGL_EXTRA_LIBRARIES}
        ${LIBIGL_OPENGL_GLFW_EXTRA_LIBRARIES}
        ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES})
