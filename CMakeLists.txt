CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

PROJECT(hfv1-cpp-project)

ADD_EXECUTABLE(hfv1-cpp)

TARGET_SOURCES(
    hfv1-cpp
    PRIVATE
    MAIN.cpp
    input/set_boundary_conditions.cpp
    input/set_simulation_parameters.cpp
    input/set_solver_parameters.cpp
    mra/decode1.cpp
    mra/decode2.cpp
    mra/decoding.cpp
    mra/encode_detail.cpp
    mra/encode_scale.cpp
    mra/encoding.cpp
    mra/extra_significance.cpp
    mra/get_max_scale_coefficients.cpp
    mra/load_fine_scale_coefficients.cpp
    mra/regularisation.cpp
    mra/thresholding.cpp
    mra/traverse_tree_decode.cpp
    mra/get_modal_values.cpp
    mra/get_nodal_values.cpp
    operator/add_ghost_cells.cpp
    operator/bed_data_c_property.cpp
    operator/bed_data_triangle.cpp
    operator/fluxHLL.cpp
    operator/friction_implicit.cpp
    operator/friction_update.cpp
    operator/fv1_operator.cpp
    operator/get_bar_values.cpp
    operator/get_dt_CFL.cpp
    operator/get_face_values.cpp
    operator/get_positivity_preserving_nodes.cpp
    operator/get_wet_dry_cells.cpp
    operator/h_init_c_property.cpp
    operator/h_init_overtop.cpp
    output/write_solution_to_file.cpp
)

TARGET_INCLUDE_DIRECTORIES(
    hfv1-cpp
    PRIVATE
    classes
	input
	mra
	operator
	output
)
