add_executable(p3_main_seidel 
    main.cpp
    solver.cpp
)

target_include_directories(p3_main_seidel PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p3_main_seidel PRIVATE
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)