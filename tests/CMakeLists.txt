set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(test_cglass.exe test_cglass.cpp)

if (GRAPH)
  find_package(glfw3 REQUIRED)
  find_package(glew REQUIRED)
  find_package(OpenGL REQUIRED)
set(LIB ${LIB} GLEW::GLEW ${GLFW3_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw)
set(INCLUDES ${INCLUDES} ${GLEW_INCLUDE_DIRS} ${GLFW3_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
else()
  add_definitions(-DNOGRAPH=TRUE)
endif()
add_definitions(-DDEBUG=TRUE)
add_definitions(-DTRACE=TRUE)
add_definitions(-DTESTS=TRUE)

target_link_libraries(test_cglass.exe PRIVATE cglass)

add_test(NAME cglass_test COMMAND test_cglass.exe)
