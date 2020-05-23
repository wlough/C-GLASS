# Add the cglass executable
add_executable(cglass.exe "main.cpp")
target_link_libraries(cglass.exe cglass)
install(TARGETS cglass.exe DESTINATION bin)
