# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/wangziqi/Desktop/USC Spring/Code/Supporter"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release"

# Include any dependencies generated for this target.
include libigl/CMakeFiles/igl_opengl_glfw.dir/depend.make

# Include the progress variables for this target.
include libigl/CMakeFiles/igl_opengl_glfw.dir/progress.make

# Include the compile flags for this target's objects.
include libigl/CMakeFiles/igl_opengl_glfw.dir/flags.make

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o: libigl/CMakeFiles/igl_opengl_glfw.dir/flags.make
libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o: ../include/igl/opengl/glfw/map_texture.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o -c "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/include/igl/opengl/glfw/map_texture.cpp"

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.i"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/include/igl/opengl/glfw/map_texture.cpp" > CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.i

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.s"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/include/igl/opengl/glfw/map_texture.cpp" -o CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.s

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.requires:

.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.requires

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.provides: libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.requires
	$(MAKE) -f libigl/CMakeFiles/igl_opengl_glfw.dir/build.make libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.provides.build
.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.provides

libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.provides.build: libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o


# Object files for target igl_opengl_glfw
igl_opengl_glfw_OBJECTS = \
"CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o"

# External object files for target igl_opengl_glfw
igl_opengl_glfw_EXTERNAL_OBJECTS =

libigl/libigl_opengl_glfw.a: libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o
libigl/libigl_opengl_glfw.a: libigl/CMakeFiles/igl_opengl_glfw.dir/build.make
libigl/libigl_opengl_glfw.a: libigl/CMakeFiles/igl_opengl_glfw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libigl_opengl_glfw.a"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && $(CMAKE_COMMAND) -P CMakeFiles/igl_opengl_glfw.dir/cmake_clean_target.cmake
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/igl_opengl_glfw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libigl/CMakeFiles/igl_opengl_glfw.dir/build: libigl/libigl_opengl_glfw.a

.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/build

libigl/CMakeFiles/igl_opengl_glfw.dir/requires: libigl/CMakeFiles/igl_opengl_glfw.dir/__/__/include/igl/opengl/glfw/map_texture.cpp.o.requires

.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/requires

libigl/CMakeFiles/igl_opengl_glfw.dir/clean:
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" && $(CMAKE_COMMAND) -P CMakeFiles/igl_opengl_glfw.dir/cmake_clean.cmake
.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/clean

libigl/CMakeFiles/igl_opengl_glfw.dir/depend:
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/shared/cmake" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/CMakeFiles/igl_opengl_glfw.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : libigl/CMakeFiles/igl_opengl_glfw.dir/depend

