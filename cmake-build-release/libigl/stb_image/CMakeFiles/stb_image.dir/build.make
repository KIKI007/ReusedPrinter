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
include libigl/stb_image/CMakeFiles/stb_image.dir/depend.make

# Include the progress variables for this target.
include libigl/stb_image/CMakeFiles/stb_image.dir/progress.make

# Include the compile flags for this target's objects.
include libigl/stb_image/CMakeFiles/stb_image.dir/flags.make

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o: libigl/stb_image/CMakeFiles/stb_image.dir/flags.make
libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o: ../external/stb_image/stb_image.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stb_image.dir/stb_image.cpp.o -c "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/external/stb_image/stb_image.cpp"

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stb_image.dir/stb_image.cpp.i"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/external/stb_image/stb_image.cpp" > CMakeFiles/stb_image.dir/stb_image.cpp.i

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stb_image.dir/stb_image.cpp.s"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/external/stb_image/stb_image.cpp" -o CMakeFiles/stb_image.dir/stb_image.cpp.s

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.requires:

.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.requires

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.provides: libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.requires
	$(MAKE) -f libigl/stb_image/CMakeFiles/stb_image.dir/build.make libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.provides.build
.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.provides

libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.provides.build: libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o


# Object files for target stb_image
stb_image_OBJECTS = \
"CMakeFiles/stb_image.dir/stb_image.cpp.o"

# External object files for target stb_image
stb_image_EXTERNAL_OBJECTS =

libigl/stb_image/libstb_image.a: libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o
libigl/stb_image/libstb_image.a: libigl/stb_image/CMakeFiles/stb_image.dir/build.make
libigl/stb_image/libstb_image.a: libigl/stb_image/CMakeFiles/stb_image.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libstb_image.a"
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && $(CMAKE_COMMAND) -P CMakeFiles/stb_image.dir/cmake_clean_target.cmake
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/stb_image.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libigl/stb_image/CMakeFiles/stb_image.dir/build: libigl/stb_image/libstb_image.a

.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/build

libigl/stb_image/CMakeFiles/stb_image.dir/requires: libigl/stb_image/CMakeFiles/stb_image.dir/stb_image.cpp.o.requires

.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/requires

libigl/stb_image/CMakeFiles/stb_image.dir/clean:
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" && $(CMAKE_COMMAND) -P CMakeFiles/stb_image.dir/cmake_clean.cmake
.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/clean

libigl/stb_image/CMakeFiles/stb_image.dir/depend:
	cd "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/external/stb_image" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image" "/Users/wangziqi/Desktop/USC Spring/Code/Supporter/cmake-build-release/libigl/stb_image/CMakeFiles/stb_image.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : libigl/stb_image/CMakeFiles/stb_image.dir/depend
