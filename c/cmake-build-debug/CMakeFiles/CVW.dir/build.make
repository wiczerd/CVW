# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /home/david/Programs/CLion-2018.3.3/clion-2018.3.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/david/Programs/CLion-2018.3.3/clion-2018.3.3/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/david/workspace/CVW/c

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/david/workspace/CVW/c/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/CVW.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CVW.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CVW.dir/flags.make

CMakeFiles/CVW.dir/main.c.o: CMakeFiles/CVW.dir/flags.make
CMakeFiles/CVW.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/david/workspace/CVW/c/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/CVW.dir/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/CVW.dir/main.c.o   -c /home/david/workspace/CVW/c/main.c

CMakeFiles/CVW.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/CVW.dir/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/david/workspace/CVW/c/main.c > CMakeFiles/CVW.dir/main.c.i

CMakeFiles/CVW.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/CVW.dir/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/david/workspace/CVW/c/main.c -o CMakeFiles/CVW.dir/main.c.s

# Object files for target CVW
CVW_OBJECTS = \
"CMakeFiles/CVW.dir/main.c.o"

# External object files for target CVW
CVW_EXTERNAL_OBJECTS =

CVW: CMakeFiles/CVW.dir/main.c.o
CVW: CMakeFiles/CVW.dir/build.make
CVW: /usr/lib/x86_64-linux-gnu/libgsl.so
CVW: /usr/lib/x86_64-linux-gnu/libgslcblas.so
CVW: CMakeFiles/CVW.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/david/workspace/CVW/c/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable CVW"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CVW.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CVW.dir/build: CVW

.PHONY : CMakeFiles/CVW.dir/build

CMakeFiles/CVW.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CVW.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CVW.dir/clean

CMakeFiles/CVW.dir/depend:
	cd /home/david/workspace/CVW/c/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/david/workspace/CVW/c /home/david/workspace/CVW/c /home/david/workspace/CVW/c/cmake-build-debug /home/david/workspace/CVW/c/cmake-build-debug /home/david/workspace/CVW/c/cmake-build-debug/CMakeFiles/CVW.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CVW.dir/depend

