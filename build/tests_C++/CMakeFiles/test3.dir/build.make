# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/oneapi/topLevel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/oneapi/topLevel/build

# Include any dependencies generated for this target.
include tests_C++/CMakeFiles/test3.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests_C++/CMakeFiles/test3.dir/compiler_depend.make

# Include the progress variables for this target.
include tests_C++/CMakeFiles/test3.dir/progress.make

# Include the compile flags for this target's objects.
include tests_C++/CMakeFiles/test3.dir/flags.make

tests_C++/CMakeFiles/test3.dir/test3.cpp.o: tests_C++/CMakeFiles/test3.dir/flags.make
tests_C++/CMakeFiles/test3.dir/test3.cpp.o: ../tests_C++/test3.cpp
tests_C++/CMakeFiles/test3.dir/test3.cpp.o: tests_C++/CMakeFiles/test3.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/oneapi/topLevel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests_C++/CMakeFiles/test3.dir/test3.cpp.o"
	cd /home/oneapi/topLevel/build/tests_C++ && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests_C++/CMakeFiles/test3.dir/test3.cpp.o -MF CMakeFiles/test3.dir/test3.cpp.o.d -o CMakeFiles/test3.dir/test3.cpp.o -c /home/oneapi/topLevel/tests_C++/test3.cpp

tests_C++/CMakeFiles/test3.dir/test3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test3.dir/test3.cpp.i"
	cd /home/oneapi/topLevel/build/tests_C++ && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/oneapi/topLevel/tests_C++/test3.cpp > CMakeFiles/test3.dir/test3.cpp.i

tests_C++/CMakeFiles/test3.dir/test3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test3.dir/test3.cpp.s"
	cd /home/oneapi/topLevel/build/tests_C++ && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/oneapi/topLevel/tests_C++/test3.cpp -o CMakeFiles/test3.dir/test3.cpp.s

# Object files for target test3
test3_OBJECTS = \
"CMakeFiles/test3.dir/test3.cpp.o"

# External object files for target test3
test3_EXTERNAL_OBJECTS =

tests_C++/test3: tests_C++/CMakeFiles/test3.dir/test3.cpp.o
tests_C++/test3: tests_C++/CMakeFiles/test3.dir/build.make
tests_C++/test3: mole_C++/libmole_C++.a
tests_C++/test3: /usr/lib/x86_64-linux-gnu/libarmadillo.so
tests_C++/test3: /opt/OpenBLAS/lib/libopenblas.so
tests_C++/test3: /usr/local/lib/libsuperlu.a
tests_C++/test3: tests_C++/CMakeFiles/test3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/oneapi/topLevel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test3"
	cd /home/oneapi/topLevel/build/tests_C++ && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests_C++/CMakeFiles/test3.dir/build: tests_C++/test3
.PHONY : tests_C++/CMakeFiles/test3.dir/build

tests_C++/CMakeFiles/test3.dir/clean:
	cd /home/oneapi/topLevel/build/tests_C++ && $(CMAKE_COMMAND) -P CMakeFiles/test3.dir/cmake_clean.cmake
.PHONY : tests_C++/CMakeFiles/test3.dir/clean

tests_C++/CMakeFiles/test3.dir/depend:
	cd /home/oneapi/topLevel/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/oneapi/topLevel /home/oneapi/topLevel/tests_C++ /home/oneapi/topLevel/build /home/oneapi/topLevel/build/tests_C++ /home/oneapi/topLevel/build/tests_C++/CMakeFiles/test3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests_C++/CMakeFiles/test3.dir/depend

