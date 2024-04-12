# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_SOURCE_DIR = /home/airribarra/ring_tumbes/ring-basic-p-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/airribarra/ring_tumbes/ring-basic-p-master

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/airribarra/ring_tumbes/ring-basic-p-master/CMakeFiles /home/airribarra/ring_tumbes/ring-basic-p-master//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/airribarra/ring_tumbes/ring-basic-p-master/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named build-index

# Build rule for target.
build-index: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 build-index
.PHONY : build-index

# fast build rule for target.
build-index/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index.dir/build.make CMakeFiles/build-index.dir/build
.PHONY : build-index/fast

#=============================================================================
# Target rules for targets named query-index

# Build rule for target.
query-index: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 query-index
.PHONY : query-index

# fast build rule for target.
query-index/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index.dir/build.make CMakeFiles/query-index.dir/build
.PHONY : query-index/fast

#=============================================================================
# Target rules for targets named query-index-fixed

# Build rule for target.
query-index-fixed: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 query-index-fixed
.PHONY : query-index-fixed

# fast build rule for target.
query-index-fixed/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-fixed.dir/build.make CMakeFiles/query-index-fixed.dir/build
.PHONY : query-index-fixed/fast

#=============================================================================
# Target rules for targets named query-index-random

# Build rule for target.
query-index-random: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 query-index-random
.PHONY : query-index-random

# fast build rule for target.
query-index-random/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-random.dir/build.make CMakeFiles/query-index-random.dir/build
.PHONY : query-index-random/fast

#=============================================================================
# Target rules for targets named build-index-uring

# Build rule for target.
build-index-uring: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 build-index-uring
.PHONY : build-index-uring

# fast build rule for target.
build-index-uring/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index-uring.dir/build.make CMakeFiles/build-index-uring.dir/build
.PHONY : build-index-uring/fast

#=============================================================================
# Target rules for targets named query-index-uring

# Build rule for target.
query-index-uring: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 query-index-uring
.PHONY : query-index-uring

# fast build rule for target.
query-index-uring/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring.dir/build.make CMakeFiles/query-index-uring.dir/build
.PHONY : query-index-uring/fast

#=============================================================================
# Target rules for targets named query-index-uring-fixed

# Build rule for target.
query-index-uring-fixed: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 query-index-uring-fixed
.PHONY : query-index-uring-fixed

# fast build rule for target.
query-index-uring-fixed/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring-fixed.dir/build.make CMakeFiles/query-index-uring-fixed.dir/build
.PHONY : query-index-uring-fixed/fast

#=============================================================================
# Target rules for targets named check

# Build rule for target.
check: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 check
.PHONY : check

# fast build rule for target.
check/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/check.dir/build.make CMakeFiles/check.dir/build
.PHONY : check/fast

#=============================================================================
# Target rules for targets named proba

# Build rule for target.
proba: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 proba
.PHONY : proba

# fast build rule for target.
proba/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/proba.dir/build.make CMakeFiles/proba.dir/build
.PHONY : proba/fast

src/build-index-uring.o: src/build-index-uring.cpp.o
.PHONY : src/build-index-uring.o

# target to build an object file
src/build-index-uring.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index-uring.dir/build.make CMakeFiles/build-index-uring.dir/src/build-index-uring.cpp.o
.PHONY : src/build-index-uring.cpp.o

src/build-index-uring.i: src/build-index-uring.cpp.i
.PHONY : src/build-index-uring.i

# target to preprocess a source file
src/build-index-uring.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index-uring.dir/build.make CMakeFiles/build-index-uring.dir/src/build-index-uring.cpp.i
.PHONY : src/build-index-uring.cpp.i

src/build-index-uring.s: src/build-index-uring.cpp.s
.PHONY : src/build-index-uring.s

# target to generate assembly for a file
src/build-index-uring.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index-uring.dir/build.make CMakeFiles/build-index-uring.dir/src/build-index-uring.cpp.s
.PHONY : src/build-index-uring.cpp.s

src/build-index.o: src/build-index.cpp.o
.PHONY : src/build-index.o

# target to build an object file
src/build-index.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index.dir/build.make CMakeFiles/build-index.dir/src/build-index.cpp.o
.PHONY : src/build-index.cpp.o

src/build-index.i: src/build-index.cpp.i
.PHONY : src/build-index.i

# target to preprocess a source file
src/build-index.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index.dir/build.make CMakeFiles/build-index.dir/src/build-index.cpp.i
.PHONY : src/build-index.cpp.i

src/build-index.s: src/build-index.cpp.s
.PHONY : src/build-index.s

# target to generate assembly for a file
src/build-index.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/build-index.dir/build.make CMakeFiles/build-index.dir/src/build-index.cpp.s
.PHONY : src/build-index.cpp.s

src/check.o: src/check.cpp.o
.PHONY : src/check.o

# target to build an object file
src/check.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/check.dir/build.make CMakeFiles/check.dir/src/check.cpp.o
.PHONY : src/check.cpp.o

src/check.i: src/check.cpp.i
.PHONY : src/check.i

# target to preprocess a source file
src/check.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/check.dir/build.make CMakeFiles/check.dir/src/check.cpp.i
.PHONY : src/check.cpp.i

src/check.s: src/check.cpp.s
.PHONY : src/check.s

# target to generate assembly for a file
src/check.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/check.dir/build.make CMakeFiles/check.dir/src/check.cpp.s
.PHONY : src/check.cpp.s

src/proba.o: src/proba.cpp.o
.PHONY : src/proba.o

# target to build an object file
src/proba.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/proba.dir/build.make CMakeFiles/proba.dir/src/proba.cpp.o
.PHONY : src/proba.cpp.o

src/proba.i: src/proba.cpp.i
.PHONY : src/proba.i

# target to preprocess a source file
src/proba.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/proba.dir/build.make CMakeFiles/proba.dir/src/proba.cpp.i
.PHONY : src/proba.cpp.i

src/proba.s: src/proba.cpp.s
.PHONY : src/proba.s

# target to generate assembly for a file
src/proba.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/proba.dir/build.make CMakeFiles/proba.dir/src/proba.cpp.s
.PHONY : src/proba.cpp.s

src/query-index-uring.o: src/query-index-uring.cpp.o
.PHONY : src/query-index-uring.o

# target to build an object file
src/query-index-uring.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring.dir/build.make CMakeFiles/query-index-uring.dir/src/query-index-uring.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring-fixed.dir/build.make CMakeFiles/query-index-uring-fixed.dir/src/query-index-uring.cpp.o
.PHONY : src/query-index-uring.cpp.o

src/query-index-uring.i: src/query-index-uring.cpp.i
.PHONY : src/query-index-uring.i

# target to preprocess a source file
src/query-index-uring.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring.dir/build.make CMakeFiles/query-index-uring.dir/src/query-index-uring.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring-fixed.dir/build.make CMakeFiles/query-index-uring-fixed.dir/src/query-index-uring.cpp.i
.PHONY : src/query-index-uring.cpp.i

src/query-index-uring.s: src/query-index-uring.cpp.s
.PHONY : src/query-index-uring.s

# target to generate assembly for a file
src/query-index-uring.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring.dir/build.make CMakeFiles/query-index-uring.dir/src/query-index-uring.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-uring-fixed.dir/build.make CMakeFiles/query-index-uring-fixed.dir/src/query-index-uring.cpp.s
.PHONY : src/query-index-uring.cpp.s

src/query-index.o: src/query-index.cpp.o
.PHONY : src/query-index.o

# target to build an object file
src/query-index.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index.dir/build.make CMakeFiles/query-index.dir/src/query-index.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-fixed.dir/build.make CMakeFiles/query-index-fixed.dir/src/query-index.cpp.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-random.dir/build.make CMakeFiles/query-index-random.dir/src/query-index.cpp.o
.PHONY : src/query-index.cpp.o

src/query-index.i: src/query-index.cpp.i
.PHONY : src/query-index.i

# target to preprocess a source file
src/query-index.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index.dir/build.make CMakeFiles/query-index.dir/src/query-index.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-fixed.dir/build.make CMakeFiles/query-index-fixed.dir/src/query-index.cpp.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-random.dir/build.make CMakeFiles/query-index-random.dir/src/query-index.cpp.i
.PHONY : src/query-index.cpp.i

src/query-index.s: src/query-index.cpp.s
.PHONY : src/query-index.s

# target to generate assembly for a file
src/query-index.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index.dir/build.make CMakeFiles/query-index.dir/src/query-index.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-fixed.dir/build.make CMakeFiles/query-index-fixed.dir/src/query-index.cpp.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/query-index-random.dir/build.make CMakeFiles/query-index-random.dir/src/query-index.cpp.s
.PHONY : src/query-index.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... build-index"
	@echo "... build-index-uring"
	@echo "... check"
	@echo "... proba"
	@echo "... query-index"
	@echo "... query-index-fixed"
	@echo "... query-index-random"
	@echo "... query-index-uring"
	@echo "... query-index-uring-fixed"
	@echo "... src/build-index-uring.o"
	@echo "... src/build-index-uring.i"
	@echo "... src/build-index-uring.s"
	@echo "... src/build-index.o"
	@echo "... src/build-index.i"
	@echo "... src/build-index.s"
	@echo "... src/check.o"
	@echo "... src/check.i"
	@echo "... src/check.s"
	@echo "... src/proba.o"
	@echo "... src/proba.i"
	@echo "... src/proba.s"
	@echo "... src/query-index-uring.o"
	@echo "... src/query-index-uring.i"
	@echo "... src/query-index-uring.s"
	@echo "... src/query-index.o"
	@echo "... src/query-index.i"
	@echo "... src/query-index.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

