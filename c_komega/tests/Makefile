# Makefile for C Komega Library Tests
# This Makefile builds and runs the C language tests that correspond
# to the Fortran test programs.

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g -fPIC
LDFLAGS = -lm

# Directories
SRC_DIR = ../src
INCLUDE_DIR = ../include
BUILD_DIR = build
LIB_DIR = ../lib

# Library name
LIBRARY_NAME = libkomega
STATIC_LIB = $(LIB_DIR)/$(LIBRARY_NAME).a
SHARED_LIB = $(LIB_DIR)/$(LIBRARY_NAME).so

# Test executables
TEST_CC_EXEC = $(BUILD_DIR)/test_c_solve_cc
TEST_RC_EXEC = $(BUILD_DIR)/test_c_solve_rc
TEST_CR_EXEC = $(BUILD_DIR)/test_c_solve_cr
TEST_RR_EXEC = $(BUILD_DIR)/test_c_solve_rr

# Test sources
TEST_CC_SOURCE = test_c_solve_cc.c
TEST_RC_SOURCE = test_c_solve_rc.c
TEST_CR_SOURCE = test_c_solve_cr.c
TEST_RR_SOURCE = test_c_solve_rr.c

# Default target
all: $(TEST_CC_EXEC) $(TEST_RC_EXEC) $(TEST_CR_EXEC) $(TEST_RR_EXEC)

# Create directories
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# Build test executables
$(TEST_CC_EXEC): $(TEST_CC_SOURCE) $(STATIC_LIB) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(STATIC_LIB) -o $@ $(LDFLAGS)
	@echo "Complex-Complex test executable created: $@"

$(TEST_RC_EXEC): $(TEST_RC_SOURCE) $(STATIC_LIB) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(STATIC_LIB) -o $@ $(LDFLAGS)
	@echo "Real-Complex test executable created: $@"

$(TEST_CR_EXEC): $(TEST_CR_SOURCE) $(STATIC_LIB) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(STATIC_LIB) -o $@ $(LDFLAGS)
	@echo "Complex-Real test executable created: $@"

$(TEST_RR_EXEC): $(TEST_RR_SOURCE) $(STATIC_LIB) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $< $(STATIC_LIB) -o $@ $(LDFLAGS)
	@echo "Real-Real test executable created: $@"

# Run tests
run-test: $(TEST_CC_EXEC) $(TEST_RC_EXEC) $(TEST_CR_EXEC) $(TEST_RR_EXEC)
	@echo "Running C Komega Library Tests..."
	@echo "=================================="
	@echo ""
	@echo "Running Complex-Complex Test..."
	@echo "==============================="
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_CC_EXEC)
	@echo ""
	@echo "Running Real-Complex Test..."
	@echo "============================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_RC_EXEC)
	@echo ""
	@echo "Running Complex-Real Test..."
	@echo "============================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_CR_EXEC)
	@echo ""
	@echo "Running Real-Real Test..."
	@echo "========================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_RR_EXEC)

# Run individual tests
run-cc: $(TEST_CC_EXEC)
	@echo "Running Complex-Complex Test..."
	@echo "==============================="
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_CC_EXEC)

run-rc: $(TEST_RC_EXEC)
	@echo "Running Real-Complex Test..."
	@echo "============================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_RC_EXEC)

run-cr: $(TEST_CR_EXEC)
	@echo "Running Complex-Real Test..."
	@echo "============================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_CR_EXEC)

run-rr: $(TEST_RR_EXEC)
	@echo "Running Real-Real Test..."
	@echo "========================"
	LD_LIBRARY_PATH=$(LIB_DIR):$$LD_LIBRARY_PATH $(TEST_RR_EXEC)

# Clean build files
clean:
	@echo "Cleaning test build files..."
	@rm -rf $(BUILD_DIR)
	@echo "Test clean complete!"

# Distclean (remove all generated files)
distclean: clean
	@echo "Cleaning all test files..."
	@rm -f *.o *.a *.so
	@echo "Test distclean complete!"

# Show help
help:
	@echo "C Komega Library Tests Makefile"
	@echo "================================"
	@echo "Available targets:"
	@echo "  all        - Build all test executables"
	@echo "  run-test   - Run all tests"
	@echo "  run-cc     - Run Complex-Complex test"
	@echo "  run-rc     - Run Real-Complex test"
	@echo "  run-cr     - Run Complex-Real test"
	@echo "  run-rr     - Run Real-Real test"
	@echo "  clean      - Clean test build files"
	@echo "  distclean  - Clean everything"
	@echo "  help       - Show this help"

# Phony targets
.PHONY: all run-test run-cc run-rc run-cr run-rr clean distclean help

# Dependencies
$(TEST_CC_EXEC): $(STATIC_LIB)
$(TEST_RC_EXEC): $(STATIC_LIB)
$(TEST_CR_EXEC): $(STATIC_LIB)
$(TEST_RR_EXEC): $(STATIC_LIB)