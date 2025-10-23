#!/bin/bash

# C Komega Library Test Runner
# This script runs all the C language tests that correspond to the Fortran tests.

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to run a test
run_test() {
    local test_name="$1"
    local test_exec="$2"
    local test_desc="$3"
    
    print_status "Running $test_desc..."
    echo "=================================="
    
    if [ -f "$test_exec" ]; then
        if LD_LIBRARY_PATH="../lib:$LD_LIBRARY_PATH" "$test_exec"; then
            print_success "$test_desc completed successfully"
        else
            print_error "$test_desc failed with exit code $?"
            return 1
        fi
    else
        print_error "$test_desc executable not found: $test_exec"
        return 1
    fi
    
    echo ""
}

# Main function
main() {
    echo "C Komega Library Test Runner"
    echo "============================"
    echo ""
    
    # Check if we're in the right directory
    if [ ! -f "Makefile" ]; then
        print_error "Makefile not found. Please run this script from the tests directory."
        exit 1
    fi
    
    # Build all tests
    print_status "Building all test executables..."
    if make all; then
        print_success "All test executables built successfully"
    else
        print_error "Failed to build test executables"
        exit 1
    fi
    
    echo ""
    
    # Run all tests
    print_status "Running all C Komega Library tests..."
    echo ""
    
    local failed_tests=0
    
    # Run Complex-Complex test
    if ! run_test "Complex-Complex" "build/test_c_solve_cc" "Complex-Complex Solver Test"; then
        ((failed_tests++))
    fi
    
    # Run Real-Complex test
    if ! run_test "Real-Complex" "build/test_c_solve_rc" "Real-Complex Solver Test"; then
        ((failed_tests++))
    fi
    
    # Run Complex-Real test
    if ! run_test "Complex-Real" "build/test_c_solve_cr" "Complex-Real Solver Test"; then
        ((failed_tests++))
    fi
    
    # Run Real-Real test
    if ! run_test "Real-Real" "build/test_c_solve_rr" "Real-Real Solver Test"; then
        ((failed_tests++))
    fi
    
    echo "=================================="
    
    # Summary
    if [ $failed_tests -eq 0 ]; then
        print_success "All tests completed successfully!"
        echo ""
        print_status "Test Summary:"
        echo "  Complex-Complex Test: PASSED"
        echo "  Real-Complex Test:    PASSED"
        echo "  Complex-Real Test:    PASSED"
        echo "  Real-Real Test:       PASSED"
        echo ""
        print_success "C Komega Library tests are working correctly!"
    else
        print_error "$failed_tests test(s) failed!"
        echo ""
        print_status "Test Summary:"
        echo "  Complex-Complex Test: $([ -f "build/test_c_solve_cc" ] && echo "PASSED" || echo "FAILED")"
        echo "  Real-Complex Test:    $([ -f "build/test_c_solve_rc" ] && echo "PASSED" || echo "FAILED")"
        echo "  Complex-Real Test:    $([ -f "build/test_c_solve_cr" ] && echo "PASSED" || echo "FAILED")"
        echo "  Real-Real Test:       $([ -f "build/test_c_solve_rr" ] && echo "PASSED" || echo "FAILED")"
        echo ""
        print_error "Some tests failed. Please check the output above for details."
        exit 1
    fi
}

# Run main function
main "$@"