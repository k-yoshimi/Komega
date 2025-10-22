#!/bin/bash

# Test runner script for Komega library
# This script runs all available tests and verifies their results

set -e  # Exit on any error

echo "=========================================="
echo "Running Komega Library Tests"
echo "=========================================="

# Function to run a test and check its success
run_test() {
    local test_name="$1"
    local input_file="$2"
    local expected_output="$3"
    
    echo "Running $test_name..."
    
    if [ ! -f "$input_file" ]; then
        echo "ERROR: Input file $input_file not found"
        return 1
    fi
    
    if [ ! -f "./$test_name" ]; then
        echo "ERROR: Test executable $test_name not found"
        return 1
    fi
    
    # Run the test
    if ./$test_name < "$input_file" > "${test_name}.log" 2>&1; then
        echo "✓ $test_name completed successfully"
        
        # Check for convergence indicators in the log
        if grep -q "Converged\|converged\|SUCCESS" "${test_name}.log"; then
            echo "✓ $test_name converged properly"
        else
            echo "⚠ $test_name may not have converged (check log)"
        fi
        
        # Check for error indicators
        if grep -q "ERROR\|error\|FAILED\|failed" "${test_name}.log"; then
            echo "⚠ $test_name may have encountered errors (check log)"
        fi
        
        return 0
    else
        echo "✗ $test_name failed"
        echo "Last few lines of output:"
        tail -10 "${test_name}.log"
        return 1
    fi
}

# Function to check if a test produces expected output files
check_output_files() {
    local test_name="$1"
    local expected_files="$2"
    
    echo "Checking output files for $test_name..."
    
    for file in $expected_files; do
        if [ -f "$file" ]; then
            echo "✓ Found $file"
        else
            echo "⚠ Missing expected file: $file"
        fi
    done
}

# Build all test programs
echo "Building test programs..."
make clean
make -j$(nproc)

# Verify all executables were built
required_tests=("solve_rr.x" "solve_rc.x" "solve_cr.x" "solve_cc.x" "solve_rc_qmr1.x" "solve_rc_qmr2.x")
for test in "${required_tests[@]}"; do
    if [ ! -f "$test" ]; then
        echo "ERROR: Failed to build $test"
        exit 1
    fi
done

echo "All test programs built successfully"
echo ""

# Run individual tests
test_results=()

echo "Running solver tests..."
echo "----------------------"

# Test 1: Real-Real solver
if run_test "solve_rr.x" "real_freq.in" "converged"; then
    test_results+=("solve_rr.x: PASS")
else
    test_results+=("solve_rr.x: FAIL")
fi

# Test 2: Real-Complex solver  
if run_test "solve_rc.x" "real_freq.in" "converged"; then
    test_results+=("solve_rc.x: PASS")
else
    test_results+=("solve_rc.x: FAIL")
fi

# Test 3: Complex-Real solver
if run_test "solve_cr.x" "real_freq.in" "converged"; then
    test_results+=("solve_cr.x: PASS")
else
    test_results+=("solve_cr.x: FAIL")
fi

# Test 4: Complex-Complex solver
if run_test "solve_cc.x" "complex_freq.in" "converged"; then
    test_results+=("solve_cc.x: PASS")
else
    test_results+=("solve_cc.x: FAIL")
fi

# Test 5: QMR1 solver
if run_test "solve_rc_qmr1.x" "real_freq.in" "converged"; then
    test_results+=("solve_rc_qmr1.x: PASS")
else
    test_results+=("solve_rc_qmr1.x: FAIL")
fi

# Test 6: QMR2 solver
if run_test "solve_rc_qmr2.x" "real_freq.in" "converged"; then
    test_results+=("solve_rc_qmr2.x: PASS")
else
    test_results+=("solve_rc_qmr2.x: FAIL")
fi

echo ""
echo "=========================================="
echo "Test Results Summary"
echo "=========================================="

# Count results
passed=0
failed=0

for result in "${test_results[@]}"; do
    echo "$result"
    if [[ "$result" == *"PASS" ]]; then
        ((passed++))
    else
        ((failed++))
    fi
done

echo ""
echo "Total: $((passed + failed)) tests"
echo "Passed: $passed"
echo "Failed: $failed"

# Check for restart files
echo ""
echo "Checking for restart files..."
restart_files=("restart.dat")
for file in "${restart_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ Found $file"
    else
        echo "⚠ Missing restart file: $file"
    fi
done

# Performance summary
echo ""
echo "Performance Summary:"
echo "-------------------"
for test in "${required_tests[@]}"; do
    if [ -f "${test}.log" ]; then
        echo "$test:"
        # Extract timing information if available
        if grep -q "time\|Time\|TIME" "${test}.log"; then
            grep -i "time" "${test}.log" | head -3
        fi
        # Extract iteration count if available
        if grep -q "iteration\|iter\|Iteration" "${test}.log"; then
            grep -i "iteration" "${test}.log" | tail -1
        fi
        echo ""
    fi
done

# Final result
if [ $failed -eq 0 ]; then
    echo "🎉 All tests passed!"
    exit 0
else
    echo "❌ Some tests failed!"
    exit 1
fi
