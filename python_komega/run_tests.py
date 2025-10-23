#!/usr/bin/env python3
"""
Test Runner for Komega Python Library

This script runs all available tests for the Komega Python library,
including basic functionality tests and Fortran compatibility tests.

Usage:
    python run_tests.py [--test-type TYPE] [--verbose]

Test Types:
    basic       - Basic functionality tests
    fortran     - Fortran compatibility tests
    detailed    - Detailed Fortran compatibility tests
    all         - Run all tests (default)

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import argparse
import os
import sys
from typing import List, Optional

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def run_basic_tests() -> bool:
    """Run basic functionality tests."""
    print("=" * 60)
    print("Running Basic Functionality Tests")
    print("=" * 60)

    try:
        from test_komega import main as test_main

        result = test_main()
        return result == 0
    except Exception as e:
        print(f"Basic tests failed: {e}")
        return False


def run_fortran_tests() -> bool:
    """Run Fortran compatibility tests."""
    print("=" * 60)
    print("Running Fortran Compatibility Tests")
    print("=" * 60)

    try:
        from test_fortran_compatibility import main as test_main

        result = test_main()
        return result == 0
    except Exception as e:
        print(f"Fortran compatibility tests failed: {e}")
        return False


def run_detailed_tests() -> bool:
    """Run detailed Fortran compatibility tests."""
    print("=" * 60)
    print("Running Detailed Fortran Compatibility Tests")
    print("=" * 60)

    try:
        from test_detailed_compatibility import main as test_main

        result = test_main()
        return result == 0
    except Exception as e:
        print(f"Detailed compatibility tests failed: {e}")
        return False


def run_all_tests() -> bool:
    """Run all available tests."""
    print("=" * 80)
    print("Komega Python Library - Complete Test Suite")
    print("=" * 80)
    print()

    test_results = []

    # Run basic tests
    print("1. Basic Functionality Tests")
    basic_success = run_basic_tests()
    test_results.append(("Basic Functionality", basic_success))
    print()

    # Run Fortran compatibility tests
    print("2. Fortran Compatibility Tests")
    fortran_success = run_fortran_tests()
    test_results.append(("Fortran Compatibility", fortran_success))
    print()

    # Run detailed tests
    print("3. Detailed Fortran Compatibility Tests")
    detailed_success = run_detailed_tests()
    test_results.append(("Detailed Compatibility", detailed_success))
    print()

    # Summary
    print("=" * 80)
    print("Test Summary")
    print("=" * 80)

    all_passed = True
    for test_name, success in test_results:
        status = "PASSED" if success else "FAILED"
        print(f"{test_name:25s}: {status}")
        if not success:
            all_passed = False

    print()
    if all_passed:
        print("🎉 All tests passed successfully!")
        print("The Python implementation is working correctly.")
    else:
        print("❌ Some tests failed.")
        print("Please check the error messages above.")

    print("=" * 80)

    return all_passed


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(
        description="Test Runner for Komega Python Library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_tests.py                    # Run all tests
  python run_tests.py --test-type basic  # Run only basic tests
  python run_tests.py --test-type fortran # Run only Fortran compatibility tests
  python run_tests.py --test-type detailed # Run only detailed tests
  python run_tests.py --verbose          # Run with verbose output
        """,
    )

    parser.add_argument(
        "--test-type",
        choices=["basic", "fortran", "detailed", "all"],
        default="all",
        help="Type of tests to run (default: all)",
    )

    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    # Set up verbose output if requested
    if args.verbose:
        import logging

        logging.basicConfig(level=logging.DEBUG)

    # Run selected tests
    success = False

    if args.test_type == "basic":
        success = run_basic_tests()
    elif args.test_type == "fortran":
        success = run_fortran_tests()
    elif args.test_type == "detailed":
        success = run_detailed_tests()
    elif args.test_type == "all":
        success = run_all_tests()
    else:
        print(f"Unknown test type: {args.test_type}")
        return 1

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
