#!/usr/bin/env python3
"""
Test result reporter for RAMSES-CPP test suite.
Adapted from legacy RAMSES utilities for C++ version.

Usage:
    python3 report_test_results.py <test_log_file>
"""

import sys
import re
from collections import defaultdict
from pathlib import Path

def parse_test_results(log_file):
    """Parse test suite log file and extract PASSED/FAILED results."""
    results = {
        'passed': [],
        'failed': [],
        'total': 0
    }

    with open(log_file, 'r') as f:
        for line in f:
            # Match lines like: "Test hydro/sod-tube passed            [ OK ]"
            if re.search(r'passed\s+\[\s*OK\s*\]', line):
                match = re.search(r'Test (.+?)\s+passed', line)
                if match:
                    results['passed'].append(match.group(1))
                    results['total'] += 1
            # Match lines like: "Test hydro/mixing-scalar failed!      [FAIL]"
            elif re.search(r'failed!\s+\[FAIL\]', line):
                match = re.search(r'Test (.+?)\s+failed!', line)
                if match:
                    results['failed'].append(match.group(1))
                    results['total'] += 1

    return results

def categorize_tests(results):
    """Categorize tests by type (hydro, mhd, etc.)."""
    categories = defaultdict(lambda: {'passed': [], 'failed': []})

    for test in results['passed']:
        category = test.split('/')[0] if '/' in test else test
        categories[category]['passed'].append(test)

    for test in results['failed']:
        category = test.split('/')[0] if '/' in test else test
        categories[category]['failed'].append(test)

    return categories

def print_report(results, categories):
    """Print formatted test results report."""
    total_passed = len(results['passed'])
    total_failed = len(results['failed'])
    total_tests = results['total']

    print("\n" + "="*70)
    print("RAMSES-CPP TEST SUITE RESULTS")
    print("="*70)

    print(f"\n📊 OVERALL SUMMARY:")
    print(f"  Total Tests:  {total_tests}")
    print(f"  ✅ Passed:    {total_passed}")
    print(f"  ❌ Failed:    {total_failed}")
    print(f"  Success Rate: {100*total_passed//total_tests if total_tests > 0 else 0}%")

    print(f"\n📋 RESULTS BY CATEGORY:")
    for category in sorted(categories.keys()):
        cat_passed = len(categories[category]['passed'])
        cat_failed = len(categories[category]['failed'])
        cat_total = cat_passed + cat_failed
        print(f"\n  {category.upper()}:")
        print(f"    {cat_passed}/{cat_total} passed")

        if categories[category]['failed']:
            print(f"    Failed tests:")
            for test in categories[category]['failed']:
                print(f"      - {test}")

    print("\n" + "="*70)

    if total_failed == 0:
        print("🎉 ALL TESTS PASSED! 🎉")
    else:
        print(f"⚠️  {total_failed} test(s) need attention")

    print("="*70 + "\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 report_test_results.py <test_log_file>")
        sys.exit(1)

    log_file = sys.argv[1]

    if not Path(log_file).exists():
        print(f"Error: Log file not found: {log_file}")
        sys.exit(1)

    results = parse_test_results(log_file)
    categories = categorize_tests(results)
    print_report(results, categories)

if __name__ == "__main__":
    main()
