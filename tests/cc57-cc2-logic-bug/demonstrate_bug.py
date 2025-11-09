#!/usr/bin/env python3
"""
Demonstration script showing the boolean logic bug in t2.cc:48

This script demonstrates why the buggy condition is always true.
"""

def buggy_logic(wfn):
    """
    Simulates the BUGGY logic from t2.cc:48
    Returns True if CCSD terms should be executed
    """
    # BUGGY: if (params_.wfn != "CC2" || params_.wfn != "EOM_CC2")
    return wfn != "CC2" or wfn != "EOM_CC2"

def correct_logic(wfn):
    """
    Simulates the CORRECT logic (with AND instead of OR)
    Returns True if CCSD terms should be executed
    """
    # CORRECT: if (params_.wfn != "CC2" && params_.wfn != "EOM_CC2")
    return wfn != "CC2" and wfn != "EOM_CC2"

def expected_behavior(wfn):
    """
    What we WANT: Skip CCSD terms for CC2/EOM_CC2, execute for others
    """
    return wfn not in ["CC2", "EOM_CC2"]

def main():
    print("=" * 70)
    print("Boolean Logic Bug Demonstration: t2.cc:48")
    print("=" * 70)
    print()
    print("The bug: if (params_.wfn != \"CC2\" || params_.wfn != \"EOM_CC2\")")
    print("The fix: if (params_.wfn != \"CC2\" && params_.wfn != \"EOM_CC2\")")
    print()
    print("Testing different wavefunction types:")
    print("-" * 70)
    print(f"{'Method':<12} | {'Expected':<8} | {'Buggy OR':<9} | {'Fixed AND':<10} | Status")
    print("-" * 70)

    test_methods = ["CC2", "EOM_CC2", "CCSD", "CC3", "BCCD", "CCSD_T"]

    for wfn in test_methods:
        expected = expected_behavior(wfn)
        buggy = buggy_logic(wfn)
        fixed = correct_logic(wfn)

        # Check if buggy logic gives wrong result
        if buggy != expected:
            status = "❌ BUG"
        else:
            status = "✓ OK"

        # Also check if fixed logic matches expected
        if fixed != expected:
            status += " (FIX WRONG!)"
        else:
            if status == "❌ BUG":
                status += " → ✓ FIXED"

        print(f"{wfn:<12} | {str(expected):<8} | {str(buggy):<9} | {str(fixed):<10} | {status}")

    print("-" * 70)
    print()
    print("Explanation:")
    print()
    print("For CC2:")
    print("  wfn != 'CC2' → False")
    print("  wfn != 'EOM_CC2' → True")
    print("  False OR True → True  ❌ (should be False!)")
    print("  False AND True → False ✓ (correct!)")
    print()
    print("For EOM_CC2:")
    print("  wfn != 'CC2' → True")
    print("  wfn != 'EOM_CC2' → False")
    print("  True OR False → True  ❌ (should be False!)")
    print("  True AND False → False ✓ (correct!)")
    print()
    print("For CCSD:")
    print("  wfn != 'CC2' → True")
    print("  wfn != 'EOM_CC2' → True")
    print("  True OR True → True ✓ (correct)")
    print("  True AND True → True ✓ (also correct)")
    print()
    print("=" * 70)
    print("Conclusion: The buggy OR logic is ALWAYS TRUE!")
    print("=" * 70)
    print()
    print("Impact:")
    print("  - If t2_build() is called for CC2/EOM_CC2, it will incorrectly")
    print("    execute CCSD-specific terms (FaetT2, FmitT2, WmnijT2, etc.)")
    print("  - This breaks the CC2 approximation (T1-only intermediates)")
    print("  - Results would be incorrect (neither proper CC2 nor CCSD)")
    print()
    print("Current status:")
    print("  - Normal CC2 calculations use cc2_t2_build() → Not affected ✓")
    print("  - Edge cases (e.g., one_step with just_residuals) → Could be affected ❌")
    print()

if __name__ == "__main__":
    main()
