# Fix for File Handle Leak in psi4.set_output_file()

## Problem

Using `psi4.set_output_file()` in a loop causes file handle leaks that eventually lead to a segmentation fault when the number of open files exceeds the OS limit (ulimit).

### Root Cause

The issue occurs because:

1. **Python layer** (`psi4.set_output_file()` in `psi4/extras.py`):
   - Creates a new logging `FileHandler` for each `.log` file
   - Adds this handler to `psi4.logger`
   - **Never removes or closes old handlers**

2. **C++ layer** (`psi4.core.close_outfile()` in `psi4/src/core.cc`):
   - Only closes the C++ `PsiOutStream` (for `.dat`/`.out` ASCII files)
   - **Does NOT affect Python logging handlers**

3. **Result**:
   - Each call to `set_output_file()` adds a new `FileHandler`
   - Old `FileHandlers` remain open indefinitely
   - File handles accumulate until hitting OS limit
   - System runs out of file descriptors → segfault

## Solution

Modified `psi4.set_output_file()` in `psi4/extras.py` to:

1. **Remove all old FileHandlers** before adding a new one
2. **Properly close** each FileHandler before removal
3. **Only affect FileHandlers**, not other handler types (StreamHandler, etc.)

### Code Changes

**File**: `psi4/extras.py`
**Function**: `set_output_file()` (lines 250-322)

**Before** (lines 309-317):
```python
if execute:
    core.set_output_file(str(out), append)
    if print_header is True or (print_header is None and not append):
        _print_header()
        core.print_out("Addons:     " + textwrap.fill(", ".join(addons()), width=95, initial_indent='', subsequent_indent='                ') + "\n")
    # Warning: baseFilename is not part of the documented API for the logging module and could change.
    filenames = [handle.baseFilename for handle in logger.handlers]
    if not f_handler.baseFilename in filenames:
        logger.addHandler(f_handler)
```

**After** (lines 309-321):
```python
if execute:
    core.set_output_file(str(out), append)
    if print_header is True or (print_header is None and not append):
        _print_header()
        core.print_out("Addons:     " + textwrap.fill(", ".join(addons()), width=95, initial_indent='', subsequent_indent='                ') + "\n")
    # Remove old FileHandlers to prevent file handle leaks
    # Warning: baseFilename is not part of the documented API for the logging module and could change.
    for handler in logger.handlers[:]:  # Iterate over a copy
        if isinstance(handler, logging.FileHandler):
            handler.close()
            logger.removeHandler(handler)
    # Add the new handler
    logger.addHandler(f_handler)
```

### Key Changes

1. **Removed conditional check**: No longer checks if handler already exists by filename
2. **Added cleanup loop**: Iterates through all existing handlers
3. **Closes FileHandlers**: Calls `handler.close()` to release file handles
4. **Removes handlers**: Calls `logger.removeHandler()` to clean up
5. **Iterates over copy**: Uses `logger.handlers[:]` to avoid modifying list during iteration
6. **Type-specific**: Only affects `FileHandler` instances, not other handler types
7. **Unconditional add**: Always adds the new handler (no duplicates possible after cleanup)

## Testing

Created comprehensive test suite in `tests/pytests/test_output_file_handler_leak.py`:

1. **test_output_file_handler_cleanup()**: Verifies handlers don't accumulate in a loop
2. **test_output_file_handler_cleanup_with_computation()**: Tests with actual SCF calculations
3. **test_output_file_switching()**: Tests switching between different files
4. **test_output_file_append_mode()**: Tests append mode doesn't leak handlers

All tests verify that:
- Number of FileHandlers never exceeds `initial + 1`
- Total handler count doesn't grow significantly
- File handles are properly released

## Benefits

1. **Prevents segfaults**: No more file descriptor exhaustion
2. **Loop-safe**: Can call `set_output_file()` indefinitely in loops
3. **Minimal impact**: Only affects FileHandlers for .log files
4. **Backward compatible**: No API changes, existing code works unchanged
5. **Simple implementation**: Clean, maintainable code

## Related Files

- **Modified**: `psi4/extras.py` - Main fix
- **Added**: `tests/pytests/test_output_file_handler_leak.py` - Test suite
- **Unchanged**: `psi4/src/core.cc` - C++ layer works as-is
- **Unchanged**: `psi4/src/psi4/libpsi4util/PsiOutStream.{h,cc}` - C++ output stream

## Usage Example

```python
import psi4

mol = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
""")

psi4.set_options({'basis': 'sto-3g'})

# This now works correctly without leaking file handles
for i in range(1000):
    psi4.set_output_file(f"output_{i}.dat")
    psi4.energy('scf')
    psi4.core.close_outfile()
```

## Notes

- The comment about `baseFilename` not being part of the documented API is preserved
- The fix maintains the existing behavior for all other aspects of `set_output_file()`
- No changes needed to C++ code or `close_outfile()` function
- The Python logger maintains at most ONE FileHandler at any time
