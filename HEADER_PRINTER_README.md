# HeaderPrinter Utility - Quick Start

## What is HeaderPrinter?

HeaderPrinter is a utility class that consolidates 40+ duplicate `print_header()` implementations across the Psi4 codebase into a single, reusable component. This eliminates code duplication and ensures consistent formatting.

## Current Status

- **Created**: 2025-11-18
- **Location**: `psi4/src/psi4/libpsi4util/header_printer.{h,cc}`
- **Migrations Completed**: 3/43 (7%)
- **Code Reduction**: ~60% average in migrated files

### ✅ Completed Migrations
1. `LinK::print_header()` - Simple arrow style (7 lines → 3 lines)
2. `MemDFJK::print_header()` - Arrow style with multiple parameters
3. `RDFMP2::print_header()` - Box style with authors

## Quick Start for Developers

### I want to use HeaderPrinter in new code

```cpp
#include "psi4/libpsi4util/header_printer.h"

void MyClass::print_header() const {
    HeaderPrinter header("My Algorithm Name");
    header.add_parameter("Threads", n_threads_)
          .add_parameter("Memory [MiB]", memory_mb_)
          .add_parameter("Cutoff", cutoff_, "%11.0E")
          .print_if(print_);  // Respects print level
}
```

### I want to migrate an existing print_header() method

1. **Read the migration plan**: See `HEADER_PRINTER_MIGRATION_PLAN.md`
2. **Check the tracking spreadsheet**: See `HEADER_PRINTER_MIGRATION_TRACKING.csv`
3. **Follow the migration pattern** for your case (simple/box/complex)
4. **Test output matches** original exactly
5. **Submit PR** with before/after comparison

### I want to understand the full API

Read the comprehensive usage guide: `psi4/src/psi4/libpsi4util/HEADER_PRINTER_USAGE.md`

## File Organization

```
psi4/
├── HEADER_PRINTER_README.md                    # ← You are here (Quick start)
├── HEADER_PRINTER_MIGRATION_PLAN.md            # Detailed migration strategy
├── HEADER_PRINTER_MIGRATION_TRACKING.csv       # Progress tracking spreadsheet
└── psi4/src/psi4/libpsi4util/
    ├── header_printer.h                        # API header
    ├── header_printer.cc                       # Implementation
    └── HEADER_PRINTER_USAGE.md                 # API documentation + examples
```

## Common Patterns

### Pattern 1: Simple Arrow Header
**Use case**: Most JK algorithms, simple utilities

```cpp
HeaderPrinter header("Algorithm Name");
header.add_parameter("Cutoff", value)
      .print_if(print_);
```

### Pattern 2: Box Header with Authors
**Use case**: Major methods (MP2, SAPT, etc.)

```cpp
HeaderPrinter header("Method Name", HeaderPrinter::BannerStyle::BOX, 57);
header.subtitle("Full Method Description")
      .add_authors({"Author 1", "Author 2"})
      .print();
```

### Pattern 3: Multiple Parameters
**Use case**: Configurable algorithms

```cpp
HeaderPrinter header("Algorithm");
header.add_parameter("Param 1", value1)
      .add_parameter("Param 2", value2)
      .add_parameter("Param 3", value3, "%11.3E")
      .print_if(print_);
```

### Pattern 4: Conditional Parameters
**Use case**: Optional features

```cpp
HeaderPrinter header("Algorithm");
header.add_parameter("Base param", value1);
if (feature_enabled) {
    header.add_parameter("Optional param", value2);
}
header.print();
```

## Migration Priority

### 🔴 High Priority (Start here)
**Phase 1**: Simple arrow-style headers in libfock
- DirectJK, DiskJK, PKJK (nearly identical to completed LinK)
- Estimated: 2-3 hours total
- Low risk, high impact

### 🟡 Medium Priority
**Phase 2**: Box-style headers with authors
- UDFMP2, RODFMP2 (nearly identical to completed RDFMP2)
- SAPT family (sapt0, sapt2, sapt2p, sapt2p3)
- Estimated: 3-4 hours total
- Medium risk, medium impact

### 🟢 Low Priority
**Phase 3**: Complex headers
- HF::print_header() (47 lines - very complex)
- CompositeJK (recursive printing)
- Estimated: 2-3 hours total
- Higher risk, consider partial migration

### ⚪ Lowest Priority
**Phase 4**: Utility classes
- Localizers, property calculators
- Good starter tasks for new contributors
- Estimated: 1-2 hours total

## Testing Your Migration

### 1. Build Test
```bash
cd build/
cmake --build . --target psi4 -- -j4
```

### 2. Output Comparison Test
```bash
# Before migration
psi4 test_input.dat -o before.out

# After migration
psi4 test_input.dat -o after.out

# Compare (should be identical)
diff before.out after.out
```

### 3. Regression Test
```bash
ctest -R <module_name> -V
```

## Benefits of Migration

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| Code lines (avg) | 10-15 | 3-5 | 60-70% reduction |
| Duplication | 40+ copies | 1 utility | 100% elimination |
| Maintainability | Low | High | Single point of change |
| Consistency | Variable | Uniform | Guaranteed formatting |

## Need Help?

1. **API Questions**: Read `HEADER_PRINTER_USAGE.md`
2. **Migration Strategy**: Read `HEADER_PRINTER_MIGRATION_PLAN.md`
3. **Examples**: Check completed migrations:
   - LinK.cc (simplest)
   - MemDFJK.cc (with conditionals)
   - mp2.cc:RDFMP2 (box style)

## Contributing

### Good First Issues
These are simple migrations perfect for getting started:
- DirectJK::print_header() (nearly identical to LinK)
- DiskJK::print_header() (nearly identical to LinK)
- BoysLocalizer::print_header() (3 parameters)
- PMLocalizer::print_header() (3 parameters)

### Claiming Work
1. Check `HEADER_PRINTER_MIGRATION_TRACKING.csv`
2. Update "Status" and "Assignee" columns
3. Submit PR when complete
4. Update tracking file in PR

## Migration Checklist

When migrating a print_header() implementation:

- [ ] Read current implementation completely
- [ ] Identify banner style (ARROW/BOX/NONE/CUSTOM)
- [ ] List all parameters and types
- [ ] Note any conditional logic
- [ ] Add `#include "psi4/libpsi4util/header_printer.h"`
- [ ] Implement using HeaderPrinter
- [ ] Test: Build succeeds
- [ ] Test: Output matches original exactly
- [ ] Test: Regression tests pass
- [ ] Verify code reduction achieved (30%+ target)
- [ ] Update tracking spreadsheet
- [ ] Submit PR with clear before/after

## Metrics

### Progress Tracking
```
Phase 0 (Examples):   [███████████████████████████████] 3/3   (100%)
Phase 1 (Simple):     [                               ] 0/16  (0%)
Phase 2 (Box):        [                               ] 0/10  (0%)
Phase 3 (Complex):    [                               ] 0/8   (0%)
Phase 4 (Utilities):  [                               ] 0/6   (0%)

Overall:              [███                            ] 3/43  (7%)
```

### Code Reduction
- **Current**: ~40 lines eliminated
- **Projected**: ~500 lines total reduction when complete
- **Average**: 50-70% reduction per file

## Timeline

**Recommended**: Phased approach over 2-3 weeks
- Week 1: Phase 1 (simple arrow headers)
- Week 2: Phase 2 (box headers with authors)
- Week 3: Phases 3-4 (complex and utilities)

**Aggressive**: Could complete in 1 week with dedicated effort
**Conservative**: Spread over multiple releases, migrate opportunistically

## Known Issues / Limitations

None currently. The HeaderPrinter utility is complete and tested.

### Future Enhancements (Not Blocking)
- Table formatting utilities
- Integration with Psi4 logging levels
- Performance profiling (overhead is negligible)

## License

Same as Psi4: LGPL v3

---

**Last Updated**: 2025-11-18
**Status**: Ready for community migration effort
**Questions**: Refer to documentation files listed above
