# HeaderPrinter Utility - Usage Guide

## Overview

The `HeaderPrinter` class consolidates the 30+ duplicate `print_header()` implementations across the Psi4 codebase into a single, flexible utility. This reduces code duplication and provides a consistent interface for printing formatted headers.

## Location

- **Header**: `psi4/libpsi4util/header_printer.h`
- **Implementation**: `psi4/libpsi4util/header_printer.cc`
- **Include**: Already included via `psi4/libpsi4util/libpsi4util.h`

## Basic Usage

### 1. Simplest Case - No Parameters

For headers with no parameters, use the temporary object pattern:

```cpp
#include "psi4/libpsi4util/header_printer.h"

void MyClass::print_header() const {
    HeaderPrinter("My Algorithm").print();
}

// Or with conditional printing:
void MyClass::print_header() const {
    HeaderPrinter("My Algorithm").print_if(print_);
}
```

**Output:**
```
  ==> My Algorithm <==

```

### 2. Simple Arrow-Style Header with Parameters

```cpp
void MyClass::print_header() const {
    HeaderPrinter header("My Algorithm");
    header.add_parameter("Cutoff", cutoff_value_)
          .print_if(print_);
}
```

**Output:**
```
  ==> My Algorithm <==

    Cutoff:                  1.0E-12
```

### 3. Header with Multiple Parameters

```cpp
void MyJK::print_header() const {
    HeaderPrinter header("My JK Builder");
    header.add_parameter("J tasked", do_J_)              // bool -> "Yes"/"No"
          .add_parameter("K tasked", do_K_)              // bool -> "Yes"/"No"
          .add_parameter("Threads", nthreads_)           // int
          .add_parameter("Memory [MiB]", memory_mb_)     // long
          .add_parameter("Cutoff", cutoff_, "%11.0E")    // double with format
          .print_if(print_);
}
```

**Output:**
```
  ==> My JK Builder <==

    J tasked:                    Yes
    K tasked:                     No
    Threads:                       4
    Memory [MiB]:               1024
    Cutoff:                 1.0E-12
```

### 4. Box-Style Header with Authors (DF-MP2 Style)

```cpp
void RDFMP2::print_header() {
    int nthread = Process::environment.get_n_threads();
    char thread_info[64];
    snprintf(thread_info, sizeof(thread_info), "RMP2 Wavefunction, %3d Threads", nthread);

    HeaderPrinter header("DF-MP2", HeaderPrinter::BannerStyle::BOX, 57);
    header.subtitle("2nd-Order Density-Fitted Moller-Plesset Theory")
          .add_line("", true)  // blank centered line
          .add_line(thread_info, true)
          .add_line("", true)
          .add_authors({
              "Rob Parrish, Justin Turney, Andy Simmonett,",
              "Ed Hohenstein, and C. David Sherrill"
          })
          .print();
}
```

**Output:**
```
	 ---------------------------------------------------------
	                          DF-MP2
	      2nd-Order Density-Fitted Moller-Plesset Theory

	              RMP2 Wavefunction,   4 Threads

	        Rob Parrish, Justin Turney, Andy Simmonett,
	           Ed Hohenstein, and C. David Sherrill
	 ---------------------------------------------------------
```

## Banner Styles

### BannerStyle::ARROW (Default)
Standard `==> Title <==` format used by most algorithms.

### BannerStyle::BOX
Full box with dashes, centered text, and author support (like DF-MP2).

### BannerStyle::NONE
No banner, just print content.

### BannerStyle::CUSTOM
User-provided custom banner strings.

```cpp
HeaderPrinter header("Title");
header.set_custom_banner("*****************************",
                         "*****************************")
      .add_parameter("Value", 42)
      .print();
```

## Best Practices

### When to Use Temporary Object Pattern

Use the simplified form `HeaderPrinter("Title").print()` when:
- No parameters need to be added
- No conditional logic is required
- The header is printed immediately

**Good:**
```cpp
HeaderPrinter("SAP guess").print();
HeaderPrinter("DFT Potential").print_if(print_ > 0);
```

### When to Use Named Variable

Use a named variable when:
- Adding parameters (especially with conditionals)
- Building complex headers incrementally
- The code is clearer with intermediate steps

**Good:**
```cpp
HeaderPrinter header("CDJK: Cholesky-decomposed J/K Matrices");
header.add_parameter("J tasked", do_J_)
      .add_parameter("K tasked", do_K_);
if (do_wK_) {
    header.add_parameter("Omega", omega_);
}
header.print();
```

## Method Chaining

All `add_*` methods return a reference to the HeaderPrinter object, enabling method chaining:

```cpp
header.add_parameter("Param1", value1)
      .add_parameter("Param2", value2)
      .add_separator()
      .add_parameter("Param3", value3)
      .print_if(condition);
```

## Available Methods

### Construction
- `HeaderPrinter(title, style=ARROW, width=60)` - Create a header printer

### Configuration
- `subtitle(text)` - Set subtitle (BOX style only)
- `set_custom_banner(top, bottom="")` - Set custom banners

### Content
- `add_line(text, centered=false)` - Add a custom text line
- `add_parameter(key, value, key_width=20)` - Add parameter (multiple overloads)
  - `add_parameter(key, string_value, key_width)`
  - `add_parameter(key, int_value, key_width)`
  - `add_parameter(key, long_value, key_width)`
  - `add_parameter(key, double_value, format="%11.3E", key_width)`
  - `add_parameter(key, bool_value, key_width)` - Prints "Yes"/"No"
- `add_authors(vector<string>)` - Add author credits (BOX style)
- `add_separator()` - Add separator line
- `add_blank_line()` - Add blank line

### Output
- `print()` - Print to `outfile`
- `print_if(condition)` - Conditionally print
- `print_to(stream)` - Print to specific output stream

## Migration Examples

### Before: LinK
```cpp
void LinK::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> LinK: Linear Exchange K <==\n\n");
        outfile->Printf("    K Screening Cutoff:%11.0E\n", linK_ints_cutoff_);
    }
}
```

### After: LinK
```cpp
void LinK::print_header() const {
    HeaderPrinter header("LinK: Linear Exchange K");
    header.add_parameter("K Screening Cutoff", linK_ints_cutoff_)
          .print_if(print_);
}
```

### Before: MemDFJK
```cpp
void MemDFJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> MemDFJK: Density-Fitted J/K Matrices <==\n\n");
        outfile->Printf("    J tasked:           %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:           %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:          %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:              %11.3E\n", omega_);
        outfile->Printf("    OpenMP threads:     %11d\n", omp_nthread_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:          %11s\n", (dfh_->get_AO_core() ? "Core" : "Disk"));
        outfile->Printf("    Schwarz Cutoff:     %11.0E\n", cutoff_);
        outfile->Printf("    Mask sparsity (%%):  %11.4f\n", 100. * dfh_->ao_sparsity());
        outfile->Printf("    Fitting Condition:  %11.0E\n\n", condition_);
    }
}
```

### After: MemDFJK
```cpp
void MemDFJK::print_header() const {
    if (print_) {
        HeaderPrinter header("MemDFJK: Density-Fitted J/K Matrices");
        header.add_parameter("J tasked", do_J_)
              .add_parameter("K tasked", do_K_)
              .add_parameter("wK tasked", do_wK_);
        if (do_wK_) {
            header.add_parameter("Omega", omega_, "%11.3E");
        }
        header.add_parameter("OpenMP threads", omp_nthread_)
              .add_parameter("Memory [MiB]", (memory_ * 8L) / (1024L * 1024L))
              .add_parameter("Algorithm", dfh_->get_AO_core() ? "Core" : "Disk")
              .add_parameter("Schwarz Cutoff", cutoff_, "%11.0E")
              .add_parameter("Mask sparsity (%)", 100. * dfh_->ao_sparsity(), "%11.4f")
              .add_parameter("Fitting Condition", condition_, "%11.0E")
              .print();
    }
}
```

## Benefits

1. **Reduced Code Duplication**: Consolidates 30+ similar implementations
2. **Consistent Formatting**: Ensures uniform output across all modules
3. **Easier Maintenance**: Changes to formatting only need to be made once
4. **Type Safety**: Overloaded methods handle different types correctly
5. **Cleaner Code**: Method chaining creates more readable code
6. **Flexibility**: Supports multiple banner styles and custom formatting

## Remaining Implementations to Migrate

The following modules still have custom `print_header()` implementations that could benefit from using HeaderPrinter:

### libfock (18 implementations)
- DirectJK, DiskJK, DiskDFJK, CDJK, PKJK, COSK, CompositeJK
- DirectDFJ, snLinK
- VBase, SAP, RV, UV
- CPHFRHamiltonian, CGRSolver, RCPHF
- DFJKGrad, DirectJKGrad

### dfmp2 (2 more implementations)
- UDFMP2, RODFMP2, DFCorrGrad

### lib3index (2 implementations)
- DFHelper, DFTensor

### Other modules
- HF (libscf_solver)
- SAPT0, USAPT0, SAPT2, SAPT2p, SAPT2p3 (libsapt_solver)
- DLPNOMP2, MP2F12
- FISAPT, IBOLocalizer2
- CubeProperties, CubicScalarGrid
- OEProp, BoysLocalizer, PMLocalizer
- LS_THC_Computer

## Notes

- The HeaderPrinter uses `outfile` by default (via `psi4::outfile`)
- Custom output streams can be specified with `print_to(stream)`
- All formatting is consistent with existing Psi4 conventions
- The utility is fully backward compatible with existing output

## Further Improvements

Future enhancements could include:
- Table formatting utilities for data tables
- Automatic width adjustment based on content
- Template specializations for common parameter types
- Integration with Psi4's logging levels
