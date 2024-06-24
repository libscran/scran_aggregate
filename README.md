# Aggregate expression values across cells

![Unit tests](https://github.com/libscran/aggregate_across_cells/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/aggregate_across_cells/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/aggregate_across_cells/graph/badge.svg?token=JWV0I4WJX2)](https://codecov.io/gh/libscran/aggregate_across_cells)

## Overview

This repository contains a function to aggregate statistics for groups of cells from a gene-by-cell matrix of expression values.
It was primarily developed for computing pseudo-bulk expression profiles for clusters of cells,
which can then be used for differential expression analysis.
The code itself was originally derived from the [**scran** R package](https://bioconductor.org/packages/scran),
factored out into a separate C++ library for easier re-use.

## Quick start

Not too much to say here.
Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami) and an array of group assignments,
the `aggregate_across_cells::compute()` function will compute the aggregate statistics across all genes for each group.

```cpp
#include "scran/aggregate_across_cells.hpp"

tatami::Matrix<double, int>* ptr = some_data_source();
std::vector<int> groupings = some_groupings();

auto res = scran::aggregate_across_cells::compute(ptr, groupings.data());
res.sums; // vector of vectors of per-group sums across genes.
res.sums[0]; // vector of sums for the first group across all genes.
res.detected; // vector of vectors of the number of detected cells. 
```

The array of groupings should contain integer assignments to groups 0, 1, 2, etc.
For more complex groupings defined from combinations of multiple factors, 
the `combine_factors::compute()` utility will create group assignments from unique combinations of those factors:

```cpp
#include "scran/combine_factors.hpp"

std::vector<int> grouping1 { 0, 0, 1, 1, 2, 2 };
std::vector<int> grouping2 { 0, 1, 0, 1, 0, 1 };

std::vector<int> combined(grouping1.size()); 
auto res = scran::combine_factors::compute(
    grouping1.size(), 
    std::vector<int*>{ grouping1.data(), grouping2.data() },
    combined.data()
);

combined; // defines unique combinations of (grouping1, grouping2).
res.factors[0]; // values of grouping1 for each unique combination.
res.factors[1]; // values of grouping2 for each unique combination.
```

Check out the [reference documentation](https://libscran.github.io/aggregate_across_cells) for more details.

## Building projects

This repository is part of the broader [**libscran**](https://github.com/libscran/libscran) library,
so users are recommended to use the latter in their projects.
Nonetheless, if this library is to be used directly, we can use CMake with `FetchContent`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_aggregate_across_cells 
  GIT_REPOSITORY https://github.com/libscran/aggregate_across_cells
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_aggregate_across_cells)
```

Then we can link to **aggregate_across_cells** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe scran::aggregate_across_cells)

# For libaries
target_link_libraries(mylib INTERFACE scran::aggregate_across_cells)
```

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
