# Aggregate expression values across cells

![Unit tests](https://github.com/libscran/scran_aggregate/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/libscran/scran_aggregate/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/libscran/scran_aggregate/graph/badge.svg?token=JWV0I4WJX2)](https://codecov.io/gh/libscran/scran_aggregate)

## Overview

This repository contains a function to aggregate statistics for groups of cells from a gene-by-cell matrix of expression values.
It was primarily developed for computing pseudo-bulk expression profiles for clusters of cells,
which can then be used for differential expression analysis.
The code itself was originally derived from the [**scran** R package](https://bioconductor.org/packages/scran),
factored out into a separate C++ library for easier re-use.

## Quick start

Given a [`tatami::Matrix`](https://github.com/tatami-inc/tatami) and an array of group assignments,
the `aggregate_across_cells()` function will compute the aggregate statistics across all genes for each group.

```cpp
#include "scran_aggregate/scran_aggregate.hpp"

const tatami::Matrix<double, int>& mat = some_data_source();
std::vector<int> groupings = some_groupings();

scran_aggregate::AggregateAcrossCellsOptions opt;
auto res = scran_aggregate::aggregate_across_cells(mat, groupings.data(), opt);

res.sums; // vector of vectors of per-group sums across genes.
res.sums[0]; // vector of sums for the first group across genes.
res.detected; // vector of vectors of the number of detected cells per gene.
```

The array of groupings should contain integer assignments to groups 0, 1, 2, etc.
For more complex groupings defined from combinations of multiple factors, 
the `combine_factors()` utility will create group assignments from unique combinations of those factors:

```cpp
std::vector<int> grouping1 { 0, 0, 1, 1, 2, 2 };
std::vector<int> grouping2 { 0, 1, 0, 1, 0, 1 };

std::vector<int> combined(grouping1.size()); 
auto res = scran_aggregate::combine_factors(
    grouping1.size(), 
    std::vector<int*>{ grouping1.data(), grouping2.data() },
    combined.data()
);

combined; // defines unique combinations of (grouping1, grouping2).
res.factors[0]; // values of grouping1 for each unique combination.
res.factors[1]; // values of grouping2 for each unique combination.
```

We can also use the `aggregate_across_genes()` function to sum expression values across gene sets, e.g., to compute the activity of a gene signature.
This can be done with any number of gene sets, possibly with a different weight for each gene in each set.

```cpp
std::vector<std::tuple<size_t, const int*, const double*> > gene_sets;

std::vector<int> set1 { 0, 5, 10, 20 };
gene_sets.emplace_back(set1.size(), set1.data(), static_cast<double*>(NULL)); // no weight

std::vector<int> set2 { 0, 2, 4, 6, 8, 10 };
std::vector<double> weight2 { 0.1, 0.3, 0.3, 0.2, 0.1, 0.05 };
gene_sets.emplace_back(set2.size(), set2.data(), weight2.data()); // weighted

scran_aggregate::AggregateAcrossGenesOptions g_opt;
auto g_res = scran_aggregate::aggregate_across_genes(
    mat,
    gene_sets,
    g_opt
);
g_res.sum[0]; // vector of sums for set 1 in each cell.
```

Check out the [reference documentation](https://libscran.github.io/scran_aggregate) for more details.

## Building projects

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  scran_aggregate
  GIT_REPOSITORY https://github.com/libscran/scran_aggregate
  GIT_TAG master # or any version of interest
)

FetchContent_MakeAvailable(scran_aggregate)
```

Then you can link to **scran_aggregate** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe libscran::scran_aggregate)

# For libaries
target_link_libraries(mylib INTERFACE libscran::scran_aggregate)
```

### CMake with `find_package()`

```cmake
find_package(libscran_scran_aggregate CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE libscran::scran_aggregate)
```

To install the library, use:

```sh
mkdir build && cd build
cmake .. -DSCRAN_AGGREGATE_TESTS=OFF
cmake --build . --target install
```

By default, this will use `FetchContent` to fetch all external dependencies.
If you want to install them manually, use `-DSCRAN_AGGREGATE_FETCH_EXTERN=OFF`.
See the tags in [`extern/CMakeLists.txt`](extern/CMakeLists.txt) to find compatible versions of each dependency.

### Manual

If you're not using CMake, the simple approach is to just copy the files in `include/` - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt), which also need to be made available during compilation.
