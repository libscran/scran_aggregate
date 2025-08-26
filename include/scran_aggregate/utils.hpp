#ifndef SCRAN_AGGREGATE_UTILS_HPP
#define SCRAN_AGGREGATE_UTILS_HPP

#include <type_traits>

namespace scran_aggregate {

template<typename Input_>
std::remove_cv_t<std::remove_reference_t<Input_> > I(Input_ x) {
    return x;
}

}

#endif
