#ifndef UTILS_H
#define UTILS_H

// This generates a different initial value for each vector() that might store results.
// The aim is to check that our functions ignore the initial contents of each output array.
// Otherwise, we might miss subtle bugs from an incorrect assumption of default-initialized contents.
// We make the values vary across calls rather than using a constant,
// as this allows us to catch bugs via tests that check for consistency between calls.
inline int initial_value() {
    static int counter = 0;
    if (counter == 255) {
        counter = 1;
    } else {
        ++counter;
    }
    return counter;
}

#endif
