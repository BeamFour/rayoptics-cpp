// C++ port of org.redukti.util.ArrayIndex2D
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_UTIL_ARRAYINDEX2D_H
#define REDUKTI_UTIL_ARRAYINDEX2D_H

namespace redukti::util {

/** Row-major index arithmetic for a 2D array held in one flat vector. */
class ArrayIndex2D {
public:
    int rowSize;
    int colSize;

    ArrayIndex2D(int rowSize_, int colSize_) : rowSize(rowSize_), colSize(colSize_) {}

    int i(int row, int col) const { return colSize * row + col; }
};

} // namespace redukti::util

#endif // REDUKTI_UTIL_ARRAYINDEX2D_H
