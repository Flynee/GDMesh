
#ifdef _DEBUG
#define _ITERATOR_DEBUG_LEVEL 2	// speedup iterator operations while debugging
#endif

#include <iostream>

#include "test/test.h"

int main() {
    std::min(1, 2);
    test_start();
}