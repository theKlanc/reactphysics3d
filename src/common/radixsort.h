/********************************************************************************
* ReactPhysics3D physics library, http://code.google.com/p/reactphysics3d/      *
* Copyright (c) 2010-2013 Daniel Chappuis                                       *
*********************************************************************************
*                                                                               *
* This software is provided 'as-is', without any express or implied warranty.   *
* In no event will the authors be held liable for any damages arising from the  *
* use of this software.                                                         *
*                                                                               *
* Permission is granted to anyone to use this software for any purpose,         *
* including commercial applications, and to alter it and redistribute it        *
* freely, subject to the following restrictions:                                *
*                                                                               *
* 1. The origin of this software must not be misrepresented; you must not claim *
*    that you wrote the original software. If you use this software in a        *
*    product, an acknowledgment in the product documentation would be           *
*    appreciated but is not required.                                           *
*                                                                               *
* 2. Altered source versions must be plainly marked as such, and must not be    *
*    misrepresented as being the original software.                             *
*                                                                               *
* 3. This notice may not be removed or altered from any source distribution.    *
*                                                                               *
********************************************************************************/

#ifndef REACTPHYSICS3D_RADIXSORT_H
#define REACTPHYSICS3D_RADIXSORT_H

// Libraries
#include "../configuration.h"
#include <vector>
#include <iostream>

namespace reactphysics3d {

/// Function that sorts an array of objects using the radix sort algorithm.
template<typename T>
void radixSort(std::vector<T*>& array, uint32 nbElements) {

    assert(array.size() == nbElements);

    // Find the maximum key in the array
    uint32 maxKey = array[0]->key();
    for (uint32 i = 1; i < nbElements; i++) {
        if (array[i]->key() > maxKey) {
            maxKey = array[i]->key();
        }
    }

    // Perform a counting sort pass for every digit
    for (uint32 exp=1; (maxKey / exp) > 0; exp *= 10) {
        countingSort(array, nbElements, exp);
    }
}

/// Function that sorts an array of integers according to the digit represented by "exp"
template<typename T>
void countingSort(std::vector<T*>& array, uint32 nbElements,
                                  uint32 exp) {

    assert(array.size() == nbElements);

    std::vector<T*> outputArray;
    uint32 count[10] = {0};
    int i;

    // Compute the histogram of the digit occurences in count[]
    for (i=0; i<nbElements; i++) {
        count[(array[i]->key() / exp) % 10]++;
        outputArray.push_back(NULL);
    }

    // Change count[] so that it contains the actual position of the digits in the output array
    for (i=1; i<10; i++) {
        count[i] += count[i - 1];
    }

    // Fill-in the output array
    for (i=nbElements-1; i >= 0; i--) {
        outputArray[count[(array[i]->key() / exp) % 10] - 1] = array[i];
        count[(array[i]->key() / exp) % 10]--;
    }

    // Copy the output array into the input array
    for (i=0; i<nbElements; i++) {
        array[i] = outputArray[i];
    }
}

}

#endif
