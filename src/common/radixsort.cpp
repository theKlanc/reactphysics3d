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

// Libraries
#include "radixsort.h"

using namespace reactphysics3d;

// Function that sorts an array of integer using the radix sort algorithm.
void reactphysics3d::radixSort(uint32 array[], uint32 nbElements) {

    // Find the maximum element in the array
    uint32 maxElement = array[0];
    for (uint32 i = 1; i < nbElements; i++) {
        if (array[i] > maxElement) {
            maxElement = array[i];
        }
    }

    // Perform a counting sort pass for every digit
    for (uint32 exp=1; (maxElement / exp) > 0; exp *= 10) {
        countingSort(array, nbElements, exp);
    }
}

// Function that sorts an array of integers according to the digit represented by "exp"
void reactphysics3d::countingSort(uint32 array[], uint32 nbElements, uint32 exp) {
    uint32* outputArray = new uint32[nbElements];
    uint32 count[10] = {0};
    int i;

    // Compute the histogram of the digit occurences in count[]
    for (i=0; i<nbElements; i++) {
        count[(array[i]/exp) % 10]++;
    }

    // Change count[] so that it contains the actual position of the digits in the output array
    for (i=1; i<10; i++) {
        count[i] += count[i - 1];
    }

    // Fill-in the output array
    for (i=nbElements-1; i >= 0; i--) {
        outputArray[count[(array[i]/exp) % 10] - 1] = array[i];
        count[(array[i]/exp) % 10]--;
    }

    // Copy the output array into the input array
    for (i=0; i<nbElements; i++) {
        array[i] = outputArray[i];
    }

    delete[] outputArray;
}
