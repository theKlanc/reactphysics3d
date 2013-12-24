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
#include "FluidSolver.h"
#include "../common/radixsort.h"
#include <cmath>


// TODO : Add assert() verification in the code

using namespace reactphysics3d;

// Initialization of static variables
const decimal FluidSolver::SPH_SUPPORT_RADIUS = decimal(0.02);
const decimal FluidSolver::SPH_SUPPORT_RADIUS_SQUARE = SPH_SUPPORT_RADIUS * SPH_SUPPORT_RADIUS;
const uint FluidSolver::GRID_DIMENSION = 64;
const decimal FluidSolver::PARTICLE_MASS = decimal(0.02);
const uint FluidSolver::BLOCk_NB_CELLS = 4;
const decimal FluidSolver::POLY6_FACTOR = decimal(315.0) /
                                          decimal(64.0 * PI * pow(SPH_SUPPORT_RADIUS, 9));

// Constructor
FluidSolver::FluidSolver(std::set<Fluid*>& fluids) : mFluids(fluids) {

}

// Destructor
FluidSolver::~FluidSolver() {

}

// Solve the fluid simulation
void FluidSolver::solve() {

    // For all fluid of the world
    for (std::set<Fluid*>::iterator it = mFluids.begin(); it != mFluids.end(); ++it) {

        Fluid* fluid = *it;

        // Compute the number of blocks of particles
        uint32 nbBlocksLine = GRID_DIMENSION / BLOCk_NB_CELLS;
        mNbBlocks = nbBlocksLine * nbBlocksLine * nbBlocksLine;

        // Allocate memory for the z-index sorted particles array
        // TODO : Use the rp3d memory allocator for this
        mSortedZIndexParticles = new uint32[fluid->getNbParticles()];
        mBlockParticles = new BlockParticles[mNbBlocks];
        mDensities = new decimal[fluid->getNbParticles()];
        assert(mSortedZIndexParticles != NULL);
        assert(mBlockParticles != NULL);

        // Compute the z-index of the particles and sort them according to this z-index
        computeZIndexAndSortParticles(fluid);

        // Compute the blocks of particles on the grid
        computeBlocks(fluid);

        // Compute the density at the particles location
        computeDensity(fluid);

        // Release allocated memory
        delete[] mSortedZIndexParticles;
        delete[] mBlockParticles;
        delete[] mDensities;
    }
}

// Compute the z-index of the particles and sort them according to this z-index
void FluidSolver::computeZIndexAndSortParticles(Fluid* fluid) {

    const Vector3 startGrid = fluid->mPosition - decimal(0.5) * fluid->mDimension;
    const decimal cellSize = fluid->mDimension.x / decimal(GRID_DIMENSION);

    // For each particle of the fluid
    for (uint i=0; i<fluid->getNbParticles(); i++) {

        // Compute its (x,y,z) grid coordinate
        uint32 x = static_cast<uint32>((fluid->mParticles[i].position.x - startGrid.x) / cellSize);
        uint32 y = static_cast<uint32>((fluid->mParticles[i].position.y - startGrid.y) / cellSize);
        uint32 z = static_cast<uint32>((fluid->mParticles[i].position.z - startGrid.z) / cellSize);

        // Compute the corresponding z-index coordinate
        uint32 zIndex = computeZIndex(x, y, z);

        // Add the mapping between the z-index and the particle index
        mMapZIndexToParticleIndex.insert(std::pair<uint32,uint32>(zIndex, i));

        // Add the z-index into the z-index array
        mSortedZIndexParticles[i] = zIndex;
    }

    // Sort the z-index array using radix sort
    radixSort(mSortedZIndexParticles, fluid->getNbParticles());
}

// Compute the particles blocks on the grid
void FluidSolver::computeBlocks(Fluid* fluid) {

    const decimal cellSize = fluid->mDimension.x / decimal(GRID_DIMENSION);
    const uint32 nbBlocks = GRID_DIMENSION / BLOCk_NB_CELLS;
    const Vector3 startGrid = fluid->mPosition - decimal(0.5) * fluid->mDimension;

    // For each particles of the sorted z-index array
    for (uint i=0; i<fluid->getNbParticles(); i++) {

        // Get the particle
        uint32 indexParticle = mMapZIndexToParticleIndex.at(mSortedZIndexParticles[i]);
        const FluidParticle& particle = fluid->mParticles[indexParticle];

        // Compute its (x,y,z) grid coordinate
        // TODO : Maybe we can store this in the particle at the previous step
        uint32 x = static_cast<uint32>((particle.position.x - startGrid.x) / cellSize);
        uint32 y = static_cast<uint32>((particle.position.y - startGrid.y) / cellSize);
        uint32 z = static_cast<uint32>((particle.position.z - startGrid.z) / cellSize);

        // Compute the (x,y,z) block coordinates
        uint32 xBlock = x / BLOCk_NB_CELLS;
        uint32 yBlock = y / BLOCk_NB_CELLS;
        uint32 zBlock = z / BLOCk_NB_CELLS;

        // Compute the index of the block in the blocks array
        uint32 indexBlock = zBlock * (nbBlocks * nbBlocks) + yBlock * nbBlocks + xBlock;

        // Compute the starting z-index of the block
        if (i < mBlockParticles[indexBlock].firstParticle) {
            mBlockParticles[indexBlock].firstParticle = i;
        }

        // Increment the number of particles in the block
        mBlockParticles[indexBlock].nbParticles++;
    }
}

// Compute the density at the particles locations
void FluidSolver::computeDensity(Fluid* fluid) {

    const uint32 nbBlocks = GRID_DIMENSION / BLOCk_NB_CELLS;
    const uint32 nbBlocksSquare = nbBlocks * nbBlocks;

    // For each block
    for (uint32 i=0; i<mNbBlocks; i++) {

        // Compute the (x,y,z) block coordinates
        const int xBlock = (i % nbBlocksSquare) % nbBlocks;
        const int yBlock = (i % nbBlocksSquare) / nbBlocks;
        const int zBlock = i / nbBlocksSquare;

        // For each particle of the block
        for (uint32 p=0; p<mBlockParticles[i].nbParticles; p++) {

            const uint32 index = mBlockParticles[i].firstParticle + p;
            const uint32 indexParticle = mMapZIndexToParticleIndex[mSortedZIndexParticles[index]];
            const FluidParticle& particle = fluid->mParticles[indexParticle];

            mDensities[index] = decimal(0.0);

            // For each of the 27 neighboring block
            for (int b=0; b<27; b++) {

                // Compute the (x,y,z) coordinates of the neighbor block
                const int xDiff = ((b % 9) % 3) - 2 + xBlock;
                const int yDiff = ((b % 9) / 3) - 2 + yBlock;
                const int zDiff = (b / 9) - 2 + zBlock;

                // If the neighboring block is out of the grid
                if (xDiff < 0 || xDiff > nbBlocks ||
                    yDiff < 0 || yDiff > nbBlocks ||
                    zDiff < 0 || zDiff > nbBlocks) continue;

                uint32 indexBlock = zDiff * nbBlocksSquare + yDiff * nbBlocks + xDiff;

                // For each particle of the neighboring block
                for (uint32 q=0; q < mBlockParticles[indexBlock].nbParticles; q++) {

                    uint32 indexArray = mBlockParticles[indexBlock].firstParticle + q;
                    uint32 indexNeighbor = mMapZIndexToParticleIndex.at(
                                              mSortedZIndexParticles[indexArray]);

                    const FluidParticle& neighbor = fluid->mParticles[indexNeighbor];

                    // Compute the square distance from particle to its neighbor
                    decimal distSquare = (particle.position - neighbor.position).lengthSquare();

                    // If the neighbor is not inside the SPH support radius, we go to the next one
                    if (distSquare > SPH_SUPPORT_RADIUS_SQUARE) continue;

                    // Increase the density value
                    mDensities[index] += fluid->mMassParticle * kernelPoly6(distSquare);
                }
            }

        }
    }
}
