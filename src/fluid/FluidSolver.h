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

#ifndef REACTPHYSICS3D_FLUID_SOLVER_H
#define REACTPHYSICS3D_FLUID_SOLVER_H

// Libraries
#include "Fluid.h"
#include <set>
#include <map>
#include <cassert>

// TODO : Make possible to use non-cubic fluid area (now grid only allows same x,y,z dimensions)

/// Namespace ReactPhysics3D
namespace reactphysics3d {

// Structure BlockParticles
/**
 * This structure represent a block of particles in the grid.
 */
struct BlockParticles {

    /// Index of the first particle in the array of z-index
    uint32 firstParticle;

    /// Number of particles in the block
    uint32 nbParticles;

    // Constructor
    BlockParticles() {
        firstParticle = std::numeric_limits<uint32>::max();
        nbParticles = 0;
    }

};

// Class FluidSolver
/**
 * This class is used to simulate fluids using the Smoothed Particle Hydrodynamics (SPH)
 * technique.
 */
class FluidSolver {

    private:

        // -------------------- Constants -------------------- //

        /// Global support radius (in meters) of the SPH solver
        static const decimal SPH_SUPPORT_RADIUS;

        /// Square of the SPH support radius
        static const decimal SPH_SUPPORT_RADIUS_SQUARE;

        /// Number of cells in the grid
        static const uint GRID_DIMENSION;

        /// Number of cells in the side of a block (matching the size of the support radius)
        static const uint BLOCk_NB_CELLS;

        /// Mass of a particle in the fluid (in kg)
        static const decimal PARTICLE_MASS;

        /// Factor for the Poly6 smoothing kernel
        static const decimal POLY6_FACTOR;

        // -------------------- Attributes -------------------- //

        /// Reference to the set of all fluids of the world
        std::set<Fluid*> mFluids;

        /// Map the z-index of a particle with its index in the fluid
        std::map<uint32, uint32> mMapZIndexToParticleIndex;

        /// Array with the sorted z-index of the particles
        uint32* mSortedZIndexParticles;

        /// Array with the densities at particles locations
        decimal* mDensities;

        /// Array with the blocks of particles
        BlockParticles* mBlockParticles;

        /// Number of blocks
        uint32 mNbBlocks;


        // -------------------- Methods -------------------- //

        /// Compute the z-index of the particles and sort them according to this z-index
        void computeZIndexAndSortParticles(Fluid* fluid);

        /// Compute the linear z-index coordinate given a 3D grid coordinate
        uint32 computeZIndex(uint32 x, uint32 y, uint32 z) const;

        /// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit.
        uint32 expandBits(uint32 x) const;

        /// Compute the particles blocks on the grid
        void computeBlocks(Fluid *fluid);

        /// Compute the density at the particles locations
        void computeDensity(Fluid* fluid);

        /// Smoothing kernel Poly6 for the SPH simulation
        decimal kernelPoly6(decimal distanceSquare);


    public:

        // -------------------- Methods -------------------- //

        /// Constructor
        FluidSolver(std::set<Fluid*>& fluids);

        /// Destructor
        ~FluidSolver();

        /// Solve the fluid simulation
        void solve();

};

// Compute the linear z-index coordinate given a 3D grid coordinate
inline uint32 FluidSolver::computeZIndex(uint32 x, uint32 y, uint32 z) const {

    // TODO : Instead of computing the z-index each time, store it in a lookup-table.

    uint32 xx = expandBits((uint32)x);
    uint32 yy = expandBits((uint32)y);
    uint32 zz = expandBits((uint32)z);
    return xx * 4 + yy * 2 + zz;
}

// Expands a 10-bit integer into 30 bits by inserting 2 zeros after each bit.
inline uint32 FluidSolver::expandBits(uint32 x) const {
    x = (x * 0x00010001u) & 0xFF0000FFu;
    x = (x * 0x00000101u) & 0x0F00F00Fu;
    x = (x * 0x00000011u) & 0xC30C30C3u;
    x = (x * 0x00000005u) & 0x49249249u;
    return x;
}

// Smoothing kernel Poly6 for the SPH simulation
inline decimal FluidSolver::kernelPoly6(decimal distanceSquare) {
    decimal value = SPH_SUPPORT_RADIUS_SQUARE - distanceSquare;
    assert(value >= decimal(0.0));
    return POLY6_FACTOR * value * value * value;
}

}

#endif
