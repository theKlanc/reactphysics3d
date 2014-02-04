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
const uint FluidSolver::GRID_DIMENSION = 16;
const decimal FluidSolver::PARTICLE_MASS = decimal(0.02);
const uint FluidSolver::BLOCK_NB_CELLS = 2;
const decimal FluidSolver::POLY6_FACTOR = decimal(315.0) /
                                          decimal(64.0 * PI * pow(SPH_SUPPORT_RADIUS, 9));
const decimal FluidSolver::SPIKY_FACTOR = decimal(45.0) / decimal(PI * pow(SPH_SUPPORT_RADIUS, 6));

// Constructor
FluidSolver::FluidSolver(std::set<ParticleFluid*>& fluids, Vector3& gravity, decimal& timestep)
            : mFluids(fluids), mGravity(gravity), mTimestep(timestep) {

}

// Destructor
FluidSolver::~FluidSolver() {

}

// Solve the fluid simulation
void FluidSolver::solve() {

    PROFILE("FluidSolver::solve()");

    // For all fluid of the world
    for (std::set<ParticleFluid*>::iterator it = mFluids.begin(); it != mFluids.end(); ++it) {

        ParticleFluid* fluid = *it;

        // If the fluid is not active, do not simulate it
        if (!fluid->mIsActive) continue;

        // Compute the number of blocks of particles
        uint32 nbBlocksLine = GRID_DIMENSION / BLOCK_NB_CELLS;
        mNbBlocks = nbBlocksLine * nbBlocksLine * nbBlocksLine;

        // Allocate memory for the z-index sorted particles array
        // TODO : Use the rp3d memory allocator for this
        mSortedZIndexParticles.clear();
        mBlockParticles = new BlockParticles[mNbBlocks];
        mDensities = new decimal[fluid->getNbParticles()];
        mForces = new Vector3[fluid->getNbParticles()];
        assert(mBlockParticles != NULL);
        assert(mDensities != NULL);
        assert(mForces != NULL);

        // Compute the z-index of the particles and sort them according to this z-index
        computeZIndexAndSortParticles(fluid);

        // Compute the blocks of particles on the grid
        computeBlocks(fluid);

        // Compute the density at the particles location
        computeDensity(fluid);

        // Initialize the velocity of particles
        initParticlesVelocity(fluid);

        // Compute the forces on the particles
        computeForcesOnParticles(fluid);

        // Update the position and velocity of the particles
        updateParticlesPosition(fluid);

        // TODO : Remove this
        // Collision detection and response
        computeCollisionDetection(fluid);

        // Release allocated memory
        delete[] mBlockParticles;
        delete[] mDensities;
        delete[] mForces;
    }
}

// Compute the z-index of the particles and sort them according to this z-index
void FluidSolver::computeZIndexAndSortParticles(ParticleFluid* fluid) {

    PROFILE("FluidSolver::computeZIndexAndSortParticles()");

    assert(mSortedZIndexParticles.empty());

    const Vector3 startGrid = fluid->mPosition - decimal(0.5) * fluid->mDimension;
    const decimal cellSize = fluid->mDimension.x / decimal(GRID_DIMENSION);

    // For each particle of the fluid
    for (uint i=0; i<fluid->mNbParticles; i++) {

        // Get the particle pointer
        FluidParticle* particle = &fluid->mParticles[i];

        // Compute its (x,y,z) grid coordinate
        uint32 x = static_cast<uint32>((fluid->mParticles[i].position.x - startGrid.x) / cellSize);
        uint32 y = static_cast<uint32>((fluid->mParticles[i].position.y - startGrid.y) / cellSize);
        uint32 z = static_cast<uint32>((fluid->mParticles[i].position.z - startGrid.z) / cellSize);

        // Compute the corresponding z-index coordinate
        particle->zIndex = computeZIndex(x, y, z);

        // Add the particle pointer into the array to be sorted
        mSortedZIndexParticles.push_back(particle);
    }

    // Sort the particles array according to their z-index using radix sort
    radixSort<FluidParticle>(mSortedZIndexParticles, fluid->getNbParticles());
}

// Compute the particles blocks on the grid
void FluidSolver::computeBlocks(ParticleFluid* fluid) {

    PROFILE("FluidSolver::computeBlocks()");

    const decimal cellSize = fluid->mDimension.x / decimal(GRID_DIMENSION);
    const uint32 nbBlocks = GRID_DIMENSION / BLOCK_NB_CELLS;
    const Vector3 startGrid = fluid->mPosition - decimal(0.5) * fluid->mDimension;

    // For each particles of the sorted z-index array
    for (uint i=0; i<fluid->getNbParticles(); i++) {

        // Get the particle
        const FluidParticle* particle = mSortedZIndexParticles[i];

        // Compute its (x,y,z) grid coordinate
        // TODO : Maybe we can store this in the particle at the previous step
        uint32 x = static_cast<uint32>((particle->position.x - startGrid.x) / cellSize);
        uint32 y = static_cast<uint32>((particle->position.y - startGrid.y) / cellSize);
        uint32 z = static_cast<uint32>((particle->position.z - startGrid.z) / cellSize);

        // Compute the (x,y,z) block coordinates
        uint32 xBlock = x / BLOCK_NB_CELLS;
        uint32 yBlock = y / BLOCK_NB_CELLS;
        uint32 zBlock = z / BLOCK_NB_CELLS;

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
void FluidSolver::computeDensity(ParticleFluid* fluid) {

    PROFILE("FluidSolver::computeDensity()");

    /*
    timeval timeValue;
    gettimeofday(&timeValue, NULL);
    double timeStart = (timeValue.tv_sec + (timeValue.tv_usec / 1000000.0));
    */

    //std::cout << "---------- COMPUTE DENSITY ---------" << std::endl;

    const uint32 nbBlocks = GRID_DIMENSION / BLOCK_NB_CELLS;
    const uint32 nbBlocksSquare = nbBlocks * nbBlocks;
    uint nbParticles = 0;
    uint nbLoop = 0;

    // For each block
    for (uint32 i=0; i<mNbBlocks; i++) {

        PROFILE("FluidSolver::computeDensity()::ForEachBlock");

        // Compute the (x,y,z) block coordinates
        const int xBlock = (i % nbBlocksSquare) % nbBlocks;
        const int yBlock = (i % nbBlocksSquare) / nbBlocks;
        const int zBlock = i / nbBlocksSquare;

        //std::cout << "_Block " << i << " with " << mBlockParticles[i].nbParticles << " particles" << std::endl;

        // For each particle of the block
        for (uint32 p=0; p<mBlockParticles[i].nbParticles; p++) {

            PROFILE("FluidSolver::computeDensity()::ForEachParticleOfBlock");

            const uint32 indexParticle = mBlockParticles[i].firstParticle + p;
            const FluidParticle* particle = mSortedZIndexParticles[indexParticle];

            //std::cout << "___Particle" << indexParticle << std::endl;

            decimal sumDensityFactor = decimal(0.0);
            nbParticles++;

            // For each of the 27 neighboring block
            for (int b=0; b<27; b++) {

                PROFILE("FluidSolver::computeDensity()::ForEachNeighboringBlock");

                // Compute the (x,y,z) coordinates of the neighbor block
                const int xDiff = ((b % 9) % 3) - 2 + xBlock;
                const int yDiff = ((b % 9) / 3) - 2 + yBlock;
                const int zDiff = (b / 9) - 2 + zBlock;

                // If the neighboring block is out of the grid
                if (xDiff < 0 || xDiff > nbBlocks ||
                    yDiff < 0 || yDiff > nbBlocks ||
                    zDiff < 0 || zDiff > nbBlocks) continue;

                uint32 indexBlock = zDiff * nbBlocksSquare + yDiff * nbBlocks + xDiff;

                //std::cout << "______Neighbor block " << b << " with " << mBlockParticles[indexBlock].nbParticles << " particles" << std::endl;

                // For each particle of the neighboring block
                for (uint32 q=0; q < mBlockParticles[indexBlock].nbParticles; q++) {
                    nbLoop++;

                    PROFILE("FluidSolver::computeDensity()::ForEachParticleinNeighboringBlock");

                    uint32 indexNeighbor = mBlockParticles[indexBlock].firstParticle + q;
                    const FluidParticle* neighbor = mSortedZIndexParticles[indexNeighbor];

                    //std::cout << "_________Neighbor particle " << indexNeighbor << std::endl;

                    // Compute the square distance from particle to its neighbor
                    decimal distSquare = (particle->position - neighbor->position).lengthSquare();

                    // If the neighbor is not inside the SPH support radius, we go to the next one
                    if (distSquare > SPH_SUPPORT_RADIUS_SQUARE) continue;

                    // Increase the density value
                    const decimal value = SPH_SUPPORT_RADIUS_SQUARE - distSquare;
                    assert(value >= decimal(0.0));
                    sumDensityFactor += value * value * value;
                    //mDensities[indexParticle] += kernelPoly6(distSquare);
                }
            }

            mDensities[indexParticle] = sumDensityFactor * fluid->mMassParticle * POLY6_FACTOR;

            //std::cout << "Particle " << indexParticle << ", density = " << mDensities[indexParticle] << std::endl;
            assert(mDensities[indexParticle] > MACHINE_EPSILON);
        }
    }

    /*
    //std::cout << "Nb Particles density : " << nbParticles << std::endl;
    timeval timeValue2;
    gettimeofday(&timeValue2, NULL);
    double timeEnd = (timeValue.tv_sec + (timeValue.tv_usec / 1000000.0));
    std::cout << "End time : " << timeEnd << std::endl;
    std::cout << "Density time : " << (timeEnd - timeStart) << std::endl;
    */
    //std::cout << "Nb loops : " << nbLoop << std::endl;
}

// Initialize the velocity of particles
void FluidSolver::initParticlesVelocity(ParticleFluid* fluid) {

    const Vector3 gravity = fluid->mIsGravityEnabled ? mGravity : Vector3(0, 0, 0);

    // For each particle of the fluid
    for (uint i=0; i<fluid->mNbParticles; i++) {

        FluidParticle* particle = mSortedZIndexParticles[i];

        // If the particle has not been simulated yet
        if (particle->isNotSimulatedYet) {

            decimal test = mDensities[i];   // TODO : DELETE THIS

            // Compute the initial velocity of the particle (for Leap-Frog integration scheme)
            Vector3 initAcceleration = gravity + fluid->mExternalForce / mDensities[i];
            particle->velocity -= decimal(0.5) * mTimestep * initAcceleration;
            particle->velocityEvaluation = particle->velocity;

            particle->isNotSimulatedYet = false;
        }
    }
}

// Compute the forces on the particles and update their positions
void FluidSolver::computeForcesOnParticles(ParticleFluid* fluid) {

    PROFILE("FluidSolver::computeForcesOnParticles()");

    const uint32 nbBlocks = GRID_DIMENSION / BLOCK_NB_CELLS;
    const uint32 nbBlocksSquare = nbBlocks * nbBlocks;

    // For each block
    for (uint32 i=0; i<mNbBlocks; i++) {

        // Compute the (x,y,z) block coordinates
        const int xBlock = (i % nbBlocksSquare) % nbBlocks;
        const int yBlock = (i % nbBlocksSquare) / nbBlocks;
        const int zBlock = i / nbBlocksSquare;

        // For each particle of the block
        for (uint32 p=0; p<mBlockParticles[i].nbParticles; p++) {

            const uint32 indexParticle = mBlockParticles[i].firstParticle + p;
            FluidParticle* particle = mSortedZIndexParticles[indexParticle];

            // Compute the pression at the current particle position
            decimal pressionParticle = fluid->mGasStiffness *
                                       (mDensities[indexParticle] - fluid->mRestDensity);

            Vector3 sumForce(0, 0, 0);

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

                    uint32 indexNeighbor = mBlockParticles[indexBlock].firstParticle + q;
                    const FluidParticle* neighbor = mSortedZIndexParticles[indexNeighbor];

                    // Compute the square distance from particle to its neighbor
                    const Vector3 r = particle->position - neighbor->position;
                    decimal distSquare = r.lengthSquare();

                    // If the neighbor is not inside the SPH support radius, we go to the next one
                    if (distSquare > SPH_SUPPORT_RADIUS_SQUARE) continue;

                    // Compute the distance between the particle and its neighbor
                    decimal distance = std::sqrt(distSquare);

                    // Compute the pression at the current neihbor position
                    decimal pressionNeighbor = fluid->mGasStiffness *
                                               (mDensities[indexNeighbor] - fluid->mRestDensity);

                    // Increase the pression force
                    const decimal value = SPH_SUPPORT_RADIUS - distance;
                    assert(value >= decimal(0.0));
                    const decimal factor = distance < MACHINE_EPSILON ? decimal(0.0) : decimal(1.0) / distance;
                    Vector3 pressionForce = decimal(0.5) * (pressionParticle + pressionNeighbor) * r * factor * value * value;

                    // Increase the viscosity force
                    Vector3 viscosityForce = (neighbor->velocityEvaluation - particle->velocityEvaluation) * value;

                    sumForce += (pressionForce + fluid->mViscosity * viscosityForce) / mDensities[indexNeighbor];
                }
            }

            // Compute the total force on the particle
            mForces[indexParticle] = fluid->mMassParticle * SPIKY_FACTOR * sumForce;
            //mForces[indexParticle] = pressionForce + fluid->mViscosity * viscosityForce +
            //                         fluid->mExternalForce;
        }
    }
}

// Update the particles positions and velocities
void FluidSolver::updateParticlesPosition(ParticleFluid* fluid) {

    PROFILE("FluidSolver::updateParticlesPosition()");

    const Vector3 gravity = fluid->mIsGravityEnabled ? mGravity : Vector3(0, 0, 0);

    const uint32 nbBlocks = GRID_DIMENSION / BLOCK_NB_CELLS;
    const uint32 nbBlocksSquare = nbBlocks * nbBlocks;

    // For each block
    for (uint32 i=0; i<mNbBlocks; i++) {

        // Compute the (x,y,z) block coordinates
        const int xBlock = (i % nbBlocksSquare) % nbBlocks;
        const int yBlock = (i % nbBlocksSquare) / nbBlocks;
        const int zBlock = i / nbBlocksSquare;

        // For each particle of the block
        for (uint32 p=0; p<mBlockParticles[i].nbParticles; p++) {

            const uint32 indexParticle = mBlockParticles[i].firstParticle + p;
            FluidParticle* particle = mSortedZIndexParticles[indexParticle];

            // Compute the acceleration of the particle and add gravity
            // TODO : Uncomment this
            const Vector3 acceleration = (mForces[indexParticle] / mDensities[indexParticle]) + gravity;
            //const Vector3 acceleration = gravity;

            Vector3 oldVelocity = particle->velocity;

            // Integrate the velocity of the particle using Leap-Frog integration
            particle->velocity += acceleration * mTimestep;

            // Compute the velocity used for evaluation at velocity mid-step
            particle->velocityEvaluation += decimal(0.5) * (oldVelocity + particle->velocity);

            // Integrate the position of the particle using Leap-Frog integration
            particle->position += particle->velocity * mTimestep;

        }
    }
}

// TODO : Delete this
void FluidSolver::computeCollisionDetection(ParticleFluid* fluid) {

    PROFILE("FluidSolver::computeCollisionDetection()");

    Vector3 fluidExtent = decimal(0.5) * fluid->mDimension;

    // For each particle
    for (uint i=0; i<fluid->getNbParticles(); ++i) {

        FluidParticle& particle = fluid->mParticles[i];
        Vector3 particleExtent = particle.position - fluid->mPosition;
        Vector3 particleAbsExtent(std::abs(particleExtent.x),
                                  std::abs(particleExtent.y),
                                  std::abs(particleExtent.z));

        Vector3 initParticlePos = particle.position;    // TODO : DELETE THIS


        // If the particle is inside the fluid box, go to the next one
        if (particleAbsExtent.x < fluidExtent.x && particleAbsExtent.y < fluidExtent.y &&
            particleAbsExtent.z < fluidExtent.z) continue;

        // Compute the contact point
        Vector3 localParticlePos = particle.position - fluid->mPosition;
        Vector3 contactPoint(std::min(fluidExtent.x, std::max(-fluidExtent.x, localParticlePos.x)),
                             std::min(fluidExtent.y, std::max(-fluidExtent.y, localParticlePos.y)),
                             std::min(fluidExtent.z, std::max(-fluidExtent.z, localParticlePos.z)));
        Vector3 worldContactPoint = contactPoint + fluid->mPosition;
        Vector3 diff = worldContactPoint - particle.position;
        decimal depth = diff.length();
        Vector3 normal = (depth > MACHINE_EPSILON) ? diff.getUnit() : -particle.velocityEvaluation.getUnit();

        // Collision response
        decimal restitution = decimal(0.01);
        particle.position += (depth + 0.001) * normal;
        particle.velocity -= (decimal(1.0) + restitution * depth / (mTimestep * particle.velocityEvaluation.length())) *
                             (particle.velocityEvaluation.dot(normal)) * normal;

        particleExtent = particle.position - fluid->mPosition;
        particleAbsExtent = Vector3(std::abs(particleExtent.x),
                                    std::abs(particleExtent.y),
                                    std::abs(particleExtent.z));
        assert(particleAbsExtent.x <= fluidExtent.x);
        assert(particleAbsExtent.y <= fluidExtent.y);
        assert(particleAbsExtent.z <= fluidExtent.z);
    }
}
