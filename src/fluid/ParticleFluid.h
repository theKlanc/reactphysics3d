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

#ifndef REACTPHYSICS3D_FLUID_H
#define REACTPHYSICS3D_FLUID_H

// Libraries
#include "../mathematics/Vector3.h"
#include <vector>

/// Namespace ReactPhysics3D
namespace reactphysics3d {

// TODO : Add get/set methods for the different variables of the fluid

// Structure FluidParticle
/**
 * This structure represents particle of a fluid.
 */
struct FluidParticle {

    // -------------------- Attributes -------------------- //

    /// Position of the particle
    Vector3 position;

    /// Velocity of the particle
    Vector3 velocity;

    /// True if the particle has not been simulated yet
    bool isNotSimulatedYet;

    // -------------------- Methods -------------------- //

    /// Constructor
    FluidParticle(const Vector3& initPosition, const Vector3& initVelocity)
                 :position(initPosition), velocity(initVelocity), isNotSimulatedYet(true) {

    }

};

// Structure ParticleFluidInfo
/**
 * This structure is used to gather the information needed to create a particle fluid.
 */
struct ParticleFluidInfo {

    /// Dimension of the fluid grid
    Vector3 dimension;

    /// Position of the center of the fluid grid
    Vector3 position;

    /// Construction
    ParticleFluidInfo(const Vector3& initPosition, const Vector3& initDimension)
                     : position(initPosition), dimension(initDimension) {

    }
};

// Class ParticleFluid
/**
 * This class represents a fluid made of several particles.
 */
class ParticleFluid {

    private:
        // -------------------- Attributes -------------------- //

        /// Default particle mass (in kg)
        const static decimal DEFAULT_PARTICLE_MASS;

        /// Default fluid rest density in kg/m^3 (density of water)
        const static decimal DEFAULT_FLUID_REST_DENSITY;

        /// Default gas stiffness constant
        const static decimal DEFAULT_GAS_STIFFNESS;

        /// Default fluid viscosity coefficient
        const static decimal DEFAULT_VISCOSITY;

        // -------------------- Attributes -------------------- //

        /// World-space fluid dimensions in the x, y, z directions
        Vector3 mDimension;

        /// Center position of the fluid in world-space coordinates
        Vector3 mPosition;

        /// Number of partices in the fluid
        uint32 mNbParticles;

        /// True if the simulation of this fluid is enabled, false otherwise
        bool mIsActive;

        /// Array with all the partices of the fluid
        // TODO : Do not use a std::vector here but replace with C-array using rp3d memmory
        //        allocator
        std::vector<FluidParticle> mParticles;

        /// Mass of each particle in the fluid
        decimal mMassParticle;

        /// Rest density of the fluid
        decimal mRestDensity;

        /// Gas stiffness constant
        decimal mGasStiffness;

        /// Fluid viscosity coefficient
        decimal mViscosity;

        /// True if the gravity needs to be applied to this fluid
        bool mIsGravityEnabled;

        /// Current external force acting on the fluid
        Vector3 mExternalForce;

        /// True if this fluid has not been simulated yet
        bool mIsNotSimulatedYet;

    public:

        // -------------------- Methods -------------------- //

        /// Constructor
        ParticleFluid(const ParticleFluidInfo& particleFluidInfo);

        /// Destructor
        ~ParticleFluid();

        /// Return the number of particles in the fluid
        uint32 getNbParticles() const;

        /// Create and add a new particle into the fluid
        void createParticle(const Vector3& position, const Vector3& velocity);

        /// Remove all the particles of the fluid
        void removeAllParticles();

        /// Return a given particle of the fluid
        const FluidParticle& getParticle(uint index);

        // -------------------- Friendship -------------------- //

        friend class FluidSolver;

};

// Return the number of particles in the fluid
inline uint32 ParticleFluid::getNbParticles() const {
    return mNbParticles;
}

// Return a given particle of the fluid
inline const FluidParticle& ParticleFluid::getParticle(uint index) {
    assert(index >= 0 && index < mNbParticles);
    return mParticles[index];
}

}

#endif
