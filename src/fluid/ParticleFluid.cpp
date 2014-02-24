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
#include "ParticleFluid.h"
#include "../engine/DynamicsWorld.h"

using namespace reactphysics3d;

// Initialization of static variables
const decimal ParticleFluid::DEFAULT_PARTICLE_MASS = decimal(0.0002);
const decimal ParticleFluid::DEFAULT_FLUID_REST_DENSITY = decimal(1);
const decimal ParticleFluid::DEFAULT_GAS_STIFFNESS = decimal(0.5);
const decimal ParticleFluid::DEFAULT_VISCOSITY = decimal(0.2);

// Constructor
ParticleFluid::ParticleFluid(const ParticleFluidInfo& particleFluidInfo)
      : mPosition(particleFluidInfo.position), mDimension(particleFluidInfo.dimension),
        mNbParticles(0), mIsActive(true), mMassParticle(DEFAULT_PARTICLE_MASS),
        mRestDensity(DEFAULT_FLUID_REST_DENSITY), mGasStiffness(DEFAULT_GAS_STIFFNESS),
        mViscosity(DEFAULT_VISCOSITY), mIsGravityEnabled(true), mIsNotSimulatedYet(true) {

}

// Destructor
ParticleFluid::~ParticleFluid() {

}

// Create and add a new particle into the fluid
void ParticleFluid::createParticle(const Vector3& position, const Vector3& velocity) {

    // Create a new particle
    FluidParticle particle(position, velocity);

    // Add the particle into the fluid
    mParticles.push_back(particle);

    mNbParticles++;
}

// Remove all the particles of the fluid
void ParticleFluid::removeAllParticles() {
    mParticles.clear();
    mNbParticles = 0;
}
