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
#include "Scene.h"
#include <cstdlib>

// Namespaces
using namespace openglframework;

// Constructor
Scene::Scene(Viewer* viewer, const std::string& shaderFolderPath,
             const std::string& meshFolderPath) : mViewer(viewer), mLight0(0),
                               mPhongShader(shaderFolderPath + "phong.vert",
                                            shaderFolderPath + "phong.frag"), mIsRunning(false) {

    // Move the light 0
    mLight0.translateWorld(Vector3(7, 15, 15));

    // Compute the radius and the center of the scene
    float radiusScene = 0.8f;
    openglframework::Vector3 center(0, 0, 0);

    // Set the center of the scene
    mViewer->setScenePosition(center, radiusScene);

    // Gravity vector in the dynamics world
    rp3d::Vector3 gravity(0, rp3d::decimal(-9.81), 0);

    // Time step for the physics simulation
    // TODO : Try different timesteps (check with 1/60 also)
    rp3d::decimal timeStep = 1.0f / 350.0f;

    // Create the dynamics world for the physics simulation
    mDynamicsWorld = new rp3d::DynamicsWorld(gravity, timeStep);

    // Set the number of iterations of the constraint solver
    mDynamicsWorld->setNbIterationsVelocitySolver(15);

    // Create the particle fluid
    createFluid();

    // Create the floor
    openglframework::Vector3 floorPosition(0, 0, 0);
    mFloor = new Box(FLOOR_SIZE, floorPosition, FLOOR_MASS, mDynamicsWorld);

    // The floor must be a static rigid body
    mFloor->getRigidBody()->setType(rp3d::STATIC);

    // Change the material properties of the floor rigid body
    rp3d::Material& material = mFloor->getRigidBody()->getMaterial();
    material.setBounciness(rp3d::decimal(0.3));

    // Create the sphere used to render the particles
    mSphere = new Sphere(PARTICLE_SPHERE_RADIUS, meshFolderPath);

    // Start the simulation
    startSimulation();
}

// Create the fluid
void Scene::createFluid() {

    // Position and dimension of the fluid grid
    const rp3d::Vector3 position(0, 0.0, 0);
    const rp3d::Vector3 dimension(0.24, 0.24, 0.24);

    // Create the fluid info object
    rp3d::ParticleFluidInfo fluidInfo(position, dimension);

    // Create the fluid in the world
    mFluid = mDynamicsWorld->createParticleFluid(fluidInfo);

    uint nbParticles = 700;
    float startDim = dimension.x * 0.5f;
    float vertDim = dimension.y * 0.5f;
    assert(startDim < dimension.x && startDim < dimension.y && startDim < dimension.z);
    /*
    for (uint p=0; p<nbParticles; p++) {

        uint i = rand() % 100 + 1;
        uint j = rand() % 100 + 1;
        uint k = rand() % 100 + 1;
        const rp3d::Vector3 particlePosition(position.x - 0.5f * startDim + (float(i) / 100.0) * startDim,
                                             position.y - 0.5f * vertDim + (float(j) / 100.0) * vertDim,
                                             position.z - 0.5f * startDim + (float(k) / 100.0) * startDim);

        // Create a new particle
        mFluid->createParticle(particlePosition, rp3d::Vector3(0.0, 0.0, 0));
    }
    */
    float epsilon = 0.0001f;
    int nbParts1D = (int) std::pow(float(nbParticles), 1.0f / 3.0f);
    for (uint i=0; i<nbParts1D; i++) {
        for (uint j=0; j<nbParts1D; j++) {
             for (uint k=0; k<nbParts1D; k++) {

                 const rp3d::Vector3 particlePosition(position.x - 0.5f * dimension.x + epsilon + i * 0.5f * (dimension.x / nbParts1D),
                                                position.y - 0.5f * dimension.y + epsilon + j * 0.5f * (dimension.y / nbParts1D),
                                                position.z - 0.5f * dimension.z + epsilon + k * 0.5f * (dimension.z / nbParts1D));

                 // Create a new particle
                 mFluid->createParticle(particlePosition, rp3d::Vector3(0.0, 0.0, 0.0));

             }
        }
    }
}

// Destructor
Scene::~Scene() {

    // Stop the physics simulation
    stopSimulation();

    // Destroy the shader
    mPhongShader.destroy();

    // Destroy the fluid from the world
    mDynamicsWorld->destroyParticleFluid(mFluid);

    // Destroy the rigid body of the floor
    mDynamicsWorld->destroyRigidBody(mFloor->getRigidBody());

    // Destroy the floor
    delete mFloor;

    // Destroy the dynamics world
    delete mDynamicsWorld;
}

// Take a step for the simulation
void Scene::simulate() {

    // If the physics simulation is running
    if (mIsRunning) {

        // Take a simulation step
        mDynamicsWorld->update();

        mFloor->updateTransform();
    }
}

// Render the scene
void Scene::render() {

    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_CULL_FACE);

    // Get the world-space to camera-space matrix
    const Camera& camera = mViewer->getCamera();
    const openglframework::Matrix4 worldToCameraMatrix = camera.getTransformMatrix().getInverse();

    // Bind the shader
    mPhongShader.bind();

    // Set the variables of the shader
    mPhongShader.setMatrix4x4Uniform("projectionMatrix", camera.getProjectionMatrix());
    mPhongShader.setVector3Uniform("light0PosCameraSpace",worldToCameraMatrix * mLight0.getOrigin());
    mPhongShader.setVector3Uniform("lightAmbientColor", Vector3(0.3f, 0.3f, 0.3f));
    const Color& diffCol = mLight0.getDiffuseColor();
    const Color& specCol = mLight0.getSpecularColor();
    mPhongShader.setVector3Uniform("light0DiffuseColor", Vector3(diffCol.r, diffCol.g, diffCol.b));
    mPhongShader.setVector3Uniform("light0SpecularColor", Vector3(specCol.r, specCol.g, specCol.b));
    mPhongShader.setFloatUniform("shininess", 60.0f);

    // Render the fluid
    renderFluid(mPhongShader, worldToCameraMatrix);

    // Render the floor
    // TODO : Uncomment this
    //mFloor->render(mPhongShader, worldToCameraMatrix);

    // Unbind the shader
    mPhongShader.unbind();
}

// Render the fluid
void Scene::renderFluid(openglframework::Shader& shader,
                        const openglframework::Matrix4& worldToCameraMatrix) {

    // For each particle of the fluid
    for (uint i=0; i<mFluid->getNbParticles(); i++) {

        const rp3d::FluidParticle& particle = mFluid->getParticle(i);

        // Translate it back to the origin
        Matrix4 matrix = mSphere->getTransformMatrix();
        matrix.m[0][3] = 0.0; matrix.m[1][3] = 0.0; matrix.m[2][3] = 0.0;
        mSphere->setTransformMatrix(matrix);

        // Translate the particle to its position
        mSphere->translateWorld(openglframework::Vector3(particle.position.x,
                                                         particle.position.y,
                                                         particle.position.z));
        mSphere->render(shader, worldToCameraMatrix);
    }
}
