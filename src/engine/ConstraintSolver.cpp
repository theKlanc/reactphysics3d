/********************************************************************************
* ReactPhysics3D physics library, http://code.google.com/p/reactphysics3d/      *
* Copyright (c) 2010-2012 Daniel Chappuis                                       *
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
#include "ConstraintSolver.h"
#include "DynamicsWorld.h"
#include "../body/RigidBody.h"

using namespace reactphysics3d;
using namespace std;

// Constants initialization
const decimal ConstraintSolver::BETA = 0.2;
const decimal ConstraintSolver::BETA_SPLIT_IMPULSE = 0.2;
const decimal ConstraintSolver::SLOP = 0.01;

// Constructor
ConstraintSolver::ConstraintSolver(DynamicsWorld* world)
    :world(world), nbConstraints(0), mNbIterations(10), mContactConstraints(0),
      mLinearVelocities(0), mAngularVelocities(0), mIsWarmStartingActive(true),
      mIsSplitImpulseActive(true), mIsSolveFrictionAtContactManifoldCenterActive(true) {

}

// Destructor
ConstraintSolver::~ConstraintSolver() {

}

// Initialize the constraint solver
void ConstraintSolver::initialize() {

    nbConstraints = 0;

    // TODO : Use better allocation here
    mContactConstraints = new ContactConstraint[world->getNbContactManifolds()];

    mNbContactConstraints = 0;

    // For each contact manifold of the world
    vector<ContactManifold>::iterator it;
    for (it = world->getContactManifoldsBeginIterator(); it != world->getContactManifoldsEndIterator(); ++it) {
        ContactManifold contactManifold = *it;

        ContactConstraint& constraint = mContactConstraints[mNbContactConstraints];

        assert(contactManifold.nbContacts > 0);

        RigidBody* body1 = contactManifold.contacts[0]->getBody1();
        RigidBody* body2 = contactManifold.contacts[0]->getBody2();

        // Fill in the body number maping
        mMapBodyToIndex.insert(make_pair(body1, mMapBodyToIndex.size()));
        mMapBodyToIndex.insert(make_pair(body2, mMapBodyToIndex.size()));

        // Add the two bodies of the constraint in the constraintBodies list
        mConstraintBodies.insert(body1);
        mConstraintBodies.insert(body2);

        Vector3 x1 = body1->getTransform().getPosition();
        Vector3 x2 = body2->getTransform().getPosition();

        constraint.indexBody1 = mMapBodyToIndex[body1];
        constraint.indexBody2 = mMapBodyToIndex[body2];
        constraint.inverseInertiaTensorBody1 = body1->getInertiaTensorInverseWorld();
        constraint.inverseInertiaTensorBody2 = body2->getInertiaTensorInverseWorld();
        constraint.isBody1Moving = body1->getIsMotionEnabled();
        constraint.isBody2Moving = body2->getIsMotionEnabled();
        constraint.massInverseBody1 = body1->getMassInverse();
        constraint.massInverseBody2 = body2->getMassInverse();
        constraint.nbContacts = contactManifold.nbContacts;
        constraint.restitutionFactor = computeMixRestitutionFactor(body1, body2);
        constraint.contactManifold = &contactManifold;

        // If we solve the friction constraints at the center of the contact manifold
        if (mIsSolveFrictionAtContactManifoldCenterActive) {
            constraint.frictionPointBody1 = Vector3(0.0, 0.0, 0.0);
            constraint.frictionPointBody2 = Vector3(0.0, 0.0, 0.0);
        }

        // For each contact point of the contact manifold
        for (uint c=0; c<contactManifold.nbContacts; c++) {

            ContactPointConstraint& contactPointConstraint = constraint.contacts[c];

            // Get a contact point
            Contact* contact = contactManifold.contacts[c];

            Vector3 p1 = contact->getWorldPointOnBody1();
            Vector3 p2 = contact->getWorldPointOnBody2();

            contactPointConstraint.contact = contact;
            contactPointConstraint.normal = contact->getNormal();
            contactPointConstraint.r1 = p1 - x1;
            contactPointConstraint.r2 = p2 - x2;
            contactPointConstraint.penetrationDepth = contact->getPenetrationDepth();
            contactPointConstraint.isRestingContact = contact->getIsRestingContact();
            contact->setIsRestingContact(true);
            contactPointConstraint.oldFrictionVector1 = contact->getFrictionVector1();
            contactPointConstraint.oldFrictionVector2 = contact->getFrictionVector2();
            contactPointConstraint.penetrationImpulse = 0.0;
            contactPointConstraint.friction1Impulse = 0.0;
            contactPointConstraint.friction2Impulse = 0.0;

            // If we solve the friction constraints at the center of the contact manifold
            if (mIsSolveFrictionAtContactManifoldCenterActive) {
                constraint.frictionPointBody1 += p1;
                constraint.frictionPointBody2 += p2;
            }
        }

        // If we solve the friction constraints at the center of the contact manifold
        if (mIsSolveFrictionAtContactManifoldCenterActive) {
            constraint.frictionPointBody1 /= constraint.nbContacts;
            constraint.frictionPointBody2 /= constraint.nbContacts;
            constraint.r1Friction = constraint.frictionPointBody1 - x1;
            constraint.r2Friction = constraint.frictionPointBody2 - x2;
            constraint.oldFrictionVector1 = contactManifold.frictionVector1;
            constraint.oldFrictionVector2 = contactManifold.frictionVector2;

            if (mIsWarmStartingActive) {
                constraint.friction1Impulse = contactManifold.friction1Impulse;
                constraint.friction2Impulse = contactManifold.friction2Impulse;
                constraint.frictionTwistImpulse = contactManifold.frictionTwistImpulse;
            }
            else {
                constraint.friction1Impulse = 0.0;
                constraint.friction2Impulse = 0.0;
                constraint.frictionTwistImpulse = 0.0;
            }
        }

        mNbContactConstraints++;
    }

    // Compute the number of bodies that are part of some active constraint
    nbBodies = mConstraintBodies.size();

    mLinearVelocities = new Vector3[nbBodies];
    mAngularVelocities = new Vector3[nbBodies];
    mSplitLinearVelocities = new Vector3[nbBodies];
    mSplitAngularVelocities = new Vector3[nbBodies];

    assert(mMapBodyToIndex.size() == nbBodies);
}

// Initialize the constrained bodies
void ConstraintSolver::initializeBodies() {

    // For each current body that is implied in some constraint
    RigidBody* rigidBody;
    for (set<RigidBody*>::iterator it = mConstraintBodies.begin(); it != mConstraintBodies.end(); ++it) {
        rigidBody = *it;
        uint bodyNumber = mMapBodyToIndex[rigidBody];

        // TODO : Use polymorphism and remove this downcasting
        assert(rigidBody);

        mLinearVelocities[bodyNumber] = rigidBody->getLinearVelocity() + mTimeStep * rigidBody->getMassInverse() * rigidBody->getExternalForce();
        mAngularVelocities[bodyNumber] = rigidBody->getAngularVelocity() + mTimeStep * rigidBody->getInertiaTensorInverseWorld() * rigidBody->getExternalTorque();
        mSplitLinearVelocities[bodyNumber] = Vector3(0, 0, 0);
        mSplitAngularVelocities[bodyNumber] = Vector3(0, 0, 0);
    }
}

// Initialize the contact constraints before solving the system
void ConstraintSolver::initializeContactConstraints() {
    
    // For each contact constraint
    for (uint c=0; c<mNbContactConstraints; c++) {

        ContactConstraint& constraint = mContactConstraints[c];

        Matrix3x3& I1 = constraint.inverseInertiaTensorBody1;
        Matrix3x3& I2 = constraint.inverseInertiaTensorBody2;

        // If we solve the friction constraints at the center of the contact manifold
        if (mIsSolveFrictionAtContactManifoldCenterActive) {
            constraint.normal = Vector3(0.0, 0.0, 0.0);
        }

        const Vector3& v1 = mLinearVelocities[constraint.indexBody1];
        const Vector3& w1 = mAngularVelocities[constraint.indexBody1];
        const Vector3& v2 = mLinearVelocities[constraint.indexBody2];
        const Vector3& w2 = mAngularVelocities[constraint.indexBody2];

        // For each contact point constraint
        for (uint i=0; i<constraint.nbContacts; i++) {

            ContactPointConstraint& contact = constraint.contacts[i];
            Contact* realContact = contact.contact;

            Vector3 deltaV = v2 + w2.cross(contact.r2) - v1 - w1.cross(contact.r1);

            contact.r1CrossN = contact.r1.cross(contact.normal);
            contact.r2CrossN = contact.r2.cross(contact.normal);

            decimal massPenetration = 0.0;
            if (constraint.isBody1Moving) {
                massPenetration += constraint.massInverseBody1 +
                        ((I1 * contact.r1CrossN).cross(contact.r1)).dot(contact.normal);
            }
            if (constraint.isBody2Moving) {
                massPenetration += constraint.massInverseBody2 +
                        ((I2 * contact.r2CrossN).cross(contact.r2)).dot(contact.normal);
            }
            massPenetration > 0.0 ? contact.inversePenetrationMass = 1.0 / massPenetration : 0.0;

            if (!mIsSolveFrictionAtContactManifoldCenterActive) {

                // Compute the friction vectors
                computeFrictionVectors(deltaV, contact);

                contact.r1CrossT1 = contact.r1.cross(contact.frictionVector1);
                contact.r1CrossT2 = contact.r1.cross(contact.frictionVector2);
                contact.r2CrossT1 = contact.r2.cross(contact.frictionVector1);
                contact.r2CrossT2 = contact.r2.cross(contact.frictionVector2);

                // Compute the inverse mass matrix K for the friction constraints at each contact point
                decimal friction1Mass = 0.0;
                decimal friction2Mass = 0.0;
                if (constraint.isBody1Moving) {
                    friction1Mass += constraint.massInverseBody1 + ((I1 * contact.r1CrossT1).cross(contact.r1)).dot(contact.frictionVector1);
                    friction2Mass += constraint.massInverseBody1 + ((I1 * contact.r1CrossT2).cross(contact.r1)).dot(contact.frictionVector2);
                }
                if (constraint.isBody2Moving) {
                    friction1Mass += constraint.massInverseBody2 + ((I2 * contact.r2CrossT1).cross(contact.r2)).dot(contact.frictionVector1);
                    friction2Mass += constraint.massInverseBody2 + ((I2 * contact.r2CrossT2).cross(contact.r2)).dot(contact.frictionVector2);
                }
                friction1Mass > 0.0 ? contact.inverseFriction1Mass = 1.0 / friction1Mass : 0.0;
                friction2Mass > 0.0 ? contact.inverseFriction2Mass = 1.0 / friction2Mass : 0.0;
            }

            // Compute the restitution velocity bias "b". We compute this here instead
            // of inside the solve() method because we need to use the velocity difference
            // at the beginning of the contact. Note that if it is a resting contact (normal velocity
            // under a given threshold), we don't add a restitution velocity bias
            contact.restitutionBias = 0.0;
            decimal deltaVDotN = deltaV.dot(contact.normal);
            // TODO : Use a constant here
            if (deltaVDotN < 1.0f) {
                contact.restitutionBias = constraint.restitutionFactor * deltaVDotN;
            }

            // Get the cached lambda values of the constraint
            if (mIsWarmStartingActive) {
                contact.penetrationImpulse = realContact->getCachedLambda(0);
                contact.friction1Impulse = realContact->getCachedLambda(1);
                contact.friction2Impulse = realContact->getCachedLambda(2);
            }

            // Initialize the split impulses to zero
            contact.penetrationSplitImpulse = 0.0;

            // If we solve the friction constraints at the center of the contact manifold
            if (mIsSolveFrictionAtContactManifoldCenterActive) {
                constraint.normal += contact.normal;
            }
        }

        // If we solve the friction constraints at the center of the contact manifold
        if (mIsSolveFrictionAtContactManifoldCenterActive) {

            constraint.normal.normalize();

            Vector3 deltaVFrictionPoint = v2 + w2.cross(constraint.r2Friction) -
                                          v1 - w1.cross(constraint.r1Friction);

            // Compute the friction vectors
            computeFrictionVectors(deltaVFrictionPoint, constraint);

            // Compute the inverse mass matrix K for the friction constraints at the center of
            // the contact manifold
            constraint.r1CrossT1 = constraint.r1Friction.cross(constraint.frictionVector1);
            constraint.r1CrossT2 = constraint.r1Friction.cross(constraint.frictionVector2);
            constraint.r2CrossT1 = constraint.r2Friction.cross(constraint.frictionVector1);
            constraint.r2CrossT2 = constraint.r2Friction.cross(constraint.frictionVector2);
            decimal friction1Mass = 0.0;
            decimal friction2Mass = 0.0;
            if (constraint.isBody1Moving) {
                friction1Mass += constraint.massInverseBody1 + ((I1 * constraint.r1CrossT1).cross(constraint.r1Friction)).dot(constraint.frictionVector1);
                friction2Mass += constraint.massInverseBody1 + ((I1 * constraint.r1CrossT2).cross(constraint.r1Friction)).dot(constraint.frictionVector2);
            }
            if (constraint.isBody2Moving) {
                friction1Mass += constraint.massInverseBody2 + ((I2 * constraint.r2CrossT1).cross(constraint.r2Friction)).dot(constraint.frictionVector1);
                friction2Mass += constraint.massInverseBody2 + ((I2 * constraint.r2CrossT2).cross(constraint.r2Friction)).dot(constraint.frictionVector2);
            }
            friction1Mass > 0.0 ? constraint.inverseFriction1Mass = 1.0 / friction1Mass : 0.0;
            friction2Mass > 0.0 ? constraint.inverseFriction2Mass = 1.0 / friction2Mass : 0.0;
        }
    }
}

// Warm start the solver
// For each constraint, we apply the previous impulse (from the previous step)
// at the beginning. With this technique, we will converge faster towards
// the solution of the linear system
void ConstraintSolver::warmStart() {

    // For each constraint
    for (uint c=0; c<mNbContactConstraints; c++) {

        ContactConstraint& constraint = mContactConstraints[c];

        bool atLeastOneRestingContactPoint = false;

        for (uint i=0; i<constraint.nbContacts; i++) {

            ContactPointConstraint& contact = constraint.contacts[i];

            // If it is not a new contact (this contact was already existing at last time step)
            if (contact.isRestingContact) {

                atLeastOneRestingContactPoint = true;

                // --------- Penetration --------- //

                // Compute the impulse P=J^T * lambda
                const Vector3 linearImpulseBody1 = -contact.normal * contact.penetrationImpulse;
                const Vector3 angularImpulseBody1 = -contact.r1CrossN * contact.penetrationImpulse;
                const Vector3 linearImpulseBody2 = contact.normal * contact.penetrationImpulse;
                const Vector3 angularImpulseBody2 = contact.r2CrossN * contact.penetrationImpulse;
                const Impulse impulsePenetration(linearImpulseBody1, angularImpulseBody1,
                                                 linearImpulseBody2, angularImpulseBody2);

                // Apply the impulse to the bodies of the constraint
                applyImpulse(impulsePenetration, constraint);

                // If we do not solve the friction constraints at the center of the contact manifold
                if (!mIsSolveFrictionAtContactManifoldCenterActive) {

                    // Project the old friction impulses (with old friction vectors) into the new friction
                    // vectors to get the new friction impulses
                    Vector3 oldFrictionImpulse = contact.friction1Impulse * contact.oldFrictionVector1 +
                            contact.friction2Impulse * contact.oldFrictionVector2;
                    contact.friction1Impulse = oldFrictionImpulse.dot(contact.frictionVector1);
                    contact.friction2Impulse = oldFrictionImpulse.dot(contact.frictionVector2);

                    // --------- Friction 1 --------- //

                    // Compute the impulse P=J^T * lambda
                    Vector3 linearImpulseBody1Friction1 = -contact.frictionVector1 * contact.friction1Impulse;
                    Vector3 angularImpulseBody1Friction1 = -contact.r1CrossT1 * contact.friction1Impulse;
                    Vector3 linearImpulseBody2Friction1 = contact.frictionVector1 * contact.friction1Impulse;
                    Vector3 angularImpulseBody2Friction1 = contact.r2CrossT1 * contact.friction1Impulse;
                    Impulse impulseFriction1(linearImpulseBody1Friction1, angularImpulseBody1Friction1,
                                             linearImpulseBody2Friction1, angularImpulseBody2Friction1);

                    // Apply the impulses to the bodies of the constraint
                    applyImpulse(impulseFriction1, constraint);

                    // --------- Friction 2 --------- //

                    // Compute the impulse P=J^T * lambda
                    Vector3 linearImpulseBody1Friction2 = -contact.frictionVector2 * contact.friction2Impulse;
                    Vector3 angularImpulseBody1Friction2 = -contact.r1CrossT2 * contact.friction2Impulse;
                    Vector3 linearImpulseBody2Friction2 = contact.frictionVector2 * contact.friction2Impulse;
                    Vector3 angularImpulseBody2Friction2 = contact.r2CrossT2 * contact.friction2Impulse;
                    Impulse impulseFriction2(linearImpulseBody1Friction2, angularImpulseBody1Friction2,
                                             linearImpulseBody2Friction2, angularImpulseBody2Friction2);

                    // Apply the impulses to the bodies of the constraint
                    applyImpulse(impulseFriction2, constraint);
                }
            }
            else {  // If it is a new contact point

                // Initialize the accumulated impulses to zero
                contact.penetrationImpulse = 0.0;
                contact.friction1Impulse = 0.0;
                contact.friction2Impulse = 0.0;
            }
        }

        // If we solve the friction constraints at the center of the contact manifold and there is
        // at least one resting contact point in the contact manifold
        if (mIsSolveFrictionAtContactManifoldCenterActive && atLeastOneRestingContactPoint) {

            // Project the old friction impulses (with old friction vectors) into the new friction
            // vectors to get the new friction impulses
            Vector3 oldFrictionImpulse = constraint.friction1Impulse * constraint.oldFrictionVector1 +
                    constraint.friction2Impulse * constraint.oldFrictionVector2;
            constraint.friction1Impulse = oldFrictionImpulse.dot(constraint.frictionVector1);
            constraint.friction2Impulse = oldFrictionImpulse.dot(constraint.frictionVector2);

            // ------ First friction constraint at the center of the contact manifol ------ //

            // Compute the impulse P=J^T * lambda
            Vector3 linearImpulseBody1 = -constraint.frictionVector1 * constraint.friction1Impulse;
            Vector3 angularImpulseBody1 = -constraint.r1CrossT1 * constraint.friction1Impulse;
            Vector3 linearImpulseBody2 = constraint.frictionVector1 * constraint.friction1Impulse;
            Vector3 angularImpulseBody2 = constraint.r2CrossT1 * constraint.friction1Impulse;
            const Impulse impulseFriction1(linearImpulseBody1, angularImpulseBody1,
                                           linearImpulseBody2, angularImpulseBody2);

            // Apply the impulses to the bodies of the constraint
            applyImpulse(impulseFriction1, constraint);

            // ------ Second friction constraint at the center of the contact manifol ----- //

            // Compute the impulse P=J^T * lambda
            linearImpulseBody1 = -constraint.frictionVector2 * constraint.friction2Impulse;
            angularImpulseBody1 = -constraint.r1CrossT2 * constraint.friction2Impulse;
            linearImpulseBody2 = constraint.frictionVector2 * constraint.friction2Impulse;
            angularImpulseBody2 = constraint.r2CrossT2 * constraint.friction2Impulse;
            const Impulse impulseFriction2(linearImpulseBody1, angularImpulseBody1,
                                           linearImpulseBody2, angularImpulseBody2);

            // Apply the impulses to the bodies of the constraint
            applyImpulse(impulseFriction2, constraint);

            // ------ Twist friction constraint at the center of the contact manifol ------ //


            // Compute the impulse P=J^T * lambda
            linearImpulseBody1 = Vector3(0.0, 0.0, 0.0);
            angularImpulseBody1 = -constraint.normal * constraint.frictionTwistImpulse;
            linearImpulseBody2 = Vector3(0.0, 0.0, 0.0);
            angularImpulseBody2 = constraint.normal * constraint.frictionTwistImpulse;
            const Impulse impulseTwistFriction(linearImpulseBody1, angularImpulseBody1,
                                               linearImpulseBody2, angularImpulseBody2);

            // Apply the impulses to the bodies of the constraint
            applyImpulse(impulseTwistFriction, constraint);
        }
        else {  // If it is a new contact manifold

            // Initialize the accumulated impulses to zero
            constraint.friction1Impulse = 0.0;
            constraint.friction2Impulse = 0.0;
            constraint.frictionTwistImpulse = 0.0;
        }
    }
}

// Solve the contact constraints by applying sequential impulses
void ConstraintSolver::solveContactConstraints() {

    decimal deltaLambda;
    decimal lambdaTemp;
    uint iter;

    // For each iteration
    for(iter=0; iter<mNbIterations; iter++) {

        // For each constraint
        for (uint c=0; c<mNbContactConstraints; c++) {

            ContactConstraint& constraint = mContactConstraints[c];

            decimal sumPenetrationImpulse = 0.0;

            const Vector3& v1 = mLinearVelocities[constraint.indexBody1];
            const Vector3& w1 = mAngularVelocities[constraint.indexBody1];
            const Vector3& v2 = mLinearVelocities[constraint.indexBody2];
            const Vector3& w2 = mAngularVelocities[constraint.indexBody2];

            for (uint i=0; i<constraint.nbContacts; i++) {

                ContactPointConstraint& contact = constraint.contacts[i];

                // --------- Penetration --------- //

                // Compute J*v
                Vector3 deltaV = v2 + w2.cross(contact.r2) - v1 - w1.cross(contact.r1);
                decimal deltaVDotN = deltaV.dot(contact.normal);
                decimal Jv = deltaVDotN;

                // Compute the bias "b" of the constraint
                decimal beta = mIsSplitImpulseActive ? BETA_SPLIT_IMPULSE : BETA;
                decimal biasPenetrationDepth = 0.0;
                if (contact.penetrationDepth > SLOP) biasPenetrationDepth = -(beta/mTimeStep) *
                        max(0.0f, float(contact.penetrationDepth - SLOP));
                decimal b = biasPenetrationDepth + contact.restitutionBias;

                // Compute the Lagrange multiplier
                if (mIsSplitImpulseActive) {
                    deltaLambda = - (Jv + contact.restitutionBias) * contact.inversePenetrationMass;
                }
                else {
                    deltaLambda = - (Jv + b) * contact.inversePenetrationMass;
                }
                lambdaTemp = contact.penetrationImpulse;
                contact.penetrationImpulse = std::max(contact.penetrationImpulse + deltaLambda, 0.0f);
                deltaLambda = contact.penetrationImpulse - lambdaTemp;

                // Compute the impulse P=J^T * lambda
                Vector3 linearImpulseBody1 = -contact.normal * deltaLambda;
                Vector3 angularImpulseBody1 = -contact.r1CrossN * deltaLambda;
                Vector3 linearImpulseBody2 = contact.normal * deltaLambda;
                Vector3 angularImpulseBody2 = contact.r2CrossN * deltaLambda;
                const Impulse impulsePenetration(linearImpulseBody1, angularImpulseBody1,
                                                 linearImpulseBody2, angularImpulseBody2);

                // Apply the impulse to the bodies of the constraint
                applyImpulse(impulsePenetration, constraint);

                sumPenetrationImpulse += contact.penetrationImpulse;

                // If the split impulse position correction is active
                if (mIsSplitImpulseActive) {

                    // Split impulse (position correction)
                    const Vector3& v1Split = mSplitLinearVelocities[constraint.indexBody1];
                    const Vector3& w1Split = mSplitAngularVelocities[constraint.indexBody1];
                    const Vector3& v2Split = mSplitLinearVelocities[constraint.indexBody2];
                    const Vector3& w2Split = mSplitAngularVelocities[constraint.indexBody2];
                    Vector3 deltaVSplit = v2Split + w2Split.cross(contact.r2) - v1Split - w1Split.cross(contact.r1);
                    decimal JvSplit = deltaVSplit.dot(contact.normal);
                    decimal deltaLambdaSplit = - (JvSplit + biasPenetrationDepth) * contact.inversePenetrationMass;
                    decimal lambdaTempSplit = contact.penetrationSplitImpulse;
                    contact.penetrationSplitImpulse = std::max(contact.penetrationSplitImpulse + deltaLambdaSplit, 0.0f);
                    deltaLambda = contact.penetrationSplitImpulse - lambdaTempSplit;

                    // Compute the impulse P=J^T * lambda
                    linearImpulseBody1 = -contact.normal * deltaLambdaSplit;
                    angularImpulseBody1 = -contact.r1CrossN * deltaLambdaSplit;
                    linearImpulseBody2 = contact.normal * deltaLambdaSplit;
                    angularImpulseBody2 = contact.r2CrossN * deltaLambdaSplit;
                    const Impulse splitImpulsePenetration(linearImpulseBody1, angularImpulseBody1,
                                                          linearImpulseBody2, angularImpulseBody2);

                    applySplitImpulse(splitImpulsePenetration, constraint);
                }

                // If we do not solve the friction constraints at the center of the contact manifold
                if (!mIsSolveFrictionAtContactManifoldCenterActive) {

                    // --------- Friction 1 --------- //

                    // Compute J*v
                    deltaV = v2 + w2.cross(contact.r2) - v1 - w1.cross(contact.r1);
                    Jv = deltaV.dot(contact.frictionVector1);

                    deltaLambda = -Jv;
                    deltaLambda *= contact.inverseFriction1Mass;
                    decimal frictionLimit = 0.3 * contact.penetrationImpulse;   // TODO : Use constant here
                    lambdaTemp = contact.friction1Impulse;
                    contact.friction1Impulse = std::max(-frictionLimit, std::min(contact.friction1Impulse + deltaLambda, frictionLimit));
                    deltaLambda = contact.friction1Impulse - lambdaTemp;

                    // Compute the impulse P=J^T * lambda
                    linearImpulseBody1 = -contact.frictionVector1 * deltaLambda;
                    angularImpulseBody1 = -contact.r1CrossT1 * deltaLambda;
                    linearImpulseBody2 = contact.frictionVector1 * deltaLambda;
                    angularImpulseBody2 = contact.r2CrossT1 * deltaLambda;
                    const Impulse impulseFriction1(linearImpulseBody1, angularImpulseBody1,
                                                   linearImpulseBody2, angularImpulseBody2);

                    // Apply the impulses to the bodies of the constraint
                    applyImpulse(impulseFriction1, constraint);

                    // --------- Friction 2 --------- //

                    // Compute J*v
                    deltaV = v2 + w2.cross(contact.r2) - v1 - w1.cross(contact.r1);
                    Jv = deltaV.dot(contact.frictionVector2);

                    deltaLambda = -Jv;
                    deltaLambda *= contact.inverseFriction2Mass;
                    frictionLimit = 0.3 * contact.penetrationImpulse;   // TODO : Use constant here
                    lambdaTemp = contact.friction2Impulse;
                    contact.friction2Impulse = std::max(-frictionLimit, std::min(contact.friction2Impulse + deltaLambda, frictionLimit));
                    deltaLambda = contact.friction2Impulse - lambdaTemp;

                    // Compute the impulse P=J^T * lambda
                    linearImpulseBody1 = -contact.frictionVector2 * deltaLambda;
                    angularImpulseBody1 = -contact.r1CrossT2 * deltaLambda;
                    linearImpulseBody2 = contact.frictionVector2 * deltaLambda;
                    angularImpulseBody2 = contact.r2CrossT2 * deltaLambda;
                    const Impulse impulseFriction2(linearImpulseBody1, angularImpulseBody1,
                                                   linearImpulseBody2, angularImpulseBody2);

                    // Apply the impulses to the bodies of the constraint
                    applyImpulse(impulseFriction2, constraint);
                }
            }

            // If we solve the friction constraints at the center of the contact manifold
            if (mIsSolveFrictionAtContactManifoldCenterActive) {

                // ------ First friction constraint at the center of the contact manifol ------ //

                // Compute J*v
                Vector3 deltaV = v2 + w2.cross(constraint.r2Friction) - v1 - w1.cross(constraint.r1Friction);
                decimal Jv = deltaV.dot(constraint.frictionVector1);

                decimal deltaLambda = -Jv * constraint.inverseFriction1Mass;
                decimal frictionLimit = 0.3 * sumPenetrationImpulse;   // TODO : Use constant here
                lambdaTemp = constraint.friction1Impulse;
                constraint.friction1Impulse = std::max(-frictionLimit, std::min(constraint.friction1Impulse + deltaLambda, frictionLimit));
                deltaLambda = constraint.friction1Impulse - lambdaTemp;

                // Compute the impulse P=J^T * lambda
                Vector3 linearImpulseBody1 = -constraint.frictionVector1 * deltaLambda;
                Vector3 angularImpulseBody1 = -constraint.r1CrossT1 * deltaLambda;
                Vector3 linearImpulseBody2 = constraint.frictionVector1 * deltaLambda;
                Vector3 angularImpulseBody2 = constraint.r2CrossT1 * deltaLambda;
                const Impulse impulseFriction1(linearImpulseBody1, angularImpulseBody1,
                                               linearImpulseBody2, angularImpulseBody2);

                // Apply the impulses to the bodies of the constraint
                applyImpulse(impulseFriction1, constraint);

                // ------ Second friction constraint at the center of the contact manifol ----- //

                // Compute J*v
                deltaV = v2 + w2.cross(constraint.r2Friction) - v1 - w1.cross(constraint.r1Friction);
                Jv = deltaV.dot(constraint.frictionVector2);

                deltaLambda = -Jv * constraint.inverseFriction2Mass;
                frictionLimit = 0.3 * sumPenetrationImpulse;   // TODO : Use constant here
                lambdaTemp = constraint.friction2Impulse;
                constraint.friction2Impulse = std::max(-frictionLimit, std::min(constraint.friction2Impulse + deltaLambda, frictionLimit));
                deltaLambda = constraint.friction2Impulse - lambdaTemp;

                // Compute the impulse P=J^T * lambda
                linearImpulseBody1 = -constraint.frictionVector2 * deltaLambda;
                angularImpulseBody1 = -constraint.r1CrossT2 * deltaLambda;
                linearImpulseBody2 = constraint.frictionVector2 * deltaLambda;
                angularImpulseBody2 = constraint.r2CrossT2 * deltaLambda;
                const Impulse impulseFriction2(linearImpulseBody1, angularImpulseBody1,
                                               linearImpulseBody2, angularImpulseBody2);

                // Apply the impulses to the bodies of the constraint
                applyImpulse(impulseFriction2, constraint);

                // ------ Twist friction constraint at the center of the contact manifol ------ //

                // TODO : Put this in the initialization method
                decimal K = constraint.normal.dot(constraint.inverseInertiaTensorBody1 * constraint.normal) +
                           constraint.normal.dot(constraint.inverseInertiaTensorBody2 * constraint.normal);


                // Compute J*v
                deltaV = w2 - w1;
                Jv = deltaV.dot(constraint.normal);

                // TODO : Compute the inverse mass matrix here for twist friction
                deltaLambda = -Jv * (1.0 / K);
                frictionLimit = 0.3 * sumPenetrationImpulse;   // TODO : Use constant here
                lambdaTemp = constraint.frictionTwistImpulse;
                constraint.frictionTwistImpulse = std::max(-frictionLimit, std::min(constraint.frictionTwistImpulse + deltaLambda, frictionLimit));
                deltaLambda = constraint.frictionTwistImpulse - lambdaTemp;

                // Compute the impulse P=J^T * lambda
                linearImpulseBody1 = Vector3(0.0, 0.0, 0.0);
                angularImpulseBody1 = -constraint.normal * deltaLambda;
                linearImpulseBody2 = Vector3(0.0, 0.0, 0.0);;
                angularImpulseBody2 = constraint.normal * deltaLambda;
                const Impulse impulseTwistFriction(linearImpulseBody1, angularImpulseBody1,
                                                   linearImpulseBody2, angularImpulseBody2);

                // Apply the impulses to the bodies of the constraint
                applyImpulse(impulseTwistFriction, constraint);
            }
        }
    }
}

// Solve the constraints
void ConstraintSolver::solve(decimal timeStep) {

    mTimeStep = timeStep;

    // Initialize the solver
    initialize();

    initializeBodies();

    // Fill-in all the matrices needed to solve the LCP problem
    initializeContactConstraints();

    // Warm start the solver
    if (mIsWarmStartingActive) {
        warmStart();
    }

    // Solve the contact constraints
    solveContactConstraints();

    // Cache the lambda values in order to use them in the next step
    storeImpulses();
}

// Store the computed impulses to use them to
// warm start the solver at the next iteration
void ConstraintSolver::storeImpulses() {

    // For each constraint
    for (uint c=0; c<mNbContactConstraints; c++) {

        ContactConstraint& constraint = mContactConstraints[c];

        for (uint i=0; i<constraint.nbContacts; i++) {

            ContactPointConstraint& contact = constraint.contacts[i];

            contact.contact->setCachedLambda(0, contact.penetrationImpulse);
            contact.contact->setCachedLambda(1, contact.friction1Impulse);
            contact.contact->setCachedLambda(2, contact.friction2Impulse);

            contact.contact->setFrictionVector1(contact.frictionVector1);
            contact.contact->setFrictionVector2(contact.frictionVector2);
        }

        constraint.contactManifold->friction1Impulse = constraint.friction1Impulse;
        constraint.contactManifold->friction2Impulse = constraint.friction2Impulse;
        constraint.contactManifold->frictionTwistImpulse = constraint.frictionTwistImpulse;
        constraint.contactManifold->frictionVector1 = constraint.frictionVector1;
        constraint.contactManifold->frictionVector2 = constraint.frictionVector2;
    }
}

// Apply an impulse to the two bodies of a constraint
void ConstraintSolver::applyImpulse(const Impulse& impulse, const ContactConstraint& constraint) {

    // Update the velocities of the bodies by applying the impulse P
    if (constraint.isBody1Moving) {
        mLinearVelocities[constraint.indexBody1] += constraint.massInverseBody1 *
                                                    impulse.linearImpulseBody1;
        mAngularVelocities[constraint.indexBody1] += constraint.inverseInertiaTensorBody1 *
                                                     impulse.angularImpulseBody1;
    }
    if (constraint.isBody2Moving) {
        mLinearVelocities[constraint.indexBody2] += constraint.massInverseBody2 *
                                                    impulse.linearImpulseBody2;
        mAngularVelocities[constraint.indexBody2] += constraint.inverseInertiaTensorBody2 *
                                                     impulse.angularImpulseBody2;
    }
}

// Apply an impulse to the two bodies of a constraint
void ConstraintSolver::applySplitImpulse(const Impulse& impulse, const ContactConstraint& constraint) {

    // Update the velocities of the bodies by applying the impulse P
    if (constraint.isBody1Moving) {
        mSplitLinearVelocities[constraint.indexBody1] += constraint.massInverseBody1 *
                                                    impulse.linearImpulseBody1;
        mSplitAngularVelocities[constraint.indexBody1] += constraint.inverseInertiaTensorBody1 *
                                                     impulse.angularImpulseBody1;
    }
    if (constraint.isBody2Moving) {
        mSplitLinearVelocities[constraint.indexBody2] += constraint.massInverseBody2 *
                                                    impulse.linearImpulseBody2;
        mSplitAngularVelocities[constraint.indexBody2] += constraint.inverseInertiaTensorBody2 *
                                                     impulse.angularImpulseBody2;
    }
}

// Compute the two unit orthogonal vectors "t1" and "t2" that span the tangential friction plane
// for a contact point constraint. The two vectors have to be such that : t1 x t2 = contactNormal.
void ConstraintSolver::computeFrictionVectors(const Vector3& deltaVelocity,
                                              ContactPointConstraint& contact) const {

    // Update the old friction vectors
    //contact.oldFrictionVector1 = contact.frictionVector1;
    //contact.oldFrictionVector2 = contact.frictionVector2;

    assert(contact.normal.length() > 0.0);

    // Compute the velocity difference vector in the tangential plane
    Vector3 normalVelocity = deltaVelocity.dot(contact.normal) * contact.normal;
    Vector3 tangentVelocity = deltaVelocity - normalVelocity;

    // If the velocty difference in the tangential plane is not zero
    decimal lengthTangenVelocity = tangentVelocity.length();
    if (lengthTangenVelocity > 0.0) {

        // Compute the first friction vector in the direction of the tangent
        // velocity difference
        contact.frictionVector1 = tangentVelocity / lengthTangenVelocity;
    }
    else {

        // Get any orthogonal vector to the normal as the first friction vector
        contact.frictionVector1 = contact.normal.getOneUnitOrthogonalVector();
    }

    // The second friction vector is computed by the cross product of the firs
    // friction vector and the contact normal
    contact.frictionVector2 = contact.normal.cross(contact.frictionVector1).getUnit();
}

// Compute the two unit orthogonal vectors "t1" and "t2" that span the tangential friction plane
// for a contact constraint. The two vectors have to be such that : t1 x t2 = contactNormal.
void ConstraintSolver::computeFrictionVectors(const Vector3& deltaVelocity,
                                              ContactConstraint& contact) const {

    assert(contact.normal.length() > 0.0);

    // Compute the velocity difference vector in the tangential plane
    Vector3 normalVelocity = deltaVelocity.dot(contact.normal) * contact.normal;
    Vector3 tangentVelocity = deltaVelocity - normalVelocity;

    // If the velocty difference in the tangential plane is not zero
    decimal lengthTangenVelocity = tangentVelocity.length();
    if (lengthTangenVelocity > 0.0) {

        // Compute the first friction vector in the direction of the tangent
        // velocity difference
        contact.frictionVector1 = tangentVelocity / lengthTangenVelocity;
    }
    else {

        // Get any orthogonal vector to the normal as the first friction vector
        contact.frictionVector1 = contact.normal.getOneUnitOrthogonalVector();
    }

    // The second friction vector is computed by the cross product of the firs
    // friction vector and the contact normal
    contact.frictionVector2 = contact.normal.cross(contact.frictionVector1).getUnit();
}
