/****************************************************************************
 * Copyright (C) 2009      Daniel Chappuis                                  *
 ****************************************************************************
 * This file is part of ReactPhysics3D.                                     *
 *                                                                          *
 * ReactPhysics3D is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU General Public License as published by     *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 *                                                                          *
 * ReactPhysics3D is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 * GNU General Public License for more details.                             *
 *                                                                          *
 * You should have received a copy of the GNU General Public License        *
 * along with ReactPhysics3D. If not, see <http://www.gnu.org/licenses/>.   *
 ***************************************************************************/

 #ifndef RIGIDBODY_H
 #define RIGIDBODY_H

 // Libraries
 #include "Body.h"
 #include "../mathematics/mathematics.h"

// Namespace reactphysics3d
namespace reactphysics3d {

/*  -------------------------------------------------------------------
    Class RigidBody :
        This class represents a rigid body of the physics
        engine. A rigid body is a non-deformable body that
        has a constant mass.
    -------------------------------------------------------------------
*/
class RigidBody : public Body {
    private :
        Matrix3x3 inertiaTensor;                    // Inertia tensor of the body
        BodyState currentBodyState;                 // Current body state
        BodyState previousBodyState;                // Previous body state
        bool isMotionEnabled;                       // True if the body can move
        bool isCollisionEnabled;                    // True if the body can collide with others bodies

    public :
        RigidBody(double mass);                     // Constructor
        RigidBody(const RigidBody& rigidBody);      // Copy-constructor
        virtual ~RigidBody();                       // Destructor

        Matrix3x3 getInertiaTensor() const;                             // Return the inertia tensor of the body
        void setInertiaTensor(const Matrix3x3& inertiaTensor);          // Set the inertia tensor of the body
        BodyState getCurrentState() const;                              // Return the current state of the body
        BodyState getPreviousState() const;                             // Return the previous state of the body
};


}

 #endif
