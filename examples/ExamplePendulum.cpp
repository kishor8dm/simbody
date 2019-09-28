/* -------------------------------------------------------------------------- *
 *                       Simbody(tm) Example: Pendulum                        *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Ajay Seth                                                         *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */


// Define two identical double pendulums, one modeled the easy way using
// two pin mobilizers, the other modeled with free mobilizers plus a ball
// constraint, plus two "constant angle" constraints to get rid of the extra 
// rotational degrees of freedom.
// We're going to show that the resulting reaction forces are identical.


#include "Simbody.h"

using namespace SimTK;
using std::cout; using std::endl;
namespace {
bool enableVis = false;
}

int main() {
    try {
        Quaternion q0(-0.1504133023, -0.7202261496,  0.6643972918,  0.1312492688);
        Quaternion q1(-0.8802997333,  -0.4219333635,  -0.1958447339, -0.09321725385);
        Quaternion q2(0.1707965794, -0.1900937466, -0.1176534344, -0.9596095901);
        Quaternion q3(0.3201724846,  0.7169894904,  0.3068382102, -0.5378345131);
        Quaternion q4(-0.2366251647, 0.3815343678,  0.845472455,   0.28916502);
        Quaternion q5(-0.8099961877, 0.09318492562,  0.5701901317, -0.1005286202);
        Quaternion q6(0.09452920407, -0.6611310244, -0.5389370918, -0.5133388638);
        Quaternion q7(-0.06702885032, -0.7941805124,   -0.17857261,  0.5769716501);
        Quaternion eye(1, 0, 0, 0);
        Quaternion halfPi(0.7071067691, 0.7071067691, 0, 0);
        // Create the system, with subsystems for the bodies and some forces.
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        GeneralForceSubsystem forces(system);
        Force::UniformGravity gravity(forces, matter, Vec3(0, Real(0.0), 0));
        Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
        if (enableVis) {
            pendulumBody.addDecoration(Transform(),DecorativeSphere(Real(0.1)).setColor(Red));
        }

        auto local0 = Vec3(0, -1, 0);
        auto local1 = Vec3(0, 1, 0);

        MobilizedBody::Pin body0(matter.updGround(), {Rotation(q0), local0}, pendulumBody, {Rotation(q1),local1});
        MobilizedBody::Pin body1(body0, Transform(Rotation(q2), local0), pendulumBody,Transform(Rotation(q3), local1));
        MobilizedBody::Pin body2(body1, Transform(Rotation(q4), local0), pendulumBody,Transform(Rotation(q5), local1));
        MobilizedBody::Pin body3(body1, Transform(Rotation(q6), local0), pendulumBody,Transform(Rotation(q7), local1));
        MobilizedBody::Pin body4(body3, Transform(Rotation(q0), local0), pendulumBody,Transform(Rotation(q1), local1));

        std::vector<MobilizedBody::Pin*> bodyList = {&body0, &body1, &body2, &body3, &body4};

        Visualizer* viz = nullptr;
        if (enableVis) {
            viz = new Visualizer(system);
            system.addEventReporter(new Visualizer::Reporter(*viz, Real(1./30)));
        }

        // Initialize the system and state.

        system.realizeTopology();
        State state = system.getDefaultState();
        for (auto body : bodyList) {
            body->setOneQ(state, 0, Pi/4);
            body->setOneU(state, 0, 1);
        }

        system.realize(state);
        for (auto body : bodyList) {
            std:: cout << body->getBodyRotation(state).convertRotationToQuaternion() << "\n";
        }
        std::cout << "--------------\n";
        for (auto body : bodyList) {
            std:: cout << body->getBodyOriginLocation(state) << "\n";
        }
        std::cout << "--------------\n";
        for (auto body : bodyList) {
            std:: cout << body->getBodyAngularVelocity(state) << "\n";
        }
        std::cout << "--------------\n";
        for (auto body : bodyList) {
            std:: cout << body->getBodyOriginVelocity(state) << "\n";
        }
        if (enableVis) {
            viz->report(state);
            RungeKuttaMersonIntegrator integ(system);
            TimeStepper ts(system, integ);
            ts.initialize(state);
            ts.stepTo(10000000.0);
        }
    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cout << "UNKNOWN EXCEPTION\n";
        return 1;
    }
    return 0;
}
