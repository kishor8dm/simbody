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

int main() {
  try {
      Quaternion q0(-0.1504133023, -0.7202261496,  0.6643972918,  0.1312492688);
      Quaternion q1(-0.8802997333,  -0.4219333635,  -0.1958447339, -0.09321725385);
      Quaternion q2(0.1707965794, -0.1900937466, -0.1176534344, -0.9596095901);
      Quaternion q3(0.3201724846,  0.7169894904,  0.3068382102, -0.5378345131);
      Quaternion eye(1, 0, 0, 0);
      Quaternion halfPi(0.7071067691, 0.7071067691, 0, 0);
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Add gravity as a force element.
    Rotation x45(Pi/4, XAxis);
    Rotation y45(Pi/4, YAxis);
    Rotation z45(Pi/4, ZAxis);
    Force::UniformGravity gravity(forces, matter, Vec3(0, Real(0.0), 0));
    // Create the body and some artwork for it.
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), 
                               DecorativeSphere(Real(0.1)).setColor(Red));

      MobilizedBody::Pin pendulum0(matter.updGround(),{Rotation(eye), Vec3(0,-1,0)},
                                   pendulumBody,{Rotation(eye),Vec3(0, 1, 0)});
      MobilizedBody::Pin pendulum1(pendulum0,
                                   Transform(Rotation(halfPi), Vec3(0,-1,0)),
                                   pendulumBody,
                                   Transform(Rotation(eye), Vec3(0, 1, 0)));


    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, Real(1./30)));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum0.setOneQ(state, 0, 0);
    pendulum1.setOneQ(state, 0, 0);

    system.realize(state);
      std::cout << pendulum0.getBodyRotation(state).convertRotationToQuaternion() << "\n";
      std::cout << pendulum1.getBodyRotation(state).convertRotationToQuaternion() << "\n";
      std::cout << pendulum0.getBodyOriginLocation(state) << "\n";
      std::cout << pendulum1.getBodyOriginLocation(state) << "\n";
    viz.report(state);
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10000000.0);
  } catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }

    return 0;
}
