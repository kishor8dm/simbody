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

#include "Simbody.h"

namespace {
//-0.1504133023, -0.7202261496  0.6643972918  0.1312492688
//-0.8802997333,  -0.4219333635  -0.1958447339 -0.09321725385
//0.1707965794, -0.1900937466 -0.1176534344 -0.9596095901
//0.08860797748, 0.4014763757 0.3004798435 0.8606260569
//0.3201724846,  0.7169894904  0.3068382102 -0.5378345131
}

using namespace SimTK;
using std::cout; using std::endl;

int main() {
  try {   
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Quaternion q0(-0.1504133023, -0.7202261496,  0.6643972918,  0.1312492688);
    Quaternion q1(-0.8802997333,  -0.4219333635,  -0.1958447339, -0.09321725385);
    Quaternion q2(0.1707965794, -0.1900937466, -0.1176534344, -0.9596095901);
    Quaternion q3(0.3201724846,  0.7169894904,  0.3068382102, -0.5378345131);
    Quaternion q4{-0.2366251647, 0.3815343678,  0.845472455,   0.28916502};
    Quaternion q5{-0.8099961877, 0.09318492562,  0.5701901317, -0.1005286202};
    Quaternion q6{0.09452920407, -0.6611310244, -0.5389370918, -0.5133388638};
    Quaternion q7{-0.06702885032, -0.7941805124,   -0.17857261,  0.5769716501};
    Quaternion halfPi(0.7071067691, 0.7071067691, 0, 0);
    Quaternion eye(1, 0, 0, 0);

    Mat33 I(1,0,0,0,2,0,0,0,3);
    Rotation Q1(q1);
    I = (~Q1) * I * Q1;
    std::cout<<I<< "\n";

    Force::UniformGravity gravity(forces, matter, Vec3(0, Real(-10.0), 0));
    Body::Rigid pendulumBody(MassProperties(10.0, Vec3(0), Inertia(I)));

    MobilizedBody::Pin pendulum0(matter.updGround(),
                                Transform(Rotation(q0), Vec3(0,-1,0)),
                                pendulumBody, 
                                Transform(Rotation(q1),Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum1(pendulum0,
                                Transform(Rotation(q2), Vec3(0,-1,0)),
                                pendulumBody, 
                                Transform(Rotation(q3), Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum2(pendulum1,
                               Transform(Rotation(q4), Vec3(0,-1,0)),
                               pendulumBody,
                               Transform(Rotation(q5), Vec3(0, 1, 0)));
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum0.setOneQ(state, 0, Pi/4);
    pendulum1.setOneQ(state, 0, Pi/4);
    pendulum2.setOneQ(state, 0, Pi/4);

    pendulum0.setOneU(state, 0, 1);
    pendulum1.setOneU(state, 0, 1);
    pendulum2.setOneU(state, 0, 1);

    system.realize(state);
    std::cout << pendulum0.getBodyRotation(state).convertRotationToQuaternion() << "\n";
    std::cout << pendulum1.getBodyRotation(state).convertRotationToQuaternion() << "\n";
    std::cout << pendulum0.getBodyOriginLocation(state) << "\n";
    std::cout << pendulum1.getBodyOriginLocation(state) << "\n";
    std::cout << state.getUDot() << "\n";

  } catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }

    return 0;
}
