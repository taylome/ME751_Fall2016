// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mike Taylor
// =============================================================================
//
// ME751- HW8 - Problem#3
//
// http://lim.ii.udc.es/mbsbenchmark/dist/A03/A03_specification.xml
//
// Recall that Irrlicht uses a left-hand frame, so everything is rendered with
// left and right flipped.
//
// =============================================================================

#include <ostream>
#include <fstream>

#include "chrono/core/ChFileutils.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChBody.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"

#include "chrono_irrlicht/ChIrrApp.h"

#include "chrono_postprocess/ChPovRay.h"

#include "chrono/assets/ChPointPointDrawing.h"

#include "ChronoME751_config.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace irr;
using namespace postprocess;

// =============================================================================
// Local variables
//
static const std::string val_dir = "../RESULTS/";
static const std::string out_dir = val_dir + "HW8/";

// =============================================================================
// Prototypes of local functions
//
utils::CSV_writer OutStream();

// =============================================================================
//
// Main driver function for running the simulation and validating the results.
//
int main(int argc, char* argv[])
{
  bool animate = (argc > 1);
  bool save = (argc > 2);

  // Set the path to the Chrono data folder
  SetChronoDataPath(CHRONO_DATA_DIR);

  // Create output directory (if it does not already exist)
  if (ChFileutils::MakeDirectory(val_dir.c_str()) < 0) {
    std::cout << "Error creating directory " << val_dir << std::endl;
    return 1;
  }
  if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
    std::cout << "Error creating directory " << out_dir << std::endl;
    return 1;
  }

  // Set the simulation and output step sizes
  double simTimeStep = 5e-6;
  double outTimeStep = 1e-4;


  // Settings
  //---------

  // There are no units in Chrono, so values must be consistent
  // (MKS is used in this example)
  double g = 9.81*0;
  double timeRecord = 0.05;                 // simulation length
  
  // Initial Locations
  double Beta0 = -0.0620;
  ChVector<>PntOLoc(0,0,0);
  ChVector<>PntALoc(-0.06934, -0.00227, 0);
  ChVector<>PntBLoc(-0.03635, 0.03273, 0);
  ChVector<>PntCLoc(0.01400, 0.07200, 0);
  ChVector<>PntFLoc(0.007*std::cos(Beta0), 0.007*std::sin(Beta0), 0);

  double a = 0.035;
  double c = 0.028;
  ChVector<>Temp = PntFLoc - PntBLoc;
  double b = Temp.Length();
  double gamma = std::acos((a*a + b*b - c*c) / (2 * a*b)); //Angle Between BF & BE
  double alpha = atan2(Temp.y, Temp.x)- gamma; //Angle between horizontal and BE
  ChVector<>PntELoc(PntBLoc.x+0.035*std::cos(alpha), PntBLoc.y + 0.035*std::sin(alpha), PntBLoc.z + 0);

  b = std::sqrt(0.02*0.02 + 0.017*0.017);
  c = std::sqrt(0.02*0.02 + (0.035 - 0.017)*(0.035 - 0.017));
  gamma = std::acos((a*a + b*b - c*c) / (2 * a*b)); //Angle Between BD & BE
  alpha += gamma; //Angle between horizontal and BD
  ChVector<>PntDLoc(PntBLoc.x + b*std::cos(alpha), PntBLoc.y + b*std::sin(alpha), PntBLoc.z + 0);

  Temp = PntELoc - PntALoc;
  a = Temp.Length();
  b = 0.04;
  c = 0.02;
  gamma = std::acos((a*a + b*b - c*c) / (2 * a*b)); //Angle Between AE & AG
  alpha = atan2(Temp.y, Temp.x) + gamma; //Angle between horizontal and AG
  ChVector<>PntGLoc(PntALoc.x + b*std::cos(alpha), PntALoc.y + b*std::sin(alpha), PntALoc.z + 0);
  

  alpha = atan2(Temp.y, Temp.x) - gamma; //Angle between horizontal and AH
  ChVector<>PntHLoc(PntALoc.x + b*std::cos(alpha), PntALoc.y + b*std::sin(alpha), PntALoc.z + 0);
  

  std::cout << "PntO: " << PntOLoc.x << " " << PntOLoc.y << " " << PntOLoc.z << std::endl;
  std::cout << "PntA: " << PntALoc.x << " " << PntALoc.y << " " << PntALoc.z << std::endl;
  std::cout << "PntB: " << PntBLoc.x << " " << PntBLoc.y << " " << PntBLoc.z << std::endl;
  std::cout << "PntC: " << PntCLoc.x << " " << PntCLoc.y << " " << PntCLoc.z << std::endl;
  std::cout << "PntD: " << PntDLoc.x << " " << PntDLoc.y << " " << PntDLoc.z << std::endl;
  std::cout << "PntE: " << PntELoc.x << " " << PntELoc.y << " " << PntELoc.z << std::endl;
  std::cout << "PntF: " << PntFLoc.x << " " << PntFLoc.y << " " << PntFLoc.z << std::endl;
  std::cout << "PntG: " << PntGLoc.x << " " << PntGLoc.y << " " << PntGLoc.z << std::endl;
  std::cout << "PntH: " << PntHLoc.x << " " << PntHLoc.y << " " << PntHLoc.z << std::endl;



  double RevGraphicRad = 0.001;
  double RevGraphicLen = 0.003;
  double RodGraphicRad = 0.002;

  double smallInertia = 1e-8;

  //Rod1 - OF - Length 0.007
  double Rod1Length = 0.007;
  double Rod1Mass = 0.04325;
  ChVector<>Rod1Axis = PntFLoc - PntOLoc;
  Rod1Axis.Normalize();
  ChVector<>Rod1CG = PntOLoc + 9.2E-004*Rod1Axis;
  ChVector<>Rod1Ixx(smallInertia, smallInertia, 2.194e-06);
  
  //Rod2 - EF - Length 0.028
  double Rod2Length = 0.028;
  double Rod2Mass = 0.00365;
  ChVector<>Rod2Axis = PntELoc - PntFLoc;
  Rod2Axis.Normalize();
  ChVector<>Rod2CG = PntFLoc + 0.0115*Rod2Axis;
  ChVector<>Rod2Ixx(smallInertia, smallInertia, 4.41E-007);

  //Triangle3 - BDE
  double Triange3Mass = 0.02373;
  ChVector<>Triange3YAxis = PntBLoc - PntELoc;
  Triange3YAxis.Normalize();
  ChVector<>Triange3XAxis(0, 0, 0);
  Triange3XAxis.Cross(Triange3YAxis, ChVector<>(0, 0, 1));
  Triange3XAxis.Normalize();
  ChVector<>Triange3CG = PntBLoc + 0.01043*Triange3XAxis + -0.01874*Triange3YAxis;
  ChVector<>Triange3Ixx(smallInertia, smallInertia, 5.255E-006);

  //Rod4 - EG - Length 0.02
  double Rod4Mass = 0.00706;
  ChVector<>Rod4Axis = PntGLoc - PntELoc;
  Rod4Axis.Normalize();
  ChVector<>Rod4CG = PntELoc + 0.01421*Rod4Axis;
  ChVector<>Rod4Ixx(smallInertia, smallInertia, 5.667E-007);

  //Rod5 - AG - Length 0.04
  double Rod5Mass = 0.0705;
  ChVector<>Rod5XAxis = PntGLoc - PntALoc;
  Rod5XAxis.Normalize();
  ChVector<>Rod5YAxis(0, 0, 0);
  Rod5YAxis.Cross(ChVector<>(0, 0, 1), Rod5XAxis);
  Rod5YAxis.Normalize();
  ChVector<>Rod5CG = PntALoc + 0.02308*Rod5XAxis + 0.00916*Rod5YAxis;
  ChVector<>Rod5Ixx(smallInertia, smallInertia, 1.169E-005);

  //Rod6 - EH - Length 0.02
  double Rod6Mass = 0.00706;
  ChVector<>Rod6Axis = PntHLoc - PntELoc;
  Rod6Axis.Normalize();
  ChVector<>Rod6CG = PntELoc + 0.01421*Rod4Axis;
  ChVector<>Rod6Ixx(smallInertia, smallInertia, 5.667E-007);

  //Rod7 - AG - Length 0.04
  double Rod7Mass = 0.05498;
  ChVector<>Rod7YAxis = PntALoc - PntHLoc;
  Rod7YAxis.Normalize();
  ChVector<>Rod7XAxis(0, 0, 0);
  Rod7XAxis.Cross(Rod7YAxis, ChVector<>(0, 0, 1));
  Rod7XAxis.Normalize();
  ChVector<>Rod7CG = PntALoc + -0.00449*Rod7XAxis + -0.01228*Rod7YAxis;
  ChVector<>Rod7Ixx(smallInertia, smallInertia, 1.912E-005);

  // Create the mechanical system
  // ----------------------------

  // Create a ChronoENGINE physical system: all bodies and constraints will be
  // handled by this ChSystem object.

  ChSystem my_system;
  my_system.Set_G_acc(ChVector<>(0.0, -g, 0.0));

  my_system.SetIntegrationType(ChSystem::INT_EULER_IMPLICIT_LINEARIZED);
  my_system.SetMaxItersSolverSpeed(100);
  my_system.SetMaxItersSolverStab(100); //Tasora stepper uses this, Anitescu does not
  my_system.SetSolverType(ChSystem::SOLVER_SOR);
  my_system.SetTol(1e-8);
  my_system.SetTolForce(1e-6);

  // Create the ground body
  auto ground = std::make_shared<ChBody>();
  my_system.AddBody(ground);
  ground->SetBodyFixed(true);
  // Add some geometry to the ground body for visualizing the Point A revolute joint
  auto cyl_gA = std::make_shared<ChCylinderShape>();
  cyl_gA->GetCylinderGeometry().p1 = PntALoc + ChVector<>(0, 0, RevGraphicLen);
  cyl_gA->GetCylinderGeometry().p2 = PntALoc - ChVector<>(0, 0, RevGraphicLen);
  cyl_gA->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gA);
  // Add some geometry to the ground body for visualizing the Point B revolute joint
  auto cyl_gB = std::make_shared<ChCylinderShape>();
  cyl_gB->GetCylinderGeometry().p1 = PntBLoc + ChVector<>(0, 0, RevGraphicLen);
  cyl_gB->GetCylinderGeometry().p2 = PntBLoc - ChVector<>(0, 0, RevGraphicLen);
  cyl_gB->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gB);
  // Add some geometry to the ground body for visualizing the Point C revolute joint
  auto cyl_gC = std::make_shared<ChCylinderShape>();
  cyl_gC->GetCylinderGeometry().p1 = PntCLoc + ChVector<>(0, 0, RevGraphicLen);
  cyl_gC->GetCylinderGeometry().p2 = PntCLoc - ChVector<>(0, 0, RevGraphicLen);
  cyl_gC->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gC);
  // Add some geometry to the ground body for visualizing the Point O revolute joint
  auto cyl_gO = std::make_shared<ChCylinderShape>();
  cyl_gO->GetCylinderGeometry().p1 = PntOLoc + ChVector<>(0, 0, RevGraphicLen);
  cyl_gO->GetCylinderGeometry().p2 = PntOLoc - ChVector<>(0, 0, RevGraphicLen);
  cyl_gO->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gO);

  //auto cyl_gD = std::make_shared<ChCylinderShape>();
  //cyl_gD->GetCylinderGeometry().p1 = PntDLoc + ChVector<>(0, 0, RevGraphicLen);
  //cyl_gD->GetCylinderGeometry().p2 = PntDLoc - ChVector<>(0, 0, RevGraphicLen);
  //cyl_gD->GetCylinderGeometry().rad = RevGraphicRad;
  //ground->AddAsset(cyl_gD);
  //auto cyl_gE = std::make_shared<ChCylinderShape>();
  //cyl_gE->GetCylinderGeometry().p1 = PntELoc + ChVector<>(0, 0, RevGraphicLen);
  //cyl_gE->GetCylinderGeometry().p2 = PntELoc - ChVector<>(0, 0, RevGraphicLen);
  //cyl_gE->GetCylinderGeometry().rad = RevGraphicRad;
  //ground->AddAsset(cyl_gE);
  //auto cyl_gF = std::make_shared<ChCylinderShape>();
  //cyl_gF->GetCylinderGeometry().p1 = PntFLoc + ChVector<>(0, 0, RevGraphicLen);
  //cyl_gF->GetCylinderGeometry().p2 = PntFLoc - ChVector<>(0, 0, RevGraphicLen);
  //cyl_gF->GetCylinderGeometry().rad = RevGraphicRad;
  //ground->AddAsset(cyl_gF);
  //auto cyl_gG = std::make_shared<ChCylinderShape>();
  //cyl_gG->GetCylinderGeometry().p1 = PntGLoc + ChVector<>(0, 0, RevGraphicLen);
  //cyl_gG->GetCylinderGeometry().p2 = PntGLoc - ChVector<>(0, 0, RevGraphicLen);
  //cyl_gG->GetCylinderGeometry().rad = RevGraphicRad;
  //ground->AddAsset(cyl_gG);
  //auto cyl_gH = std::make_shared<ChCylinderShape>();
  //cyl_gH->GetCylinderGeometry().p1 = PntHLoc + ChVector<>(0, 0, RevGraphicLen);
  //cyl_gH->GetCylinderGeometry().p2 = PntHLoc - ChVector<>(0, 0, RevGraphicLen);
  //cyl_gH->GetCylinderGeometry().rad = RevGraphicRad;
  //ground->AddAsset(cyl_gH);



  // Create the Rod1 body in an initial configuration
  auto Rod1 = std::make_shared<ChBody>();
  my_system.AddBody(Rod1);
  Rod1->SetPos(Rod1CG);
  Rod1->SetRot(QUNIT);
  Rod1->SetMass(Rod1Mass);
  Rod1->SetInertiaXX(Rod1Ixx);
  // Add some geometry for visualization
  auto cyl_Rod1 = std::make_shared<ChCylinderShape>();
  cyl_Rod1->GetCylinderGeometry().p1 = PntOLoc- Rod1CG;
  cyl_Rod1->GetCylinderGeometry().p2 = PntFLoc- Rod1CG;
  cyl_Rod1->GetCylinderGeometry().rad = RodGraphicRad;
  Rod1->AddAsset(cyl_Rod1);

  // Create the Rod2 body in an initial configuration
  auto Rod2 = std::make_shared<ChBody>();
  my_system.AddBody(Rod2);
  Rod2->SetPos(Rod2CG);
  Rod2->SetRot(QUNIT);
  Rod2->SetMass(Rod2Mass);
  Rod2->SetInertiaXX(Rod2Ixx);
  // Add some geometry for visualization
  auto cyl_Rod2 = std::make_shared<ChCylinderShape>();
  cyl_Rod2->GetCylinderGeometry().p1 = PntFLoc- Rod2CG;
  cyl_Rod2->GetCylinderGeometry().p2 = PntELoc- Rod2CG;
  cyl_Rod2->GetCylinderGeometry().rad = RodGraphicRad;
  Rod2->AddAsset(cyl_Rod2);

  // Create the Triangle3 body in an initial configuration
  auto Triange3 = std::make_shared<ChBody>();
  my_system.AddBody(Triange3);
  Triange3->SetPos(Triange3CG);
  Triange3->SetRot(QUNIT);
  Triange3->SetMass(Triange3Mass);
  Triange3->SetInertiaXX(Triange3Ixx);
  // Add some geometry for visualization
  auto cyl_Triange3BE = std::make_shared<ChCylinderShape>();
  cyl_Triange3BE->GetCylinderGeometry().p1 = PntBLoc- Triange3CG;
  cyl_Triange3BE->GetCylinderGeometry().p2 = PntELoc- Triange3CG;
  cyl_Triange3BE->GetCylinderGeometry().rad = RodGraphicRad;
  Triange3->AddAsset(cyl_Triange3BE);
  auto cyl_Triange3DE = std::make_shared<ChCylinderShape>();
  cyl_Triange3DE->GetCylinderGeometry().p1 = PntDLoc- Triange3CG;
  cyl_Triange3DE->GetCylinderGeometry().p2 = PntELoc- Triange3CG;
  cyl_Triange3DE->GetCylinderGeometry().rad = RodGraphicRad;
  Triange3->AddAsset(cyl_Triange3DE);
  auto cyl_Triange3BD = std::make_shared<ChCylinderShape>();
  cyl_Triange3BD->GetCylinderGeometry().p1 = PntBLoc- Triange3CG;
  cyl_Triange3BD->GetCylinderGeometry().p2 = PntDLoc- Triange3CG;
  cyl_Triange3BD->GetCylinderGeometry().rad = RodGraphicRad;
  Triange3->AddAsset(cyl_Triange3BD);

  // Create the Rod4 body in an initial configuration
  auto Rod4 = std::make_shared<ChBody>();
  my_system.AddBody(Rod4);
  Rod4->SetPos(Rod4CG);
  Rod4->SetRot(QUNIT);
  Rod4->SetMass(Rod4Mass);
  Rod4->SetInertiaXX(Rod4Ixx);
  // Add some geometry for visualization
  auto cyl_Rod4 = std::make_shared<ChCylinderShape>();
  cyl_Rod4->GetCylinderGeometry().p1 = PntGLoc- Rod4CG;
  cyl_Rod4->GetCylinderGeometry().p2 = PntELoc- Rod4CG;
  cyl_Rod4->GetCylinderGeometry().rad = RodGraphicRad;
  Rod4->AddAsset(cyl_Rod4);

  // Create the Rod5 body in an initial configuration
  auto Rod5 = std::make_shared<ChBody>();
  my_system.AddBody(Rod5);
  Rod5->SetPos(Rod5CG);
  Rod5->SetRot(QUNIT);
  Rod5->SetMass(Rod5Mass);
  Rod5->SetInertiaXX(Rod5Ixx);
  // Add some geometry for visualization
  auto cyl_Rod5 = std::make_shared<ChCylinderShape>();
  cyl_Rod5->GetCylinderGeometry().p1 = PntGLoc- Rod5CG;
  cyl_Rod5->GetCylinderGeometry().p2 = PntALoc- Rod5CG;
  cyl_Rod5->GetCylinderGeometry().rad = RodGraphicRad;
  Rod5->AddAsset(cyl_Rod5);

  // Create the Rod6 body in an initial configuration
  auto Rod6 = std::make_shared<ChBody>();
  my_system.AddBody(Rod6);
  Rod6->SetPos(Rod6CG);
  Rod6->SetRot(QUNIT);
  Rod6->SetMass(Rod6Mass);
  Rod6->SetInertiaXX(Rod6Ixx);
  // Add some geometry for visualization
  auto cyl_Rod6 = std::make_shared<ChCylinderShape>();
  cyl_Rod6->GetCylinderGeometry().p1 = PntELoc- Rod6CG;
  cyl_Rod6->GetCylinderGeometry().p2 = PntHLoc- Rod6CG;
  cyl_Rod6->GetCylinderGeometry().rad = RodGraphicRad;
  Rod6->AddAsset(cyl_Rod6);

  // Create the Rod6 body in an initial configuration
  auto Rod7 = std::make_shared<ChBody>();
  my_system.AddBody(Rod7);
  Rod7->SetPos(Rod7CG);
  Rod7->SetRot(QUNIT);
  Rod7->SetMass(Rod7Mass);
  Rod7->SetInertiaXX(Rod7Ixx);
  // Add some geometry for visualization
  auto cyl_Rod7 = std::make_shared<ChCylinderShape>();
  cyl_Rod7->GetCylinderGeometry().p1 = PntALoc - Rod7CG;
  cyl_Rod7->GetCylinderGeometry().p2 = PntHLoc - Rod7CG;
  cyl_Rod7->GetCylinderGeometry().rad = RodGraphicRad;
  Rod7->AddAsset(cyl_Rod7);

  //// Create revolute or spherical joints between the parts
  //// The revolute joint's axis of rotation will be the Z axis
  //// of the specified rotation matrix.

  //auto revoluteJointO = std::make_shared<ChLinkLockRevolute>();
  //revoluteJointO->Initialize(Rod1, ground, ChCoordsys<>(PntOLoc, QUNIT));
  //my_system.AddLink(revoluteJointO);

  auto revoluteJointA5 = std::make_shared<ChLinkLockRevolute>();
  revoluteJointA5->Initialize(Rod5, ground, ChCoordsys<>(PntALoc, QUNIT));
  my_system.AddLink(revoluteJointA5);

  auto revoluteJointA7 = std::make_shared<ChLinkLockRevolute>();
  revoluteJointA7->Initialize(Rod7, ground, ChCoordsys<>(PntALoc, QUNIT));
  my_system.AddLink(revoluteJointA7);

  auto revoluteJointB = std::make_shared<ChLinkLockRevolute>();
  revoluteJointB->Initialize(Triange3, ground, ChCoordsys<>(PntBLoc, QUNIT));
  my_system.AddLink(revoluteJointB);

  auto revoluteJointE2 = std::make_shared<ChLinkLockRevolute>();
  revoluteJointE2->Initialize(Rod2, Triange3, ChCoordsys<>(PntELoc, QUNIT));
  my_system.AddLink(revoluteJointE2);

  auto revoluteJointE4 = std::make_shared<ChLinkLockRevolute>();
  revoluteJointE4->Initialize(Rod4, Triange3, ChCoordsys<>(PntELoc, QUNIT));
  my_system.AddLink(revoluteJointE4);

  auto revoluteJointE6 = std::make_shared<ChLinkLockRevolute>();
  revoluteJointE6->Initialize(Rod6, Triange3, ChCoordsys<>(PntELoc, QUNIT));
  my_system.AddLink(revoluteJointE6);

  auto revoluteJointF = std::make_shared<ChLinkLockRevolute>();
  revoluteJointF->Initialize(Rod2, Rod1, ChCoordsys<>(PntFLoc, QUNIT));
  my_system.AddLink(revoluteJointF);

  auto revoluteJointG = std::make_shared<ChLinkLockRevolute>();
  revoluteJointG->Initialize(Rod5, Rod4, ChCoordsys<>(PntGLoc, QUNIT));
  my_system.AddLink(revoluteJointG);

  auto revoluteJointH = std::make_shared<ChLinkLockRevolute>();
  revoluteJointH->Initialize(Rod7, Rod6, ChCoordsys<>(PntHLoc, QUNIT));
  my_system.AddLink(revoluteJointH);


  // Create the spring between Triange3 and ground. The spring end points are
  // specified in the ground frame.
  auto spring_1 = std::make_shared<ChLinkSpring>();
  spring_1->Initialize(Triange3, ground, false, PntDLoc, PntCLoc, false, 0.07785);
  spring_1->Set_SpringK(4530);
  spring_1->Set_SpringR(0);
  my_system.AddLink(spring_1);
  // Attach a visualization asset.
  spring_1->AddAsset(std::make_shared<ChPointPointSpring>(0.005, 80, 15));

  // Add an engine at Point O (also acts as a revolute joint)
  auto link_engine = std::make_shared<ChLinkEngine>();
  link_engine->Initialize(Rod1, ground, ChCoordsys<>(ChVector<>(0, 0, 0), QUNIT));
  link_engine->Set_shaft_mode(ChLinkEngine::ENG_SHAFT_LOCK);  // also works as revolute support
  link_engine->Set_eng_mode(ChLinkEngine::ENG_MODE_TORQUE);
  if (auto mfun = std::dynamic_pointer_cast<ChFunction_Const>(link_engine->Get_tor_funct()))
	  mfun->Set_yconst(0.033);  // Nm
  my_system.AddLink(link_engine);




  // Perform the simulation (animation with Irrlicht option)
  // -------------------------------------------------------

  if (animate)
  {
      ChVector<> CameraPnt = PntELoc;
      // Create the Irrlicht application for visualization
      ChIrrApp * application = new ChIrrApp(&my_system, L"ME751 - HW8 - Problem 3", core::dimension2d<u32>(800, 600), false, true);
      application->AddTypicalLogo();
      application->AddTypicalSky();
      application->AddTypicalLights();
      core::vector3df lookat((f32)CameraPnt.x, (f32)CameraPnt.y, (f32)CameraPnt.z);
      application->AddTypicalCamera(lookat + core::vector3df(0, -.05, -.2), lookat);

      // Now have the visulization tool (Irrlicht) create its geometry from the
      // assets defined above
      application->AssetBindAll();
      application->AssetUpdateAll();

      application->SetTimestep(simTimeStep);

      // Simulation loop
      double outTime = 0;
      int    outFrame = 1;
      double time = 0;


      // Create an exporter to POVray
      ChPovRay pov_exporter = ChPovRay(&my_system);

      if (save) {

          // ==Asset== Attach a video camera. This will be used by Irrlicht,
          // or POVray postprocessing, etc. Note that a camera can also be
          // put in a moving object.
          auto mcamera = std::make_shared<ChCamera>();
          mcamera->SetAngle(50);
          mcamera->SetPosition(ChVector<>(0, 0.05, 0.2));
          mcamera->SetAimPoint(ChVector<>(0, 0, 0));
          ground->AddAsset(mcamera);


          // Sets some file names for in-out processes.
          pov_exporter.SetTemplateFile("../data/_template_POV.pov");
          pov_exporter.SetOutputScriptFile("rendering_frames.pov");
          pov_exporter.SetOutputDataFilebase("my_state");
          pov_exporter.SetPictureFilebase("picture");

          // Even better: save the .dat files and the .bmp files
          // in two subdirectories, to avoid cluttering the current
          // directory...
          ChFileutils::MakeDirectory("output");
          ChFileutils::MakeDirectory("anim");
          pov_exporter.SetOutputDataFilebase("output/my_state");
          pov_exporter.SetPictureFilebase("anim/picture");

          // Tell to the POVray exporter that
          // he must take care of converting the shapes of
          // all items!
          pov_exporter.AddAll();

          // 1) Create the two .pov and .ini files for POV-Ray (this must be done
          //    only once at the beginning of the simulation).
          pov_exporter.ExportScript();
      }
      while (application->GetDevice()->run())
      {
          time = my_system.GetChTime();
          // End simulation
          if ((time > timeRecord) && (timeRecord > 0))
              break;

          if (save && my_system.GetChTime() >= outTime - simTimeStep / 2) {
              //char filename[100];
              //sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), outFrame);
              //utils::WriteShapesPovray(&my_system, filename);
              outTime += outTimeStep;
              outFrame++;
              pov_exporter.ExportData();
          }

          application->BeginScene();
          application->DrawAll();

          //// Draw an XZ grid at the global origin to add in visualization
          //ChIrrTools::drawGrid(
          //    application->GetVideoDriver(), 1, 1, 20, 20,
          //    ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(CH_C_PI_2)),
          //    video::SColor(255, 80, 100, 100), true);

          application->DoStep();  //Take one step in time
          application->EndScene();
      }

      return true;
  }

  // Perform the simulation (record results option)
  // ------------------------------------------------

  // Create the CSV_Writer output objects (TAB delimited)
  utils::CSV_writer out_pos_F = OutStream();
  utils::CSV_writer out_pos_D = OutStream();
  utils::CSV_writer out_vel_D = OutStream();
  utils::CSV_writer out_SprFrc = OutStream();


  // Perform a system assembly to ensure we have the correct accelerations at
  // the initial time.
  my_system.DoFullAssembly();


  // Simulation loop
  double simTime = 0;
  double outTime = 0;

  //timeRecord = .0001;
  while (simTime <= timeRecord + simTimeStep / 2)
  {
      // Ensure that the final data point is recorded.
      if (simTime >= outTime - simTimeStep / 2)
      {

          // position & velocity (expressed in global frame).
          ChVector<> &Rod1pos = Rod1->GetFrame_REF_to_abs().TransformPointLocalToParent(PntFLoc-Rod1CG);
          ChVector<> &Triange3Dpos = Triange3->GetFrame_REF_to_abs().TransformPointLocalToParent(PntDLoc-Triange3CG);
          ChVector<> &Triange3Dvel = Triange3->GetFrame_REF_to_abs().PointSpeedLocalToParent(PntDLoc - Triange3CG);

		  out_pos_F << simTime << Rod1pos << std::endl;
		  out_pos_D << simTime << Triange3Dpos << std::endl;
		  out_vel_D << simTime << Triange3Dvel << std::endl;

		  //Spring Reaction Force 
		  out_SprFrc << simTime << spring_1->Get_SpringReact() << std::endl;


          // Increment output time
          outTime += outTimeStep;
      }

      // Advance simulation by one step
      my_system.DoStepDynamics(simTimeStep);

      // Increment simulation time
      simTime += simTimeStep;
  }

  // Write output files
  out_pos_F.write_to_file(out_dir + "HW8P3_CHRONO_F_Pos.txt", "");
  out_pos_D.write_to_file(out_dir + "HW8P3_CHRONO_D_Pos.txt", "");
  out_vel_D.write_to_file(out_dir + "HW8P3_CHRONO_D_Vel.txt", "");
  out_SprFrc.write_to_file(out_dir + "HW8P3_CHRONO_SprFrc.txt", "");


  return true;
}



// =============================================================================
//
// Utility function to create a CSV output stream and set output format options.
//
utils::CSV_writer OutStream()
{
    utils::CSV_writer out("\t");

    out.stream().setf(std::ios::scientific | std::ios::showpos);
    out.stream().precision(6);

    return out;
}
