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
// ME751- HW9 - Problem#4
//
// http://lim.ii.udc.es/mbsbenchmark/dist/A04/A04_specification.xml
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
static const std::string out_dir = val_dir + "HW9/";

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
  double simTimeStep = 5e-4;
  double outTimeStep = 1e-2; 
  //double outTimeStep = 1.0 / 30.0;


  // Settings
  //---------

  // There are no units in Chrono, so values must be consistent
  // (MKS is used in this example)
  double g = 9.81;
  double timeRecord = 10.0;                 // simulation length
  
  double mass = 1.0;                     // mass of link
  double length = 1.0;                   // length of link
  ChVector<> inertiaXX(0.0833333, 0.0001, 0.0833333);  // mass moments of inertia of pendulum (centroidal frame)

  double RevGraphicRad = 0.05;
  double RevGraphicLen = 0.1;
  double RodGraphicRad = 0.05;

  double smallInertia = 1e-8;

  //Inital Joint Locations
  ChVector<> PntALoc(0,0,0);
  ChQuaternion<> PntARot;
  PntARot.Q_from_AngX(-CH_C_PI_2);

  ChVector<> PntBLoc(1, 0, 0);
  ChQuaternion<> PntBRot = QUNIT;

  ChVector<> PntCLoc(1, -1, 0);
  ChQuaternion<> PntCRot;
  PntCRot.Q_from_AngY(CH_C_PI_2);

  ChVector<> PntDLoc(1, -1, 1);
  ChQuaternion<> PntDRot;
  PntDRot.Q_from_AngX(-CH_C_PI_2);

  ChVector<> PntELoc(0, -1, 1);
  ChQuaternion<> PntERot = QUNIT;

  ChVector<> PntFLoc(0, 0, 1);
  ChQuaternion<> PntFRot;
  PntFRot.Q_from_AngY(CH_C_PI_2);


  //Inital Rod Locations & Rotations
  auto RodABCG = 0.5*(PntALoc + PntBLoc);
  ChQuaternion<> RodABRot;
  RodABRot.Q_from_AngY(CH_C_PI_2);

  auto RodBCCG = 0.5*(PntBLoc + PntCLoc);
  ChQuaternion<> RodBCRot;
  RodBCRot.Q_from_AngX(-CH_C_PI_2);

  auto RodCDCG = 0.5*(PntCLoc + PntDLoc);
  ChQuaternion<> RodCDRot = QUNIT;

  auto RodDECG = 0.5*(PntDLoc + PntELoc);
  ChQuaternion<> RodDERot;
  RodDERot.Q_from_AngY(CH_C_PI_2);

  auto RodEFCG = 0.5*(PntELoc + PntFLoc);
  ChQuaternion<> RodEFRot;
  RodEFRot.Q_from_AngX(-CH_C_PI_2);



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
  cyl_gA->GetCylinderGeometry().p1 = PntALoc + PntARot.Rotate(ChVector<>(0, 0, RevGraphicLen));
  cyl_gA->GetCylinderGeometry().p2 = PntALoc - PntARot.Rotate(ChVector<>(0, 0, RevGraphicLen));
  cyl_gA->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gA);
  // Add some geometry to the ground body for visualizing the Point F revolute joint
  auto cyl_gF = std::make_shared<ChCylinderShape>();
  cyl_gF->GetCylinderGeometry().p1 = PntFLoc + PntFRot.Rotate(ChVector<>(0, 0, RevGraphicLen));
  cyl_gF->GetCylinderGeometry().p2 = PntFLoc - PntFRot.Rotate(ChVector<>(0, 0, RevGraphicLen));
  cyl_gF->GetCylinderGeometry().rad = RevGraphicRad;
  ground->AddAsset(cyl_gF);


  // Create the RodAB body in an initial configuration
  auto RodAB = std::make_shared<ChBody>();
  my_system.AddBody(RodAB);
  RodAB->SetPos(RodABCG);
  RodAB->SetRot(RodABRot);
  RodAB->SetMass(mass);
  RodAB->SetInertiaXX(inertiaXX);
  // Add some geometry for visualization
  auto cyl_RodAB = std::make_shared<ChCylinderShape>();
  cyl_RodAB->GetCylinderGeometry().p1 = RodABRot.Rotate(PntALoc - RodABCG);
  cyl_RodAB->GetCylinderGeometry().p2 = RodABRot.Rotate(PntBLoc - RodABCG);
  cyl_RodAB->GetCylinderGeometry().rad = RodGraphicRad;
  RodAB->AddAsset(cyl_RodAB);
  auto cyl_B = std::make_shared<ChCylinderShape>();
  cyl_B->GetCylinderGeometry().p1 = ChVector<>(0, 0, .5*length) + RodABRot.Rotate(PntBRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_B->GetCylinderGeometry().p2 = ChVector<>(0, 0, .5*length) - RodABRot.Rotate(PntBRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_B->GetCylinderGeometry().rad = RodGraphicRad;
  RodAB->AddAsset(cyl_B);

  // Create the RodBC body in an initial configuration
  auto RodBC = std::make_shared<ChBody>();
  my_system.AddBody(RodBC);
  RodBC->SetPos(RodBCCG);
  RodBC->SetRot(RodBCRot);
  RodBC->SetMass(mass);
  RodBC->SetInertiaXX(inertiaXX);
  // Add some geometry for visualization
  auto cyl_RodBC = std::make_shared<ChCylinderShape>();
  cyl_RodBC->GetCylinderGeometry().p1 = RodBCRot.Rotate(PntBLoc - RodBCCG);
  cyl_RodBC->GetCylinderGeometry().p2 = RodBCRot.Rotate(PntCLoc - RodBCCG);
  cyl_RodBC->GetCylinderGeometry().rad = RodGraphicRad;
  RodBC->AddAsset(cyl_RodBC);
  auto cyl_C = std::make_shared<ChCylinderShape>();
  cyl_C->GetCylinderGeometry().p1 = ChVector<>(0, 0, -.5*length) + RodBCRot.Rotate(PntCRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_C->GetCylinderGeometry().p2 = ChVector<>(0, 0, -.5*length) - RodBCRot.Rotate(PntCRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_C->GetCylinderGeometry().rad = RodGraphicRad;
  RodBC->AddAsset(cyl_C);

  // Create the RodCD body in an initial configuration
  auto RodCD = std::make_shared<ChBody>();
  my_system.AddBody(RodCD);
  RodCD->SetPos(RodCDCG);
  RodCD->SetRot(RodCDRot);
  RodCD->SetMass(mass);
  RodCD->SetInertiaXX(inertiaXX);
  // Add some geometry for visualization
  auto cyl_RodCD = std::make_shared<ChCylinderShape>();
  cyl_RodCD->GetCylinderGeometry().p1 = RodCDRot.Rotate(PntCLoc - RodCDCG);
  cyl_RodCD->GetCylinderGeometry().p2 = RodCDRot.Rotate(PntDLoc - RodCDCG);
  cyl_RodCD->GetCylinderGeometry().rad = RodGraphicRad;
  RodCD->AddAsset(cyl_RodCD);
  auto cyl_D = std::make_shared<ChCylinderShape>();
  cyl_D->GetCylinderGeometry().p1 = ChVector<>(0, 0, .5*length) + RodCDRot.Rotate(PntDRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_D->GetCylinderGeometry().p2 = ChVector<>(0, 0, .5*length) - RodCDRot.Rotate(PntDRot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_D->GetCylinderGeometry().rad = RodGraphicRad;
  RodCD->AddAsset(cyl_D);


  // Create the RodDE body in an initial configuration
  auto RodDE = std::make_shared<ChBody>();
  my_system.AddBody(RodDE);
  RodDE->SetPos(RodDECG);
  RodDE->SetRot(RodDERot);
  RodDE->SetMass(mass);
  RodDE->SetInertiaXX(inertiaXX);
  // Add some geometry for visualization
  auto cyl_RodDE = std::make_shared<ChCylinderShape>();
  cyl_RodDE->GetCylinderGeometry().p1 = RodDERot.Rotate(PntDLoc - RodDECG);
  cyl_RodDE->GetCylinderGeometry().p2 = RodDERot.Rotate(PntELoc - RodDECG);
  cyl_RodDE->GetCylinderGeometry().rad = RodGraphicRad;
  RodDE->AddAsset(cyl_RodDE);
  auto cyl_E = std::make_shared<ChCylinderShape>();
  cyl_E->GetCylinderGeometry().p1 = ChVector<>(0, 0, -.5*length) + RodDERot.Rotate(PntERot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_E->GetCylinderGeometry().p2 = ChVector<>(0, 0, -.5*length) - RodDERot.Rotate(PntERot.Rotate(ChVector<>(0, 0, RevGraphicLen)));
  cyl_E->GetCylinderGeometry().rad = RodGraphicRad;
  RodDE->AddAsset(cyl_E);

  // Create the RodEF body in an initial configuration
  auto RodEF = std::make_shared<ChBody>();
  my_system.AddBody(RodEF);
  RodEF->SetPos(RodEFCG);
  RodEF->SetRot(RodEFRot);
  RodEF->SetMass(mass);
  RodEF->SetInertiaXX(inertiaXX);
  // Add some geometry for visualization
  auto cyl_RodEF = std::make_shared<ChCylinderShape>();
  cyl_RodEF->GetCylinderGeometry().p1 = RodEFRot.Rotate(PntELoc - RodEFCG);
  cyl_RodEF->GetCylinderGeometry().p2 = RodEFRot.Rotate(PntFLoc - RodEFCG);
  cyl_RodEF->GetCylinderGeometry().rad = RodGraphicRad;
  RodEF->AddAsset(cyl_RodEF);

  //// Create revolute joints between the parts
  //// The revolute joint's axis of rotation will be the Z axis
  //// of the specified rotation matrix.

  auto revoluteJointA = std::make_shared<ChLinkLockRevolute>();
  revoluteJointA->Initialize(RodAB, ground, ChCoordsys<>(PntALoc, PntARot));
  my_system.AddLink(revoluteJointA);

  auto revoluteJointB = std::make_shared<ChLinkLockRevolute>();
  revoluteJointB->Initialize(RodBC, RodAB, ChCoordsys<>(PntBLoc, PntBRot));
  my_system.AddLink(revoluteJointB);

  auto revoluteJointC = std::make_shared<ChLinkLockRevolute>();
  revoluteJointC->Initialize(RodCD, RodBC, ChCoordsys<>(PntCLoc, PntCRot));
  my_system.AddLink(revoluteJointC);

  auto revoluteJointD = std::make_shared<ChLinkLockRevolute>();
  revoluteJointD->Initialize(RodDE, RodCD, ChCoordsys<>(PntDLoc, PntDRot));
  my_system.AddLink(revoluteJointD);

  auto revoluteJointE = std::make_shared<ChLinkLockRevolute>();
  revoluteJointE->Initialize(RodEF, RodDE, ChCoordsys<>(PntELoc, PntERot));
  my_system.AddLink(revoluteJointE);

  auto revoluteJointF = std::make_shared<ChLinkLockRevolute>();
  revoluteJointF->Initialize(RodEF, ground, ChCoordsys<>(PntFLoc, PntFRot));
  my_system.AddLink(revoluteJointF);



  // Perform the simulation (animation with Irrlicht option)
  // -------------------------------------------------------

  if (animate)
  {
      ChVector<> CameraPnt(0.5, -1.0, -0.5);
      // Create the Irrlicht application for visualization
      ChIrrApp * application = new ChIrrApp(&my_system, L"ME751 - HW9 - Problem 3", core::dimension2d<u32>(800, 600), false, true);
      application->AddTypicalLogo();
      application->AddTypicalSky();
      application->AddTypicalLights();
      core::vector3df lookat((f32)CameraPnt.x, (f32)CameraPnt.y, (f32)CameraPnt.z);
      application->AddTypicalCamera(lookat + core::vector3df(-2, 3, 4), lookat);

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
          mcamera->SetPosition(ChVector<>(-2, 3, 4));
          mcamera->SetAimPoint(ChVector<>(0.5, -1.0, 0.5));
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
  utils::CSV_writer out_pos_C = OutStream();


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
          ChVector<> &Pnt3 = RodBC->GetFrame_REF_to_abs().TransformPointLocalToParent(ChVector<>(0, 0, -.5*length));

          out_pos_C << simTime << Pnt3 << std::endl;

          // Increment output time
          outTime += outTimeStep;
      }

      // Advance simulation by one step
      my_system.DoStepDynamics(simTimeStep);

      // Increment simulation time
      simTime += simTimeStep;
  }

  // Write output files
  out_pos_C.write_to_file(out_dir + "HW9P3_CHRONO_Pnt3_Pos.txt", "");

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
