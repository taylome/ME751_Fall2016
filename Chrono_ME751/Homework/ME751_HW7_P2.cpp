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
// ME751- HW7 - Problem#2
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

#include "ChronoME751_config.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace irr;
using namespace postprocess;

// =============================================================================
// Local variables
//
static const std::string val_dir = "../RESULTS/";
static const std::string out_dir = val_dir + "HW7/";

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
  double outTimeStep = 1. / 60.;//5e-1;


  // Settings
  //---------

  // There are no units in Chrono, so values must be consistent
  // (MKS is used in this example)

  double mass = 1.0;                     // mass of link
  double length = 1.0;                   // length of link
  ChVector<> inertiaXX(0.0833333, 0.0001, 0.0833333);  // mass moments of inertia of pendulum (centroidal frame)
  double g = 9.81;

  double timeRecord = 10;                 // simulation length

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
  my_system.SetTol(1e-6);
  my_system.SetTolForce(1e-4);

  // Create the ground body
  ChVector<>joint1Loc(0, 0, 0);
  auto joint1Rot = QUNIT;
  ChVector<>joint2Loc(0, 1, 0);
  auto joint2Rot = QUNIT;
  ChVector<>joint3Loc(1, 1, 0);
  auto joint3Rot = QUNIT;
  ChVector<>joint4Loc(1, 0, 0);
  auto joint4Rot = QUNIT;

  auto ground = std::make_shared<ChBody>();
  my_system.AddBody(ground);
  ground->SetBodyFixed(true);
  // Add some geometry to the ground body for visualizing the first revolute joint
  auto cyl_g = std::make_shared<ChCylinderShape>();
  cyl_g->GetCylinderGeometry().p1 = joint1Loc + joint1Rot.Rotate(ChVector<>(0, 0, 0.1));
  cyl_g->GetCylinderGeometry().p2 = joint1Loc + joint1Rot.Rotate(ChVector<>(0, 0, -0.1));
  cyl_g->GetCylinderGeometry().rad = 0.05;
  ground->AddAsset(cyl_g);

  // Add some geometry to the ground body for visualizing the last revolute joint
  auto cyl_g4 = std::make_shared<ChCylinderShape>();
  cyl_g4->GetCylinderGeometry().p1 = joint4Loc + joint4Rot.Rotate(ChVector<>(0, 0, 0.1));
  cyl_g4->GetCylinderGeometry().p2 = joint4Loc + joint4Rot.Rotate(ChVector<>(0, 0, -0.1));
  cyl_g4->GetCylinderGeometry().rad = 0.05;
  ground->AddAsset(cyl_g4);


  // Create the first link body in an initial configuration
  ChCoordsys<>Link1CSYS(ChVector<>(0, 0.5, 0), QUNIT);
  
  auto link1 = std::make_shared<ChBody>();
  my_system.AddBody(link1);
  link1->SetPos(Link1CSYS.pos);
  link1->SetRot(Link1CSYS.rot);
  link1->SetMass(mass);
  link1->SetInertiaXX(inertiaXX);
  // Add some geometry to the pendulum for visualization
  auto box_link1 = std::make_shared<ChBoxShape>();
  box_link1->GetBoxGeometry().Size = ChVector<>(0.05 * length, 0.5 * length, 0.05 * length);
  link1->AddAsset(box_link1);

  link1->SetPos_dt(ChVector<>(0.5, 0, 0));

  // Create the second link body in an initial configuration
  ChCoordsys<>Link2CSYS(ChVector<>(0.5, 1, 0), Q_from_AngZ(CH_C_PI_2));

  auto link2 = std::make_shared<ChBody>();
  my_system.AddBody(link2);
  link2->SetPos(Link2CSYS.pos);
  link2->SetRot(Link2CSYS.rot);
  link2->SetMass(mass);
  link2->SetInertiaXX(inertiaXX);
  // Add some geometry to the pendulum for visualization
  auto box_link2 = std::make_shared<ChBoxShape>();
  box_link2->GetBoxGeometry().Size = ChVector<>(0.05 * length, 0.5 * length, 0.05 * length);
  link2->AddAsset(box_link2);
  // Add some geometry to the second body for visualizing the 2nd and 3rd revolute joints
  auto cyl_g2 = std::make_shared<ChCylinderShape>();
  cyl_g2->GetCylinderGeometry().p1 = ChVector<>(0, 0.5, 0) + joint2Rot.Rotate(ChVector<>(0, 0, 0.1));
  cyl_g2->GetCylinderGeometry().p2 = ChVector<>(0, 0.5, 0) + joint2Rot.Rotate(ChVector<>(0, 0, -0.1));
  cyl_g2->GetCylinderGeometry().rad = 0.05;
  link2->AddAsset(cyl_g2);
  auto cyl_g3 = std::make_shared<ChCylinderShape>();
  cyl_g3->GetCylinderGeometry().p1 = ChVector<>(0, -0.5, 0) + joint3Rot.Rotate(ChVector<>(0, 0, 0.1));
  cyl_g3->GetCylinderGeometry().p2 = ChVector<>(0, -0.5, 0) + joint3Rot.Rotate(ChVector<>(0, 0, -0.1));
  cyl_g3->GetCylinderGeometry().rad = 0.05;
  link2->AddAsset(cyl_g3);

  link2->SetPos_dt(ChVector<>(1, 0, 0));

  // Create the third link body in an initial configuration
  ChCoordsys<>Link3CSYS(ChVector<>(1, 0.5, 0), QUNIT);

  auto link3 = std::make_shared<ChBody>();
  my_system.AddBody(link3);
  link3->SetPos(Link3CSYS.pos);
  link3->SetRot(Link3CSYS.rot);
  link3->SetMass(mass);
  link3->SetInertiaXX(inertiaXX);
  // Add some geometry to the pendulum for visualization
  auto box_link3 = std::make_shared<ChBoxShape>();
  box_link3->GetBoxGeometry().Size = ChVector<>(0.05 * length, 0.5 * length, 0.05 * length);
  link3->AddAsset(box_link3);

  link3->SetPos_dt(ChVector<>(0.5, 0, 0));

  // Create revolute or spherical joints between the links and/or ground at "loc" in the global
  // reference frame. The revolute joint's axis of rotation will be the Z axis
  // of the specified rotation matrix.

  auto revoluteJoint1 = std::make_shared<ChLinkLockRevolute>();
  revoluteJoint1->Initialize(link1, ground, ChCoordsys<>(joint1Loc, joint1Rot));
  my_system.AddLink(revoluteJoint1);

  auto revoluteJoint2 = std::make_shared<ChLinkLockSpherical>();
  revoluteJoint2->Initialize(link2, link1, ChCoordsys<>(joint2Loc, joint2Rot));
  my_system.AddLink(revoluteJoint2);

  auto revoluteJoint3 = std::make_shared<ChLinkLockSpherical>();
  revoluteJoint3->Initialize(link3, link2, ChCoordsys<>(joint3Loc, joint3Rot));
  my_system.AddLink(revoluteJoint3);

  auto revoluteJoint4 = std::make_shared<ChLinkLockRevolute>();
  revoluteJoint4->Initialize(link3, ground, ChCoordsys<>(joint4Loc, joint4Rot));
  my_system.AddLink(revoluteJoint4);

  // Perform the simulation (animation with Irrlicht option)
  // -------------------------------------------------------

  if (animate)
  {
      ChVector<> CameraPnt(0.5,0,0);
      // Create the Irrlicht application for visualization
      ChIrrApp * application = new ChIrrApp(&my_system, L"ME751 - HW7 - Problem 2", core::dimension2d<u32>(800, 600), false, true);
      application->AddTypicalLogo();
      application->AddTypicalSky();
      application->AddTypicalLights();
      core::vector3df lookat((f32)CameraPnt.x, (f32)CameraPnt.y, (f32)CameraPnt.z);
      application->AddTypicalCamera(lookat + core::vector3df(0, 1.5, -3), lookat);

      // Now have the visulization tool (Irrlicht) create its geometry from the
      // assets defined above
      application->AssetBindAll();
      application->AssetUpdateAll();

      application->SetTimestep(simTimeStep);

      // Simulation loop
      double outTime = 0;
      int    outFrame = 1;
      double time = 0;

      std::string pov_dir = out_dir + "POVRAY_HW7P2";
      if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
          std::cout << "Error creating directory " << pov_dir << std::endl;
          return false;
      }


      // Create an exporter to POVray
      ChPovRay pov_exporter = ChPovRay(&my_system);

      if (save) {

          // ==Asset== Attach a video camera. This will be used by Irrlicht,
          // or POVray postprocessing, etc. Note that a camera can also be
          // put in a moving object.
          auto mcamera = std::make_shared<ChCamera>();
          mcamera->SetAngle(50);
          mcamera->SetPosition(ChVector<>(-3., 4., -5.));
          mcamera->SetAimPoint(ChVector<>(0, 0.5, 0));
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

          // Draw an XZ grid at the global origin to add in visualization
          ChIrrTools::drawGrid(
              application->GetVideoDriver(), 1, 1, 20, 20,
              ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(CH_C_PI_2)),
              video::SColor(255, 80, 100, 100), true);

          application->DoStep();  //Take one step in time
          application->EndScene();
      }

      return true;
  }

  // Perform the simulation (record results option)
  // ------------------------------------------------

  // Create the CSV_Writer output objects (TAB delimited)
  utils::CSV_writer out_pos_LeftLinkTip = OutStream();
  utils::CSV_writer out_vel_LeftLinkTip = OutStream();
  utils::CSV_writer out_rfrc_LeftLinkTip = OutStream();

  utils::CSV_writer out_pos_RightLinkTip = OutStream();
  utils::CSV_writer out_vel_RightLinkTip = OutStream();
  utils::CSV_writer out_rfrc_RightLinkTip = OutStream();


  // Write headers
  //out_pos_LeftLinkTip << "Time" << "X_Pos" << "Y_Pos" << "Z_Pos" << std::endl;
  //out_vel_LeftLinkTip << "Time" << "X_Vel" << "Y_Vel" << "Z_Vel" << std::endl;
  //out_rfrc_LeftLinkTip << "Time" << "X_Force" << "Y_Force" << "Z_Force" << std::endl;

  //out_pos_RightLinkTip << "Time" << "X_Pos" << "Y_Pos" << "Z_Pos" << std::endl;
  //out_vel_RightLinkTip << "Time" << "X_Vel" << "Y_Vel" << "Z_Vel" << std::endl;
  //out_rfrc_RightLinkTip << "Time" << "X_Force" << "Y_Force" << "Z_Force" << std::endl;


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

          // CM position, velocity, and acceleration (expressed in global frame).
          ChVector<> &link1pos = link1->GetFrame_REF_to_abs().TransformPointLocalToParent(ChVector<>(0, 0.5, 0));
          ChVector<> &link1vel = link1->GetFrame_REF_to_abs().PointSpeedLocalToParent(ChVector<>(0, 0.5, 0));
          ChVector<> &link3pos = link3->GetFrame_REF_to_abs().TransformPointLocalToParent(ChVector<>(0, 0.5, 0));
          ChVector<> &link3vel = link3->GetFrame_REF_to_abs().PointSpeedLocalToParent(ChVector<>(0, 0.5, 0));

          out_pos_LeftLinkTip << simTime << link1pos << std::endl;
          out_vel_LeftLinkTip << simTime << link1vel << std::endl;
          out_pos_RightLinkTip << simTime << link3pos << std::endl;
          out_vel_RightLinkTip << simTime << link3vel << std::endl;

          //    joint frames
          ChCoordsys<> link1Coordsys = revoluteJoint2->GetLinkRelativeCoords();
          ChCoordsys<> link3Coordsys = revoluteJoint3->GetLinkRelativeCoords();

          //    reaction force
          ChVector<> reactForce2 = revoluteJoint2->Get_react_force();
          ChVector<> reactForce3 = revoluteJoint3->Get_react_force();

          //    reaction force, expressed in ground frame
          reactForce2 = link1Coordsys.TransformDirectionLocalToParent(reactForce2);
          reactForce3 = link3Coordsys.TransformDirectionLocalToParent(reactForce3);

          out_rfrc_LeftLinkTip << simTime << reactForce2 << std::endl;
          out_rfrc_RightLinkTip << simTime << reactForce3 << std::endl;


          // Increment output time
          outTime += outTimeStep;
      }

      // Advance simulation by one step
      my_system.DoStepDynamics(simTimeStep);

      // Increment simulation time
      simTime += simTimeStep;
  }

  // Write output files
  out_pos_LeftLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link1Tip_Pos.txt", "");
  out_vel_LeftLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link1Tip_Vel.txt", "");
  out_rfrc_LeftLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link1Tip_RFrc.txt", "");

  out_pos_RightLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link3Tip_Pos.txt", "");
  out_vel_RightLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link3Tip_Vel.txt", "");
  out_rfrc_RightLinkTip.write_to_file(out_dir + "HW7P2_CHRONO_Link3Tip_RFrc.txt", "");

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
