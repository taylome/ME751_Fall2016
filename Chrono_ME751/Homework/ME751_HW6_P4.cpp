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
// ME751- HW6 - Problem#4
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
static const std::string out_dir = val_dir + "HW6/";

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

  double mass = 1.0;                     // mass of pendulum
  double length = 1.0;                   // length of pendulum
  ChVector<> inertiaXX(0.001, 0.001, 0.001);  // mass moments of inertia of pendulum (centroidal frame)
  double g = 9.81;

  double timeRecord = 5;                 // simulation length

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
  ChVector<>jointLocGnd(0, 0, 0);
  ChVector<>jointRevAxis(0, 0, 1);

  auto ground = std::make_shared<ChBody>();
  my_system.AddBody(ground);
  ground->SetBodyFixed(true);
  // Add some geometry to the ground body for visualizing the revolute joint
  auto cyl_g = std::make_shared<ChCylinderShape>();
  cyl_g->GetCylinderGeometry().p1 = jointLocGnd + jointRevAxis*0.4;
  cyl_g->GetCylinderGeometry().p2 = jointLocGnd + jointRevAxis*-0.4;
  cyl_g->GetCylinderGeometry().rad = 0.05;
  ground->AddAsset(cyl_g);

  // Create the pendulum body in an initial configuration at rest, with an 
  // orientatoin that matches the specified joint orientation and a position
  // consistent with the specified joint location.
  // The pendulum CG is assumed to be at half its length.

  ChCoordsys<>PendCSYS(ChVector<>(-1, 0, 0), QUNIT);
  ChVector<>jointLocPend(-1, 0, 0);
  
  auto pendulum = std::make_shared<ChBody>();
  my_system.AddBody(pendulum);
  pendulum->SetPos(PendCSYS.pos);
  pendulum->SetRot(PendCSYS.rot);
  pendulum->SetMass(mass);
  pendulum->SetInertiaXX(inertiaXX);
  // Add some geometry to the pendulum for visualization
  auto sph_p = std::make_shared<ChSphereShape>();
  sph_p->GetSphereGeometry().center = pendulum->TransformPointParentToLocal(jointLocPend);
  sph_p->GetSphereGeometry().rad = 0.05;
  pendulum->AddAsset(sph_p);
  //auto box_p = std::make_shared<ChBoxShape>();
  //box_p->GetBoxGeometry().Size = ChVector<>(0.05 * length, 0.5 * length, 0.05 * length);
  //pendulum->AddAsset(box_p);

  // Create a Revolute-Spherical constraint between pendulum at "jointLocPend" 
  // and ground at "jointLocGnd" in the global reference frame. 
  // The constrained distance is set equal to the inital distance between
  // "jointLocPend" and "jointLocGnd".

  auto revSphericalConstraint = std::make_shared<ChLinkRevoluteSpherical>();
  revSphericalConstraint->Initialize(ground, pendulum, false, jointLocGnd, jointRevAxis, jointLocPend, true);
  my_system.AddLink(revSphericalConstraint);

  // Perform the simulation (animation with Irrlicht option)
  // -------------------------------------------------------

  if (animate)
  {

      // Create the Irrlicht application for visualization
      ChIrrApp * application = new ChIrrApp(&my_system, L"ME751 - HW6 - Problem 4", core::dimension2d<u32>(800, 600), false, true);
      application->AddTypicalLogo();
      application->AddTypicalSky();
      application->AddTypicalLights();
      core::vector3df lookat((f32)jointLocGnd.x, (f32)jointLocGnd.y, (f32)jointLocGnd.z);
      application->AddTypicalCamera(lookat + core::vector3df(0, 3, -6), lookat);

      // Now have the visulization tool (Irrlicht) create its geometry from the
      // assets defined above
      application->AssetBindAll();
      application->AssetUpdateAll();

      application->SetTimestep(simTimeStep);

      // Simulation loop
      double outTime = 0;
      int    outFrame = 1;
      double time = 0;

      std::string pov_dir = out_dir + "POVRAY_HW6P4";
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
          mcamera->SetPosition(ChVector<>(-3. / 2., 4. / 2., -5. / 2.));
          mcamera->SetAimPoint(ChVector<>(0, -0.5, 0));
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

          // Render the rev-sph massless link.
          ChIrrTools::drawSegment(application->GetVideoDriver(),
              revSphericalConstraint->GetPoint1Abs(),
              revSphericalConstraint->GetPoint2Abs(),
              video::SColor(255, 0, 20, 0),
              true);

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
  utils::CSV_writer out_pos = OutStream();
  utils::CSV_writer out_vel = OutStream();
  utils::CSV_writer out_acc = OutStream();

  utils::CSV_writer out_quat = OutStream();
  utils::CSV_writer out_avel = OutStream();
  utils::CSV_writer out_aacc = OutStream();

  utils::CSV_writer out_rfrc = OutStream();
  utils::CSV_writer out_rtrq = OutStream();
  utils::CSV_writer out_rfrc1 = OutStream();
  utils::CSV_writer out_rtrq1 = OutStream();
  utils::CSV_writer out_rfrc2 = OutStream();
  utils::CSV_writer out_rtrq2 = OutStream();

  utils::CSV_writer out_energy = OutStream();

  utils::CSV_writer out_cnstr = OutStream();

  // Write headers
  out_pos << "Time" << "X_Pos" << "Y_Pos" << "Z_Pos" << std::endl;
  out_vel << "Time" << "X_Vel" << "Y_Vel" << "Z_Vel" << std::endl;
  out_acc << "Time" << "X_Acc" << "Y_Acc" << "Z_Acc" << std::endl;

  out_quat << "Time" << "e0" << "e1" << "e2" << "e3" << std::endl;
  out_avel << "Time" << "X_AngVel" << "Y_AngVel" << "Z_AngVel" << std::endl;
  out_aacc << "Time" << "X_AngAcc" << "Y_AngAcc" << "Z_AngAcc" << std::endl;

  out_rfrc << "Time" << "X_Force" << "Y_Force" << "Z_Force" << std::endl;
  out_rtrq << "Time" << "X_Torque" << "Y_Torque" << "Z_Torque" << std::endl;

  out_rfrc1 << "Time" << "X_Force" << "Y_Force" << "Z_Force" << std::endl;
  out_rtrq1 << "Time" << "X_Torque" << "Y_Torque" << "Z_Torque" << std::endl;

  out_rfrc2 << "Time" << "X_Force" << "Y_Force" << "Z_Force" << std::endl;
  out_rtrq2 << "Time" << "X_Torque" << "Y_Torque" << "Z_Torque" << std::endl;


  out_energy << "Time" << "Transl_KE" << "Rot_KE" << "Delta_PE" << "KE+PE" << std::endl;

  out_cnstr << "Time" << "Cnstr_1" << "Cnstr_2" << std::endl;

  // Perform a system assembly to ensure we have the correct accelerations at
  // the initial time.
  my_system.DoFullAssembly();

  // Total energy at initial time.
  ChMatrix33<> inertia = pendulum->GetInertia();
  ChVector<> angVelLoc = pendulum->GetWvel_loc();
  double transKE = 0.5 * mass * pendulum->GetPos_dt().Length2();
  double rotKE = 0.5 * Vdot(angVelLoc, inertia * angVelLoc);
  double deltaPE = mass * g * (pendulum->GetPos().z - PendCSYS.pos.z);
  double totalE0 = transKE + rotKE + deltaPE;

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
          const ChVector<>& position = pendulum->GetPos();
          const ChVector<>& velocity = pendulum->GetPos_dt();
          out_pos << simTime << position << std::endl;
          out_vel << simTime << velocity << std::endl;
          out_acc << simTime << pendulum->GetPos_dtdt() << std::endl;

          // Orientation, angular velocity, and angular acceleration (expressed in
          // global frame).
          out_quat << simTime << pendulum->GetRot() << std::endl;
          out_avel << simTime << pendulum->GetWvel_par() << std::endl;
          out_aacc << simTime << pendulum->GetWacc_par() << std::endl;


          // Chrono returns the reaction force and torque on body 2 (as specified in
          // the joint Initialize() function), as applied at the joint location and
          // expressed in the joint frame. Here, the 2nd body is the pendulum.

          //    joint frame on 2nd body (pendulum), expressed in the body frame
          ChCoordsys<> linkCoordsys = revSphericalConstraint->GetLinkRelativeCoords();

          //    reaction force and torque on pendulum, expressed in joint frame
          //       at the joint frame origin (center of the revolute)
          ChVector<> reactForce = revSphericalConstraint->Get_react_force();
          ChVector<> reactTorque = revSphericalConstraint->Get_react_torque();

          //    reaction force and torque on the ground, expressed in joint frame
          //       at the revolute joint center (joint frame origin)
          ChVector<> reactForceB1 = revSphericalConstraint->Get_react_force_body1();
          ChVector<> reactTorqueB1 = revSphericalConstraint->Get_react_torque_body1();

          //    reaction force and torque on the ground, expressed in joint frame
          //       at the spherical joint center
          ChVector<> reactForceB2 = revSphericalConstraint->Get_react_force_body2();
          ChVector<> reactTorqueB2 = revSphericalConstraint->Get_react_torque_body2();

          //    Transform from the joint frame into the pendulum frame
          reactForce = linkCoordsys.TransformDirectionLocalToParent(reactForce);
          reactTorque = linkCoordsys.TransformDirectionLocalToParent(reactTorque);
          reactForceB1 = linkCoordsys.TransformDirectionLocalToParent(reactForceB1);
          reactTorqueB1 = linkCoordsys.TransformDirectionLocalToParent(reactTorqueB1);
          reactForceB2 = linkCoordsys.TransformDirectionLocalToParent(reactForceB2);
          reactTorqueB2 = linkCoordsys.TransformDirectionLocalToParent(reactTorqueB2);

          //    Transform from the joint frame into the global frame
          reactForce = pendulum->TransformDirectionLocalToParent(reactForce);
          reactTorque = pendulum->TransformDirectionLocalToParent(reactTorque);
          reactForceB1 = pendulum->TransformDirectionLocalToParent(reactForceB1);
          reactTorqueB1 = pendulum->TransformDirectionLocalToParent(reactTorqueB1);
          reactForceB2 = pendulum->TransformDirectionLocalToParent(reactForceB2);
          reactTorqueB2 = pendulum->TransformDirectionLocalToParent(reactTorqueB2);

          out_rfrc << simTime << reactForce << std::endl;
          out_rtrq << simTime << reactTorque << std::endl;
          out_rfrc1 << simTime << reactForceB1 << std::endl;
          out_rtrq1 << simTime << reactTorqueB1 << std::endl;
          out_rfrc2 << simTime << reactForceB2 << std::endl;
          out_rtrq2 << simTime << reactTorqueB2 << std::endl;

          // Conservation of Energy
          // Translational Kinetic Energy (1/2*m*||v||^2)
          // Rotational Kinetic Energy (1/2 w'*I*w)
          // Delta Potential Energy (m*g*dz)
          ChMatrix33<> inertia = pendulum->GetInertia();
          ChVector<> angVelLoc = pendulum->GetWvel_loc();
          double transKE = 0.5 * mass * velocity.Length2();
          double rotKE = 0.5 * Vdot(angVelLoc, inertia * angVelLoc);
          double deltaPE = mass * g * (position.z - PendCSYS.pos.z);
          double totalE = transKE + rotKE + deltaPE;
          out_energy << simTime << transKE << rotKE << deltaPE << totalE - totalE0 << std::endl;;

          // Constraint violations
          ChMatrix<>* C = revSphericalConstraint->GetC();
          out_cnstr << simTime
              << C->GetElement(0, 0)
              << C->GetElement(1, 0) << std::endl;


          // Increment output time
          outTime += outTimeStep;
      }

      // Advance simulation by one step
      my_system.DoStepDynamics(simTimeStep);

      // Increment simulation time
      simTime += simTimeStep;
  }

  // Write output files
  out_pos.write_to_file(out_dir + "HW6P4_CHRONO_Pos.txt", "HW6P4\n\n");
  out_vel.write_to_file(out_dir + "HW6P4_CHRONO_Vel.txt", "HW6P4\n\n");
  out_acc.write_to_file(out_dir + "HW6P4_CHRONO_Acc.txt", "HW6P4\n\n");

  out_quat.write_to_file(out_dir + "HW6P4_CHRONO_Quat.txt", "HW6P4\n\n");
  out_avel.write_to_file(out_dir + "HW6P4_CHRONO_Avel.txt", "HW6P4\n\n");
  out_aacc.write_to_file(out_dir + "HW6P4_CHRONO_Aacc.txt", "HW6P4\n\n");

  out_rfrc.write_to_file(out_dir + "HW6P4_CHRONO_Rforce.txt", "HW6P4\n\n");
  out_rtrq.write_to_file(out_dir + "HW6P4_CHRONO_Rtorque.txt", "HW6P4\n\n");
  out_rfrc1.write_to_file(out_dir + "HW6P4_CHRONO_Rforce_Body1.txt", "HW6P4\n\n");
  out_rtrq1.write_to_file(out_dir + "HW6P4_CHRONO_Rtorque_Body1.txt", "HW6P4\n\n");
  out_rfrc2.write_to_file(out_dir + "HW6P4_CHRONO_Rforce_Body2.txt", "HW6P4\n\n");
  out_rtrq2.write_to_file(out_dir + "HW6P4_CHRONO_Rtorque_Body2.txt", "HW6P4\n\n");

  out_energy.write_to_file(out_dir + "HW6P4_CHRONO_Energy.txt", "HW6P4\n\n");

  out_cnstr.write_to_file(out_dir + "HW6P4_CHRONO_Constraints.txt", "HW6P4\n\n");

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
