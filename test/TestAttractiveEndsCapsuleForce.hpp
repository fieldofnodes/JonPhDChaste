
#ifndef _TESTATTRACTIVEENDSCAPSULEFORCE_HPP_
#define _TESTATTRACTIVEENDSCAPSULEFORCE_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "CapsuleForce.hpp"
#include "AttractiveEndsCapsuleForce.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "Node.hpp"
#include "NodeAttributes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestAttractiveEndsCapsuleForce : public AbstractCellBasedTestSuite
{
public:

    void TestDistanceBetweenTwoCapsules2d() // throw(Exception)
    {
        AttractiveEndsCapsuleForce<2, 2> force;

        // Case 1: Two horizontal rods on same y-line at y=3 and total capsule length = 3. Capsules are end to end.
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{3.0, 3.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{7.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            c_vector<double, 2u> vec_a_to_b;
            double distance_centre_mass_to_contact_a;
            double distance_centre_mass_to_contact_b;

            double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

            TS_ASSERT_DELTA(vec_a_to_b[0], 1.0, 1e-9);
            TS_ASSERT_DELTA(vec_a_to_b[1], 0.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 1.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -1.0, 1e-9);

            TS_ASSERT_DELTA(overlap, -1.0 , 1e-6);

            // swap around nodes and check it all still works

            overlap = force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

			TS_ASSERT_DELTA(vec_a_to_b[0],-1.0, 1e-9);
			TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
			TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -1.0, 1e-9);
			TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 1.0, 1e-9);

			TS_ASSERT_DELTA(overlap, -1.0, 1e-6);
        }
        // Case 2: One horizontal and one vertical capsule. With vertical capsule being like a "T" shape
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{3.0, 3.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{5.0, 5.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5*M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            c_vector<double, 2u> vec_a_to_b;
            double distance_centre_mass_to_contact_a;
            double distance_centre_mass_to_contact_b;

            double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

            TS_ASSERT_DELTA(vec_a_to_b[0], 1/sqrt(2.0), 1e-9);
            TS_ASSERT_DELTA(vec_a_to_b[1], 1/sqrt(2.0), 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 1.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -1.0, 1e-9);

            TS_ASSERT_DELTA(overlap, 1.0-sqrt(2.0) , 1e-6);

            // swap around nodes and check it all still works

            overlap = force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

			TS_ASSERT_DELTA(vec_a_to_b[0],-1/sqrt(2.0), 1e-9);
			TS_ASSERT_DELTA(vec_a_to_b[1],-1/sqrt(2.0), 1e-9);
			TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -1.0, 1e-9);
			TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 1.0, 1e-9);

			TS_ASSERT_DELTA(overlap, 1.0-sqrt(2.0), 1e-6);
        }
        // Case 3: Two parallel capsules vertical
        {
              // index, {x, y}
              Node<2> node_a(0u, std::vector<double>{3.0, 3.0});
              node_a.AddNodeAttribute(0.0);

              std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
              attributes_a.resize(NA_VEC_LENGTH);
              attributes_a[NA_THETA] = 0.5*M_PI;
              attributes_a[NA_LENGTH] = 2.0;
              attributes_a[NA_RADIUS] = 0.5;

              Node<2> node_b(0u, std::vector<double>{4.0, 3.0});
              node_b.AddNodeAttribute(0.0);

              std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
              attributes_b.resize(NA_VEC_LENGTH);
              attributes_b[NA_THETA] = 0.5*M_PI;
              attributes_b[NA_LENGTH] = 2.0;
              attributes_b[NA_RADIUS] = 0.5;

              c_vector<double, 2u> vec_a_to_b;
              double distance_centre_mass_to_contact_a;
              double distance_centre_mass_to_contact_b;

              double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

              TS_ASSERT_DELTA(vec_a_to_b[0], 1.0, 1e-9);
              TS_ASSERT_DELTA(vec_a_to_b[1], 0.0, 1e-9);
              TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
              TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

              TS_ASSERT_DELTA(overlap, 0.0 , 1e-6);

              // swap around nodes and check it all still works

              overlap = force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

  			TS_ASSERT_DELTA(vec_a_to_b[0],-1.0, 1e-9);
  			TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
  			TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
  			TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

  			TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
          }
        // Case 4: Intersecting cross shape, distance 0 from each other
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{3.0, 3.0});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = 0.0;
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{3.0, 3.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.5 * M_PI;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            c_vector<double, 2u> vec_a_to_b;
            double distance_centre_mass_to_contact_a;
            double distance_centre_mass_to_contact_b;

            double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

            TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
            TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

            TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
        }
        // Case 5: 3,4,5 Triangle
        {
            // index, {x, y}
            Node<2> node_a(0u, std::vector<double>{2.6, 2.8});
            node_a.AddNodeAttribute(0.0);

            std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
            attributes_a.resize(NA_VEC_LENGTH);
            attributes_a[NA_THETA] = atan(4.0/3.0);
            attributes_a[NA_LENGTH] = 2.0;
            attributes_a[NA_RADIUS] = 0.5;

            Node<2> node_b(0u, std::vector<double>{3.0, 1.0});
            node_b.AddNodeAttribute(0.0);

            std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
            attributes_b.resize(NA_VEC_LENGTH);
            attributes_b[NA_THETA] = 0.0;
            attributes_b[NA_LENGTH] = 2.0;
            attributes_b[NA_RADIUS] = 0.5;

            c_vector<double, 2u> vec_a_to_b;
            double distance_centre_mass_to_contact_a;
            double distance_centre_mass_to_contact_b;

            double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

            TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
            TS_ASSERT_DELTA(vec_a_to_b[1],-1.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -1.0, 1e-9);
            TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -1.0, 1e-9);

            TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
        }
    }
        void TestDistanceBetweenTwoCapsules3d() // throw(Exception)
        {
            AttractiveEndsCapsuleForce<3, 3> force;

            // Case 1: Two horizontal rods a distance 2 from each other about the x-axis
            {
                // index, {x, y, z}
                Node<3> node_a(0u, std::vector<double>{1.0, 0.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_PHI] = 0.0;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                Node<3> node_b(0u, std::vector<double>{1.0, 2.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_PHI] = 0.0;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],1.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }

            // Case 2: two vertical end to end about the z-axis
            {
                // index, {x, y, z}: this capsule is tc for shortest distance
                Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 1.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_PHI] =-M_PI;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                // index v {x,y,z}: this capsule is sc for shortest distance
                Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 4.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_PHI] = -M_PI;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],1.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -attributes_b[NA_LENGTH]/2.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, attributes_b[NA_LENGTH]/2.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }
            // Case 3: One capsule in the x-y plane parallel to the x-axis, one capsule in the x-y plane parallel to y-axis
            // in the shape of a "T"
            {
                // index, {x, y, z}: this capsule is tc for shortest distance
                Node<3> node_a(0u, std::vector<double>{3.0, 1.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_PHI] =0.5*M_PI;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                // index v {x,y,z}: this capsule is sc for shortest distance
                Node<3> node_b(0u, std::vector<double>{3.0, 3.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.5*M_PI;
                attributes_b[NA_PHI] = 0.5*M_PI;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],1.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -attributes_b[NA_LENGTH]/2.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }
            // Case 4: One capsule in the y-z plane parallel to the y-axis, one capsule in the y-z plane parallel to z-axis
            // in the shape of a "L"
            {
                // index, {x, y, z}: this capsule is tc for shortest distance
                Node<3> node_a(0u, std::vector<double>{0.0, 3.0, 1.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.5*M_PI;
                attributes_a[NA_PHI] = 0.5*M_PI;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                // index v {x,y,z}: this capsule is sc for shortest distance
                Node<3> node_b(0u, std::vector<double>{0.0, 2.0, 3.5});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.5*M_PI;
                attributes_b[NA_PHI] = 0.0;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],1.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -attributes_b[NA_LENGTH]/2.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -attributes_b[NA_LENGTH]/2.0, 1e-9);

                TS_ASSERT_DELTA(overlap, -0.5, 1e-6);
            }
            // Case 6: Intersecting cross shape, distance 0 from each other
            {
                // index, {x, y, z}
                Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_PHI] = M_PI/2.0;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.5 * M_PI;
                attributes_b[NA_PHI] = M_PI/2.0;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }
            // Case 7:   Capsule A is on the y-z plane parallel to the y-axis and
            // Capsule B is on the y-z plane at a PI/4 angle to the z-axis
            {
                // index, {x, y, z}
                Node<3> node_a(0u, std::vector<double>{0.0, 4.0, 1.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = M_PI/2.0;
                attributes_a[NA_PHI] = M_PI/2.0;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                Node<3> node_b(0u, std::vector<double>{0.0, 6.0+sqrt(2.0)/2.0, 2.0+sqrt(2.0)/2.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = M_PI/2.0;
                attributes_b[NA_PHI] = M_PI/4.0;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],sqrt(2.0)/2.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],sqrt(2.0)/2.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, attributes_b[NA_LENGTH]/2.0 , 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -attributes_b[NA_LENGTH]/2.0, 1e-9);

                TS_ASSERT_DELTA(overlap, -1.0*(sqrt(2.0)-1.0), 1e-6);
            }

            // Case 8: Crossing over each other one hovering 1 unit above in z axis.
            {
                // index, {x, y, z}
                Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 1.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_PHI] = M_PI/2.0;
                attributes_a[NA_LENGTH] = 2.0;
                attributes_a[NA_RADIUS] = 0.5;

                Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.5 * M_PI;
                attributes_b[NA_PHI] = M_PI/2.0;
                attributes_b[NA_LENGTH] = 2.0;
                attributes_b[NA_RADIUS] = 0.5;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],-1.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, 0.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, 0.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }
            // Case 9:
            {
                // index, {x, y, z}
                Node<3> node_a(0u, std::vector<double>{0.0, 0.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 1.0;
                attributes_a[NA_PHI] = M_PI/2.0;
                attributes_a[NA_LENGTH] = 4.0;
                attributes_a[NA_RADIUS] = 0.12;

                Node<3> node_b(0u, std::vector<double>{0.0, 0.0, 3.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 1.0;
                attributes_b[NA_PHI] = M_PI/2.0-atan(0.75);
                attributes_b[NA_LENGTH] = 4.0;
                attributes_b[NA_RADIUS] = 0.13;

                c_vector<double, 3u> vec_a_to_b;
                double distance_centre_mass_to_contact_a;
                double distance_centre_mass_to_contact_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, distance_centre_mass_to_contact_a, distance_centre_mass_to_contact_b);

                TS_ASSERT_DELTA(vec_a_to_b[0],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[1],0.0, 1e-9);
                TS_ASSERT_DELTA(vec_a_to_b[2],1.0, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_a, -1.6, 1e-9);
                TS_ASSERT_DELTA(distance_centre_mass_to_contact_b, -attributes_b[NA_LENGTH]/2.0, 1e-9);

                TS_ASSERT_DELTA(overlap, 0.0, 1e-6);
            }

        } // This bracket will  close the current void function: void TestDistanceBetweenTwoCapsules3d()
        void TestCalculateForceDirectionAndContactPoints() // throw(Exception)
        {
            AttractiveEndsCapsuleForce<2, 2> force;

            // Two horizontal rods a distance 2 from each other
            {
                // index, {x, y}
                Node<2> node_a(0u, std::vector<double>{10.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_LENGTH] = 2.0;

                Node<2> node_b(0u, std::vector<double>{0.0, 3.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_LENGTH] = 2.0;

                double contact_dist_a;
                double contact_dist_b;
                c_vector<double, 2> vec_a_to_b;

                double overlap = force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, contact_dist_a, contact_dist_b);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], -8.0 / sqrt(73.0), 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], 3.0 / sqrt(73.0), 1e-6);
                TS_ASSERT_DELTA(overlap, -sqrt(73), 1e-9);

                // Swap a->b for coverage
                force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, contact_dist_b, contact_dist_a);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 8.0 / sqrt(73.0), 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], -3.0 / sqrt(73.0), 1e-6);
            }

            // Two horizontal rods a distance 2 from each other
            {
                // index, {x, y}
                Node<2> node_a(0u, std::vector<double>{10.0, 0.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = 0.0;
                attributes_a[NA_LENGTH] = 2.0;

                Node<2> node_b(0u, std::vector<double>{0.0, 3.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_LENGTH] = 18.0;

                double contact_dist_a;
                double contact_dist_b;
                c_vector<double, 2> vec_a_to_b;

                force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, contact_dist_a, contact_dist_b);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);

                // Swap a->b
                force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, contact_dist_b, contact_dist_a);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, 0.5 * attributes_b[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);
            }

            // 3-4-5 triangle-based orientation
            {
                // index, {x, y}
                Node<2> node_a(0u, std::vector<double>{4.0, 3.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = atan(0.75);
                attributes_a[NA_LENGTH] = 4.0;

                Node<2> node_b(0u, std::vector<double>{4.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_LENGTH] = 4.0;

                double contact_dist_a;
                double contact_dist_b;
                c_vector<double, 2> vec_a_to_b;

                force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, contact_dist_a, contact_dist_b);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);

                // Swap a->b for coverage
                force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, contact_dist_b, contact_dist_a);

                TS_ASSERT_DELTA(contact_dist_a, -0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);
            }

            // 3-4-5 triangle-based orientation, with the horizontal capsule rotated 180 from the test above
            {
                // index, {x, y}
                Node<2> node_a(0u, std::vector<double>{4.0, 3.0});
                node_a.AddNodeAttribute(0.0);

                std::vector<double> &attributes_a = node_a.rGetNodeAttributes();
                attributes_a.resize(NA_VEC_LENGTH);
                attributes_a[NA_THETA] = atan(0.75) + M_PI;
                attributes_a[NA_LENGTH] = 4.0;

                Node<2> node_b(0u, std::vector<double>{4.0, 0.0});
                node_b.AddNodeAttribute(0.0);

                std::vector<double> &attributes_b = node_b.rGetNodeAttributes();
                attributes_b.resize(NA_VEC_LENGTH);
                attributes_b[NA_THETA] = 0.0;
                attributes_b[NA_LENGTH] = 4.0;

                double contact_dist_a;
                double contact_dist_b;
                c_vector<double, 2> vec_a_to_b;

                force.CalculateForceDirectionAndContactPoints(node_a, node_b, vec_a_to_b, contact_dist_a, contact_dist_b);

                TS_ASSERT_DELTA(contact_dist_a, 0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], -1.0, 1e-6);

                // Swap a->b for coverage
                force.CalculateForceDirectionAndContactPoints(node_b, node_a, vec_a_to_b, contact_dist_b, contact_dist_a);

                TS_ASSERT_DELTA(contact_dist_a, 0.5 * attributes_a[NA_LENGTH], 1e-6);
                TS_ASSERT_DELTA(contact_dist_b, -8.0 / 5.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(vec_a_to_b[1], 1.0, 1e-6);
            }
        }
        void TestCalculateForceMagnitude() // throw(Exception)
        {
            AttractiveEndsCapsuleForce<2, 2> force;

            // Hand calculate a trivial case
            {
                double overlap = 1.0;
                double radiusA = 4.0;
                double radiusB = 4.0;

                TS_ASSERT_DELTA(force.CalculateForceMagnitude(overlap, radiusA, radiusB), -800.0 / 3.0, 1e-6);
            }

            // Hand calculate a nontrivial case
            {
                double overlap = 1.23;
                double radiusA = 2.34;
                double radiusB = 3.45;

                TS_ASSERT_DELTA(force.CalculateForceMagnitude(overlap, radiusA, radiusB), -303.731332875, 1e-6);
            }
        }
        void TestAddForceContribution()
        {
            // Create two nodes
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u,  false,  0.0, 0.0));
            nodes.push_back(new Node<2>(1u,  false,  2.0, 0.0));

            // Create mesh with massive interaction distance so all nodes interact with each other
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1e6);

            mesh.GetNode(0u)->AddNodeAttribute(0.0);
            mesh.GetNode(0u)->ClearAppliedForce();
            mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
            mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
            mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
            mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.25;
            mesh.GetNode(0u)->rGetNodeAttributes()[NA_APPLIED_THETA] = 0.0;

            mesh.GetNode(1u)->AddNodeAttribute(0.0);
            mesh.GetNode(1u)->ClearAppliedForce();
            mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
            mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.5 * M_PI;
            mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
            mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.25;
            mesh.GetNode(1u)->rGetNodeAttributes()[NA_APPLIED_THETA] = 0.0;

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

            // Create cell population
            NodeBasedCellPopulation<2> population(mesh, cells);


            AttractiveEndsCapsuleForce<2, 2> force;

            force.AddForceContribution(population);


            // Nodes 0 and 1 are too far apart to interact with each other, so no applied force or angle
            {
                TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], 0.0, 1e-6);
                TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetNodeAttributes()[NA_APPLIED_THETA], 0.0, 1e-6);

                TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], 0.0, 1e-6);
                TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.0, 1e-6);
                TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetNodeAttributes()[NA_APPLIED_THETA], 0.0, 1e-6);
            }
        }

}; // This bracket will close the class: TestAttractiveEndsCapsuleForce : public AbstractCellBasedTestSuite

#endif /*_TESTATTRACTIVEENDSCAPSULEFORCE_HPP_*/
