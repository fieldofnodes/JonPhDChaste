
#ifndef TESTCAPSULESIMULATION2DATTRACTIVEENDS_HPP_
#define TESTCAPSULESIMULATION2DATTRACTIVEENDS_HPP_

#include <cxxtest/TestSuite.h>
#include <cycle/UniformCellCycleModel.hpp>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Cell.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellIdWriter.hpp"

// Header files included in this project
#include "TypeSixSecretionEnumerations.hpp"
#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "CapsuleForce.hpp"
#include "AttractiveEndsCapsuleForce.hpp"
#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "SquareBoundaryCondition.hpp"
#include "CapsuleBasedDivisionRule.hpp"
#include "TypeSixMachineModifier.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "MachineStateCountWriter.hpp"



// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsuleSimulation2dAttractiveEnds : public AbstractCellBasedTestSuite
{
public:

	void TestCapsule2dAttractiveEndsHorizontalOverlap() //throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;


			nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
			nodes.push_back(new Node<2>(1u, Create_c_vector(1.5, 0.6)));


			/*
			 * We then convert this list of nodes to a `NodesOnlyMesh`,
			 * which doesn't do very much apart from keep track of the nodes.
			 */
			NodesOnlyMesh<2> mesh;
			mesh.ConstructNodesWithoutMesh(nodes, 150.5);

			mesh.GetNode(0u)->AddNodeAttribute(0.0);
			mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

			mesh.GetNode(1u)->AddNodeAttribute(0.0);
			mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.0;
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


			//Create cells
			std::vector<CellPtr> cells;
			auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
			CellsGenerator<NoCellCycleModel, 2> cells_generator;
			cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

			// Create cell population
			NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

			population.AddCellWriter<CellIdWriter>();
			population.AddCellWriter<CapsuleOrientationWriter>();
			population.AddCellWriter<CapsuleScalingWriter>();

			// Create simulation
			OffLatticeSimulation<2> simulator(population);
			simulator.SetOutputDirectory("TestCapsule2dAttractiveEndsHorizontalOverlap");
			simulator.SetDt(1.0/1200.0);
			simulator.SetSamplingTimestepMultiple(1u);

			auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
			simulator.SetNumericalMethod(p_numerical_method);

			/*
			 * We now create a force law and pass it to the simulation
			 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
			 */

			auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
			simulator.AddForce(p_capsule_force);

			auto p_atractive_ends_capsule_force = boost::make_shared<AttractiveEndsCapsuleForce<2>>();
			simulator.AddForce(p_atractive_ends_capsule_force);

			/* We then set an end time and run the simulation */
			simulator.SetEndTime(100.0/1200.0);
			simulator.Solve();
	}
	void TestCapsule2dAttractiveEndsLSectionOverlap() //throw (Exception)
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;


			nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
			nodes.push_back(new Node<2>(1u, Create_c_vector(1.5, 1.5)));


			/*
			 * We then convert this list of nodes to a `NodesOnlyMesh`,
			 * which doesn't do very much apart from keep track of the nodes.
			 */
			NodesOnlyMesh<2> mesh;
			mesh.ConstructNodesWithoutMesh(nodes, 150.5);

			mesh.GetNode(0u)->AddNodeAttribute(0.0);
			mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

			mesh.GetNode(1u)->AddNodeAttribute(0.0);
			mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


			//Create cells
			std::vector<CellPtr> cells;
			auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
			CellsGenerator<NoCellCycleModel, 2> cells_generator;
			cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

			// Create cell population
			NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

			population.AddCellWriter<CellIdWriter>();
			population.AddCellWriter<CapsuleOrientationWriter>();
			population.AddCellWriter<CapsuleScalingWriter>();

			// Create simulation
			OffLatticeSimulation<2> simulator(population);
			simulator.SetOutputDirectory("TestCapsule2dAttractiveEndsLSectionOverlap");
			simulator.SetDt(1.0/1200.0);
			simulator.SetSamplingTimestepMultiple(1u);

			auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
			simulator.SetNumericalMethod(p_numerical_method);

			/*
			 * We now create a force law and pass it to the simulation
			 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
			 */

			auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
			simulator.AddForce(p_capsule_force);

			auto p_atractive_ends_capsule_force = boost::make_shared<AttractiveEndsCapsuleForce<2>>();
			simulator.AddForce(p_atractive_ends_capsule_force);

			/* We then set an end time and run the simulation */
			simulator.SetEndTime(100.0/1200.0);
			simulator.Solve();
	}
};
#endif /*TESTCAPSULESIMULATION2DATTRACTIVEENDS_HPP_*/
