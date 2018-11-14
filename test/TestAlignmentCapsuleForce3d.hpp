
#ifndef _TESTALILGNMENTCAPSULEFORCE_HPP_
#define _TESTALIGNMENTCAPSULEFORCE_HPP_

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
#include "AttractiveEndsAlignmentCapsuleForce.hpp"
#include "AlignmentCapsuleForce.hpp"
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


class TestAlignmentCapsuleForce : public AbstractCellBasedTestSuite
{
public:
	
	void TestAlignmentCapsuleForceCases2DCase01() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at 0 Pi and Capsule 2 is at 1/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


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
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = .25*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase01");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestAlignmentCapsuleForceCases2DCase02() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at  Pi and Capsule 2 is at 1/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = .25*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase02");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}
	void TestAlignmentCapsuleForceCases2DCase03() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at 0 Pi and Capsule 2 is at 5/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


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
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 1.25*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase03");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestAlignmentCapsuleForceCases2DCase04() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at Pi and Capsule 2 is at 5/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 1.25*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase04");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}
	
	
	void TestAlignmentCapsuleForceCases2DCase05() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at 0 Pi and Capsule 2 is at 1/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(-1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


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
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = .75*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase05");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestAlignmentCapsuleForceCases2DCase06() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at  Pi and Capsule 2 is at 1/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(-1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = .75*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase06");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}
	void TestAlignmentCapsuleForceCases2DCase07() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at 0 Pi and Capsule 2 is at 5/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(-1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


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
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 1.75*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase07");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}

	void TestAlignmentCapsuleForceCases2DCase08() 
	{
		EXIT_IF_PARALLEL;

		//Case 01: Capsule 1 is at Pi and Capsule 2 is at 5/4 Pi
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(-1.5,1.5)));
		//nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,1.5)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 1.75*M_PI;
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
		simulator.SetOutputDirectory("TestAlignmentCapsuleForceCases2DCase08");
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

		auto p_alignment_capsule_force = boost::make_shared<AlignmentCapsuleForce<2>>();
		simulator.AddForce(p_alignment_capsule_force);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
	}
	void TestAttractiveEndsAlignmentDivision2dCase01()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0*M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


		//Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.01);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);


			p_cell->SetBirthTime(-0.9);
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);

			MAKE_PTR(TypeSixMachineProperty, p_property);
			p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_coordinates));

			p_cell->AddCellProperty(p_property);

			cells.push_back(p_cell);
		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
				 population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestAttractiveEndsAlignmentDivision2dCase01");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(30u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */

		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		simulator.AddForce(p_capsule_force);


		auto p_atractive_ends_capsule_force = boost::make_shared<AttractiveEndsCapsuleForce<2>>();
		p_atractive_ends_capsule_force->SetYoungModulus(100.0);
		simulator.AddForce(p_atractive_ends_capsule_force);

		auto p_alignment_capsule_force = boost::make_shared<AttractiveEndsAlignmentCapsuleForce<2>>();
		p_alignment_capsule_force->SetGamma(150.0);
		simulator.AddForce(p_alignment_capsule_force);



		/* We then set an end time and run the simulation */
		simulator.SetEndTime(6.0);
		simulator.Solve();
		//PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
	}

	
}; 

#endif /*_TESTALIGNMENTSCAPSULEFORCE_HPP_*/

