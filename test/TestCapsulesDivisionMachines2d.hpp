
#ifndef TESTCAPSULESDIVISIONMACHINES2D_HPP_
#define TESTCAPSULESDIVISIONMACHINES2D_HPP_

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
#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "SquareBoundaryCondition.hpp"
#include "CapsuleBasedDivisionRule.hpp"
#include "TypeSixMachineModifier.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "NodeBasedCellPopulationCapsules.hpp"
#include "NodeBasedCellPopulationCapsulesMachines.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "TypeSixMachineCellLabelledKiller.hpp"
#include "MachineStateCountWriter.hpp"
#include "CapsuleTypeLabelWriter.hpp"




// Should usually be called last.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestCapsulesDivisionMachines2d : public AbstractCellBasedTestSuite
{
public:
//__________________________________________________________________________________________________________
// Start of Test 01
//__________________________________________________________________________________________________________

	void NoTestSingleCapsuleDivision()
            		   {
		EXIT_IF_PARALLEL;

		const unsigned num_nodes = 1u;
		//auto p_rand_gen = RandomNumberGenerator::Instance();

		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
		for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
		{
			c_vector<double, 2> safe_location;

			bool safe = false;
			while(!safe)
			{
				safe = true;
				safe_location = Create_c_vector(3.0, 3.0);

				for(auto&& p_node : nodes)
				{
					if(norm_2(p_node->rGetLocation() - safe_location) < 2.0)
					{
						safe = false;
					}
				}
			}

			nodes.push_back(new Node<2>(node_idx, safe_location));
		}

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		 NodesOnlyMesh<2> mesh;
		 mesh.ConstructNodesWithoutMesh(nodes, 100.0);
		 c_vector<double, 4> domain_size;
		 domain_size[0] = -1000.0;
		 domain_size[1] = 1000.0;
		 domain_size[2] = -1000.0;
		 domain_size[3] = 1000.0;
		 mesh.SetInitialBoxCollection(domain_size, 10.0);

		 mesh.GetNode(0u)->AddNodeAttribute(0.0);
		 mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		 for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
		 {
			 mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
			 mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_THETA] = 0.0;
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
		 }

		 // Create cells
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


			 double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			 std::vector<double> machine_angles;
			 machine_angles.push_back(rand_angle);
			 MAKE_PTR(TypeSixMachineProperty, p_property);
			 p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
			 p_cell->AddCellProperty(p_property);

			 //double birth_time = -RandomNumberGenerator::Instance()->ranf();
			 p_cell->SetBirthTime(-0.1);
			 mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			 cells.push_back(p_cell);
		 }



		 // Create cell population
		 NodeBasedCellPopulationCapsules<2> population(mesh, cells);


		 // Create writers to visualise
		 population.AddCellWriter<CellIdWriter>();
		 population.AddCellWriter<CapsuleOrientationWriter>();
		 population.AddCellWriter<CapsuleScalingWriter>();


		 // Create division rules
		 boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		 population.SetCentreBasedDivisionRule(p_division_rule);

		 // Create simulation
		 OffLatticeSimulation<2> simulator(population);
		 simulator.SetOutputDirectory("TestSingleCapsuleDivision");
		 double dt = 1.0/1200.0;
		 simulator.SetDt(dt);
		 simulator.SetSamplingTimestepMultiple(30u);

		 auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		 simulator.SetNumericalMethod(p_numerical_method);

		 //Add capsule force law.
		 auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		 p_capsule_force->SetYoungModulus(200.0);
		 simulator.AddForce(p_capsule_force);
		 
		 MARK;
		 /* We then set an end time and run the simulation */
		 simulator.SetEndTime(8.0);
		 simulator.Solve();
		 MARK;
		 //Tests
		 TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 2u);
     }
//__________________________________________________________________________________________________________
// End of Test 01
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
// Start of Test 02
//__________________________________________________________________________________________________________
	void NoTestSingleCapsuleDivisionMachines()
            		   {
		EXIT_IF_PARALLEL;

		const unsigned num_nodes = 1u;
		//auto p_rand_gen = RandomNumberGenerator::Instance();

		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));
		for (unsigned node_idx = 1; node_idx < num_nodes; ++node_idx)
		{
			c_vector<double, 2> safe_location;

			bool safe = false;
			while(!safe)
			{
				safe = true;
				safe_location = Create_c_vector(3.0, 3.0);

				for(auto&& p_node : nodes)
				{
					if(norm_2(p_node->rGetLocation() - safe_location) < 2.0)
					{
						safe = false;
					}
				}
			}

			nodes.push_back(new Node<2>(node_idx, safe_location));
		}

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		 NodesOnlyMesh<2> mesh;
		 mesh.ConstructNodesWithoutMesh(nodes, 100.0);
		 c_vector<double, 4> domain_size;
		 domain_size[0] = -1000.0;
		 domain_size[1] = 1000.0;
		 domain_size[2] = -1000.0;
		 domain_size[3] = 1000.0;
		 mesh.SetInitialBoxCollection(domain_size, 10.0);

		 mesh.GetNode(0u)->AddNodeAttribute(0.0);
		 mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.0;
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		 mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		 for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
		 {
			 mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
			 mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_THETA] = 0.0;
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			 mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
		 }

		 // Create cells
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


			 double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			 std::vector<double> machine_angles;
			 machine_angles.push_back(rand_angle);
			 MAKE_PTR(TypeSixMachineProperty, p_property);
			 p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
			 p_cell->AddCellProperty(p_property);

			 //double birth_time = -RandomNumberGenerator::Instance()->ranf();
			 p_cell->SetBirthTime(-0.1);
			 mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			 cells.push_back(p_cell);
		 }



		 // Create cell population
		 NodeBasedCellPopulationCapsulesMachines<2> population(mesh, cells);


		 // Create writers to visualise
		 population.AddCellWriter<CellIdWriter>();
		 population.AddCellWriter<CapsuleOrientationWriter>();
		 population.AddCellWriter<CapsuleScalingWriter>();
		 population.AddCellWriter<MachineStateCountWriter>();



		 // Create division rules
		 boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		 population.SetCentreBasedDivisionRule(p_division_rule);

		 // Create simulation
		 OffLatticeSimulation<2> simulator(population);
		 simulator.SetOutputDirectory("TestSingleCapsuleDivisionMachines");
		 double dt = 1.0/1200.0;
		 simulator.SetDt(dt);
		 simulator.SetSamplingTimestepMultiple(10u);

		 auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		 simulator.SetNumericalMethod(p_numerical_method);

		 //Add capsule force law.
		 auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		 p_capsule_force->SetYoungModulus(200.0);
		 simulator.AddForce(p_capsule_force);

		 // Add capsule modifier with machines
		 MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		 p_modifier->SetOutputDirectory("TestSingleCapsuleDivisionMachines");
		 //Set the variable for the states
		 //p_modifier->Setk_1(0.0);
		 //p_modifier->Setk_2(0.02);
		 //p_modifier->Setk_3(0.0);
		 //p_modifier->Setk_4(0.0);
		 //p_modifier->Setk_5(0.0);
		 //p_modifier->Setk_6(0.0);
		 //p_modifier->Setk_7(0.0);
		 p_modifier->SetMachineParametersFromGercEtAl();
		 simulator.AddSimulationModifier(p_modifier);

		 /* We then set an end time and run the simulation */
		 simulator.SetEndTime(2.2);
		 simulator.Solve();

		 //Tests
		 TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 8u);

		 /*
		 //Get
		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				 cell_iter != simulator.rGetCellPopulation().End();
				 ++cell_iter)
			{
				// Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineModifier>();
				boost::shared_ptr<TypeSixMachineModifier> p_modifier_test = boost::static_pointer_cast<TypeSixMachineModifier>(collection.GetProperty());
				unsigned number_machines = p_modifier_test->GetTotalNumberOfMachines();
				PRINT_VARIABLE(number_machines);
			}
*/

     }
//__________________________________________________________________________________________________________
// End of Test 02
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
// Start of Test 03
//__________________________________________________________________________________________________________
	void NoTestDoubleCapsuleDivisionMachinesLabels()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(2.5,2.5)));


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


		 // Create cells
		 std::vector<CellPtr> cells;
		 MAKE_PTR(WildTypeCellMutationState, p_state);
		 MAKE_PTR(TransitCellProliferativeType, p_type);
		 MAKE_PTR(TypeSixMachineProperty,p_property1);
		 MAKE_PTR(TypeSixMachineProperty,p_property2);
		 for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		 {
			 UniformCellCycleModel* p_model = new UniformCellCycleModel();
			 p_model->SetMinCellCycleDuration(1.0);
			 p_model->SetMaxCellCycleDuration(1.01);
			 CellPtr p_cell(new Cell(p_state, p_model));
			 p_cell->SetCellProliferativeType(p_type);

			if (i % 2 == 0)
			{
				p_cell->SetBirthTime(-0.1);
			} else
			{
				p_cell->SetBirthTime(-0.1);
			}



			//double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			//std::vector<double> machine_angles;
			//machine_angles.push_back(rand_angle);


			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);



			if (p_cell -> GetCellId() % 2 == 0)
				{

					p_property1->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_coordinates));
					p_property1->SetCellTypeLabel(0u);
					p_cell->AddCellProperty(p_property1);
				} else
				{

					p_property2->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_coordinates));
					p_property2->SetCellTypeLabel(1u);
					p_cell->AddCellProperty(p_property2);
				}


			cells.push_back(p_cell);

		 }



		 // Create cell population
		 NodeBasedCellPopulationCapsulesMachines<2> population(mesh, cells);


		 // Create writers to visualise
		 population.AddCellWriter<CellIdWriter>();
		 population.AddCellWriter<CapsuleOrientationWriter>();
		 population.AddCellWriter<CapsuleScalingWriter>();
		 population.AddCellWriter<CapsuleTypeLabelWriter>();
		 population.AddCellWriter<MachineStateCountWriter>();



		 // Create division rules
		 boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		 population.SetCentreBasedDivisionRule(p_division_rule);

		 // Create simulation
		 OffLatticeSimulation<2> simulator(population);
		 simulator.SetOutputDirectory("TestDoubleCapsuleDivisionMachinesLabels");
		 double dt = 1.0/1200.0;
		 simulator.SetDt(dt);
		 simulator.SetSamplingTimestepMultiple(20u);

		 auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		 simulator.SetNumericalMethod(p_numerical_method);

		 //Add capsule force law.
		 auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		 p_capsule_force->SetYoungModulus(200.0);
		 simulator.AddForce(p_capsule_force);

		 // Add capsule modifier with machines
		 MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		 p_modifier->SetOutputDirectory("TestDoubleCapsuleDivisionMachinesLabels");
		 //Set the variable for the states
		 //p_modifier->Setk_1(0.0);
		 //p_modifier->Setk_2(0.02);
		 //p_modifier->Setk_3(0.0);
		 //p_modifier->Setk_4(0.0);
		 //p_modifier->Setk_5(0.0);
		 //p_modifier->Setk_6(0.0);
		 //p_modifier->Setk_7(0.0);
		 p_modifier->SetMachineParametersFromGercEtAl();
		 simulator.AddSimulationModifier(p_modifier);

		 /* We then set an end time and run the simulation */
		 simulator.SetEndTime(3.0);
		 simulator.Solve();

		 //Tests
		 TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 16u);


	 }
//__________________________________________________________________________________________________________
// End of Test 03
//__________________________________________________________________________________________________________
//__________________________________________________________________________________________________________
// Start of Test 04
//__________________________________________________________________________________________________________
	void NoTestDoubleCapsuleDivisionMachinesLabelsKillers()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(0.0,2.0)));
		nodes.push_back(new Node<2>(2u, Create_c_vector(0.0,-2.0)));
		//nodes.push_back(new Node<2>(3u, Create_c_vector(-2.0,0.0)));
		//nodes.push_back(new Node<2>(4u, Create_c_vector(2.0,0.0)));


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

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
/*
		mesh.GetNode(3u)->AddNodeAttribute(0.0);
		mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(4u)->AddNodeAttribute(0.0);
		mesh.GetNode(4u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(4u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
*/

		 // Create cells
		 std::vector<CellPtr> cells;
		 MAKE_PTR(WildTypeCellMutationState, p_state);
		 MAKE_PTR(TransitCellProliferativeType, p_type);
		 MAKE_PTR(TypeSixMachineProperty,p_property1);
		 MAKE_PTR(TypeSixMachineProperty,p_property2);
		 for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		 {
			 UniformCellCycleModel* p_model = new UniformCellCycleModel();
			 p_model->SetMinCellCycleDuration(1.0);
			 p_model->SetMaxCellCycleDuration(1.6);
			 CellPtr p_cell(new Cell(p_state, p_model));
			 p_cell->SetCellProliferativeType(p_type);
			 double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			 std::vector<double> machine_angles;
			 machine_angles.push_back(rand_angle);


			if (i % 2 == 0)
			{
				p_cell->SetBirthTime(-0.1);
			} else
			{
				p_cell->SetBirthTime(-0.1);
			}
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;
/*
			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);
*/


			if (p_cell -> GetCellId() % 2 == 0)
				{

					p_property1->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
					p_property1->SetCellTypeLabel(0u);
					p_cell->AddCellProperty(p_property1);
				} else
				{

					p_property2->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
					p_property2->SetCellTypeLabel(1u);
					p_cell->AddCellProperty(p_property2);
				}


			cells.push_back(p_cell);

		 }



		 // Create cell population
		 NodeBasedCellPopulationCapsulesMachines<2> population(mesh, cells);


		 // Create writers to visualise
		 population.AddCellWriter<CellIdWriter>();
		 population.AddCellWriter<CapsuleOrientationWriter>();
		 population.AddCellWriter<CapsuleScalingWriter>();
		 population.AddCellWriter<CapsuleTypeLabelWriter>();
		 population.AddCellWriter<MachineStateCountWriter>();



		 // Create division rules
		 boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		 population.SetCentreBasedDivisionRule(p_division_rule);

		 // Create simulation
		 OffLatticeSimulation<2> simulator(population);
		 simulator.SetOutputDirectory("TestDoubleCapsuleDivisionMachinesLabelsKillers");
		 double dt = 1.0/1200.0;
		 simulator.SetDt(dt);
		 simulator.SetSamplingTimestepMultiple(10u);

		 auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		 simulator.SetNumericalMethod(p_numerical_method);

		 // Killing
		 MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
		 simulator.AddCellKiller(p_killer);

		 //Add capsule force law.
		 auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		 p_capsule_force->SetYoungModulus(200.0);
		 simulator.AddForce(p_capsule_force);

		 // Add capsule modifier with machines
		 MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		 p_modifier->SetOutputDirectory("TestDoubleCapsuleDivisionMachinesLabelsKillers");
		 //Set the variable for the states
		 //p_modifier->Setk_1(0.0);
		 //p_modifier->Setk_2(0.02);
		 //p_modifier->Setk_3(0.0);
		 //p_modifier->Setk_4(0.0);
		 //p_modifier->Setk_5(0.0);
		 //p_modifier->Setk_6(0.0);
		 //p_modifier->Setk_7(0.0);
		 p_modifier->SetMachineParametersFromGercEtAl();
		 simulator.AddSimulationModifier(p_modifier);

		 /* We then set an end time and run the simulation */
		 simulator.SetEndTime(2.0);
		 simulator.Solve();

		 //Tests
		 //TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 8u);


	 }
//__________________________________________________________________________________________________________
// End of Test 04
//__________________________________________________________________________________________________________

	void TestSingleCapsuleSimulationWithDivisionAndMachinesKillerGerc()
                		  {
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.0,2.5)));
		nodes.push_back(new Node<2>(2u, Create_c_vector(-3.0,-2.5)));
		nodes.push_back(new Node<2>(3u, Create_c_vector(-3.0,2.5)));

		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.0);
		c_vector<double, 4> domain_size;
		domain_size[0] = -2000.0;
		domain_size[1] = 2000.0;
		domain_size[2] = -2000.0;
		domain_size[3] = 2000.0;
		mesh.SetInitialBoxCollection(domain_size, 10.0);


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

		mesh.GetNode(2u)->AddNodeAttribute(0.0);
		mesh.GetNode(2u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(2u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(3u)->AddNodeAttribute(0.0);
		mesh.GetNode(3u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_THETA] = 0.5*M_PI;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(3u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;


		 // Create cells
		 std::vector<CellPtr> cells;
		 MAKE_PTR(WildTypeCellMutationState, p_state);
		 MAKE_PTR(TransitCellProliferativeType, p_type);
		 MAKE_PTR(TypeSixMachineProperty,p_property1);
		 MAKE_PTR(TypeSixMachineProperty,p_property2);
		 for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		 {
			 UniformCellCycleModel* p_model = new UniformCellCycleModel();
			 p_model->SetMinCellCycleDuration(1.0);
			 p_model->SetMaxCellCycleDuration(1.01);
			 CellPtr p_cell(new Cell(p_state, p_model));
			 p_cell->SetCellProliferativeType(p_type);
			 double rand_angle = 0.0;

			 //double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			 std::vector<double> machine_angles;
			 machine_angles.push_back(rand_angle);


			if (i % 2 == 0)
			{
				p_cell->SetBirthTime(-0.1);
			} else
			{
				p_cell->SetBirthTime(-0.1);
			}
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;
/*
			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);
*/


			if (p_cell -> GetCellId() % 2 == 0)
				{

					p_property1->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
					p_property1->SetCellTypeLabel(0u);
					p_cell->AddCellProperty(p_property1);
				} else
				{

					p_property2->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
					p_property2->SetCellTypeLabel(1u);
					p_cell->AddCellProperty(p_property2);
				}


			cells.push_back(p_cell);

		 }



		 // Create cell population
		 NodeBasedCellPopulationCapsulesMachines<2> population(mesh, cells);


		 // Create writers to visualise
		 population.AddCellWriter<CellIdWriter>();
		 population.AddCellWriter<CapsuleOrientationWriter>();
		 population.AddCellWriter<CapsuleScalingWriter>();
		 population.AddCellWriter<CapsuleTypeLabelWriter>();
		 population.AddCellWriter<MachineStateCountWriter>();


		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKillerGerc");
		double dt = 1.0/1200.0;
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(30u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);


		//
		MAKE_PTR_ARGS(TypeSixMachineCellKiller<2>, p_killer, (&population));
		simulator.AddCellKiller(p_killer);
		/*
		 * We now create a capsuleforce law and pass it to the simulation
		 */
		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		p_capsule_force->SetYoungModulus(200.0);

		simulator.AddForce(p_capsule_force);
		//

		MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		p_modifier->SetOutputDirectory("TestSingleCapsuleWithDivisionAndMachinesKillerGerc");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(8.0); // was 1.0075
		simulator.Solve();
		MARK;
		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
                		  }


};

#endif /*TESTCAPSULESDIVISIONMACHINES2D_HPP_*/
