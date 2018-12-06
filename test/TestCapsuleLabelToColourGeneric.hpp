
#ifndef _TESTCAPSULELABELTOCOLOUR_HPP_
#define _TESTCAPSULELABELTOCOLOUR_HPP_


#include <cxxtest/TestSuite.h>
#include <cycle/UniformCellCycleModel.hpp>
#include <boost/shared_ptr.hpp>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "CheckpointArchiveTypes.hpp"


#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellMutationState.hpp"
#include "AbstractCellProperty.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "CellPropertyRegistry.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "OutputFileHandler.hpp"
#include "NoCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Cell.hpp"
#include "CellLabel.hpp"
#include "CellIdWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"



// Header files included in this project
#include "TypeSixSecretionEnumerations.hpp"
#include "ForwardEulerNumericalMethodForCapsules.hpp"
#include "CapsuleForce.hpp"
#include "AttractiveEndsCapsuleForce.hpp"
#include "AttractiveEndsAlignmentCapsuleForce.hpp"
#include "AlignmentCapsuleForce.hpp"
#include "CapsuleOrientationWriter.hpp"
#include "CapsuleScalingWriter.hpp"
#include "CapsuleTypeLabelWriter.hpp"
#include "SquareBoundaryCondition.hpp"
#include "CapsuleBasedDivisionRule.hpp"
#include "TypeSixMachineModifier.hpp"
#include "NodeBasedCellPopulationWithCapsules.hpp"
#include "TypeSixMachineProperty.hpp"
#include "TypeSixMachineCellKiller.hpp"
#include "TypeSixMachineCellLabelledKiller.hpp"
#include "MachineStateCountWriter.hpp"



// Should usually be called last.
#include "FakePetscSetup.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Testing the variety of label interactions on Capsules
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________

class TestCapsuleLabelToColour : public AbstractCellBasedTestSuite
{
private:
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Randomly label calsules  -- Not really working at the moment
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
	void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<TypeSixMachineProperty> pLabel1, boost::shared_ptr<TypeSixMachineProperty> pLabel2, double labelledRatio)
    {
		// 0 is Attacker and 1 is not attacker

    	for (unsigned i = 0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
            	pLabel1 -> SetCellTypeLabel(0u);
            	rCells[i]->AddCellProperty(pLabel1);

            }
            else
            {
            	pLabel2 -> SetCellTypeLabel(1u);
            	rCells[i]->AddCellProperty(pLabel2);

            }
        }
    }
public:
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Labelling cells -- no division
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________

	void xTestCellLabelsToColourAttacker()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(TypeSixMachineProperty, p_property);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);
        p_property1 -> SetCellTypeLabel(0u);
        cells[0]->AddCellProperty(p_property1);
        p_property2 -> SetCellTypeLabel(1u);
        cells[1]->AddCellProperty(p_property2);




		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourAttacker");
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

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
		         cell_iter != simulator.rGetCellPopulation().End();
		         ++cell_iter)
		    {
		        // Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
		        if (collection.GetSize() != 1)
		        {
		            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
		        }
		        	boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
		        	unsigned cell_type = p_property_test->GetCellTypeLabel();
		        	if (cell_iter -> GetCellId() == 0)
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
		        	} else
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
		        	}


		    }


	}

	void xTestCellLabelsToColourNotAttacker()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(TypeSixMachineProperty, p_property);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);
        p_property1 -> SetCellTypeLabel(0u);
        cells[0]->AddCellProperty(p_property1);
        p_property2 -> SetCellTypeLabel(0u);
        cells[1]->AddCellProperty(p_property2);


		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourNotAttacker");
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

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();
		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
		         cell_iter != simulator.rGetCellPopulation().End();
		         ++cell_iter)
		    {
		        // Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
		        if (collection.GetSize() != 1)
		        {
		            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
		        }
		        	boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
		        	unsigned cell_type = p_property_test->GetCellTypeLabel();
		        	if (cell_iter -> GetCellId() == 0)
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
		        	} else
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
		        	}


		    }
	}
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
// Machine no division
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Labelling cells -- no division
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________

	void xTestCellLabelsToColourAttackerMachine()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(TypeSixMachineProperty, p_property);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);
		p_property1 -> SetCellTypeLabel(0u);
		cells[0]->AddCellProperty(p_property1);
		p_property2 -> SetCellTypeLabel(1u);
		cells[1]->AddCellProperty(p_property2);




		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();
		population.AddCellWriter<MachineStateCountWriter>();


		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourAttackerMachines");
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

		MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		p_modifier->SetOutputDirectory("TestCellLabelsToColourAttackerMachines");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(100.0/1200.0);
		simulator.Solve();

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				 cell_iter != simulator.rGetCellPopulation().End();
				 ++cell_iter)
			{
				// Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
				if (collection.GetSize() != 1)
				{
					EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
				}
					boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
					unsigned cell_type = p_property_test->GetCellTypeLabel();
					if (cell_iter -> GetCellId() == 0)
					{
						TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
						PRINT_VARIABLE(cell_type);
						PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
					} else
					{
						TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
						PRINT_VARIABLE(cell_type);
						PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
					}


			}


	}


//___________________________________________________________________________________________________________________________________
// Just focus on this test with the division
// It compiles but something is wrong with the test -- memory stuff

	void NoTestCellLabelsToColourNotAttackerDivision()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.01);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);

			if (i==0)
			{
				p_cell->SetBirthTime(-0.1);
			} else
			{
				p_cell->SetBirthTime(-0.1);
			}

			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);




			//boost::shared_ptr<TypeSixMachineProperty> p_property;



			if (p_cell -> GetCellId() % 2 == 0)
				{

					p_property1->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_coordinates));
					p_property1->SetCellTypeLabel(0u);
					p_cell->AddCellProperty(p_property1);
				} else
				{

					p_property2->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_coordinates));
					p_property2->SetCellTypeLabel(1u);
					p_cell->AddCellProperty(p_property2);
				}


	        cells.push_back(p_cell);

		}


		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
				 population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourNotAttackerDivision");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(10u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */

		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		simulator.AddForce(p_capsule_force);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(8.0);
		simulator.Solve();

		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),3u);
		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
		         cell_iter != simulator.rGetCellPopulation().End();
		         ++cell_iter)
		    {
		        // Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
		        if (collection.GetSize() != 1)
		        {
		            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
		        }
		        	boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
		        	unsigned cell_type = p_property_test->GetCellTypeLabel();
		        	if (cell_iter -> GetCellId() % 2 == 0)
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
		        	} else
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
		        	}


		    }
	}
//________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________
// Just focus on this test with the division  Machines
// It compiles but something is wrong with the test -- memory stuff
//________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________



	void TestCellLabelsToColourNotAttackerMachineDivision()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = .25*M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

/*
		for (unsigned node_idx = 1; node_idx < mesh.GetNumNodes(); ++node_idx)
		{
			mesh.GetNode(node_idx)->AddNodeAttribute(0.0);
			mesh.GetNode(node_idx)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_THETA] = .25*M_PI;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
			mesh.GetNode(node_idx)->rGetNodeAttributes()[NA_RADIUS] = 0.5;
		}
*/
		//Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			UniformCellCycleModel* p_model = new UniformCellCycleModel();
			p_model->SetMinCellCycleDuration(1.0);
			p_model->SetMaxCellCycleDuration(1.01);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);

			if (i % 2 == 0)
			{
				p_cell->SetBirthTime(-0.9);
			} else
			{
				p_cell->SetBirthTime(-0.9);
			}

			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;


			double vertical_coordinate = 0.25*(mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH]);
			double azimuthal_coordinate = M_PI ;


			std::vector<double> machine_coordinates;
			machine_coordinates.push_back(vertical_coordinate);
			machine_coordinates.push_back(azimuthal_coordinate);

			double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			std::vector<double> machine_angles;
			machine_angles.push_back(rand_angle);


			if (p_cell -> GetCellId() % 2 == 0)
				{

					p_property1->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1, machine_coordinates));
					p_property1->SetCellTypeLabel(0u);
					p_cell->AddCellProperty(p_property1);
				} else
				{

					p_property2->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1, machine_coordinates));
					p_property2->SetCellTypeLabel(1u);
					p_cell->AddCellProperty(p_property2);
				}


			cells.push_back(p_cell);

		}


		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();
		population.AddCellWriter<MachineStateCountWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourNotAttackerMachineDivision");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(10u);

		auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		simulator.SetNumericalMethod(p_numerical_method);



		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */

		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		simulator.AddForce(p_capsule_force);

		MAKE_PTR(TypeSixMachineModifier<2>, p_modifier);
		p_modifier->SetOutputDirectory("TestCellLabelsToColourNotAttackerMachineDivision");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);

		/* We then set an end time and run the simulation */
		simulator.SetEndTime(3.0);
		simulator.Solve();

		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),3u);
		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				 cell_iter != simulator.rGetCellPopulation().End();
				 ++cell_iter)
			{
				// Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
				if (collection.GetSize() != 1)
				{
					EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
				}
					boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
					unsigned cell_type = p_property_test->GetCellTypeLabel();
					if (cell_iter -> GetCellId() % 2 == 0)
					{
						TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
						PRINT_VARIABLE(cell_type);
						PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
					} else
					{
						TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
						PRINT_VARIABLE(cell_type);
						PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
					}


			}
	}
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//Gerc Test with Apoptosis
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________________
	void xTestSingleCapsuleSimulationWithDivisionAndMachinesKillerGercVersionTwo()
                		  {
		EXIT_IF_PARALLEL;

		//auto p_rand_gen = RandomNumberGenerator::Instance();

		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(5.0, 5.0)));

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
			p_model->SetMaxCellCycleDuration(1.6);
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_type);


			double rand_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf()-M_PI;
			std::vector<double> machine_angles;
             machine_angles.push_back(rand_angle);
             MAKE_PTR(TypeSixMachineProperty, p_property);
             p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(1u, machine_angles));
             p_cell->AddCellProperty(p_property);

			//double birth_time = -RandomNumberGenerator::Instance()->ranf();
			p_cell->SetBirthTime(-0.9);
			mesh.GetNode(i)->rGetNodeAttributes()[NA_LENGTH] = 2.0 +3.0*p_cell->GetBirthTime()/p_model->GetCellCycleDuration(); ;

			cells.push_back(p_cell);


		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<MachineStateCountWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
		population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestSingleCapsuleSimulationWithDivisionAndMachinesKillerGercVersionTwo");
		double dt = 1.0/1200.0;
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(10);

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
		p_modifier->SetOutputDirectory("TestSingleCapsuleSimulationWithDivisionAndMachinesKillerGercVersionTwo");
		p_modifier->SetMachineParametersFromGercEtAl();
		simulator.AddSimulationModifier(p_modifier);


		/* We then set an end time and run the simulation */
		simulator.SetEndTime(1.0); // was 1.0075
		simulator.Solve();
		MARK;
		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
	  }








	//___________________________________________________________________________________________________________________________________
	void xTestCellLabelsToColourAttackerDivision()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		MAKE_PTR(TypeSixMachineProperty, p_property);
		MAKE_PTR(TypeSixMachineProperty, p_property1);
		MAKE_PTR(TypeSixMachineProperty, p_property2);

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



			p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_coordinates));
	       /*
			if (i == 0)
	        {
	        	p_property1 -> SetCellTypeLabel(0u);
	        	cells[i]->AddCellProperty(p_property1);
	        } else
	        {
		        p_property2 -> SetCellTypeLabel(1u);
		        cells[i]->AddCellProperty(p_property2);
	        }
	        */

			cells.push_back(p_cell);
		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
				 population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourAttackerDivision");
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



		/* We then set an end time and run the simulation */
		simulator.SetEndTime(6.0);
		MARK;
		simulator.Solve();

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
		         cell_iter != simulator.rGetCellPopulation().End();
		         ++cell_iter)
		    {
		        // Get this cell's type six machine property data
				CellPropertyCollection collection = cell_iter->rGetCellPropertyCollection().template GetProperties<TypeSixMachineProperty>();
		        if (collection.GetSize() != 1)
		        {
		            EXCEPTION("TypeSixMachineModifier cannot be used unless each cell has a TypeSixMachineProperty");
		        }
		        	boost::shared_ptr<TypeSixMachineProperty> p_property_test = boost::static_pointer_cast<TypeSixMachineProperty>(collection.GetProperty());
		        	unsigned cell_type = p_property_test->GetCellTypeLabel();
		        	if (cell_iter -> GetCellId() == 0)
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property1 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property1 -> GetCellTypeLabel());
		        	} else
		        	{
		        		TS_ASSERT_EQUALS(cell_type, p_property2 -> GetCellTypeLabel());
		        		PRINT_VARIABLE(cell_type);
		        		PRINT_VARIABLE(p_property2 -> GetCellTypeLabel());
		        	}


		    }
	}

	void NoTestCellLabelsToColourAttackerAndNoTAttacker()
	{
		EXIT_IF_PARALLEL;

		// Create some capsules
		std::vector<Node<2>*> nodes;

		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(2.0, 2.0)));


		/*
		 * We then convert this list of nodes to a `NodesOnlyMesh`,
		 * which doesn't do very much apart from keep track of the nodes.
		 */
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 150.5);

		mesh.GetNode(0u)->AddNodeAttribute(0.0);
		mesh.GetNode(0u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_THETA] = 0.25*M_PI;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(0u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;

		mesh.GetNode(1u)->AddNodeAttribute(0.0);
		mesh.GetNode(1u)->rGetNodeAttributes().resize(NA_VEC_LENGTH);
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_THETA] = 0.25*M_PI;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_LENGTH] = 2.0;
		mesh.GetNode(1u)->rGetNodeAttributes()[NA_RADIUS] = 0.5;





		//Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(WildTypeCellMutationState, p_state);
        // Randomly label some cells

		// 0 is Attacker and 1 is not attacker
		//p_property -> SetCellTypeLabel(1);


		auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
		CellsGenerator<NoCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);
		MAKE_PTR(TypeSixMachineProperty, pLabel1);
		MAKE_PTR(TypeSixMachineProperty, pLabel2);
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<TypeSixMachineProperty>());
        pLabel1 -> SetCellTypeLabel(1u);
        cells[0]->AddCellProperty(pLabel1);
        pLabel2 -> SetCellTypeLabel(0u);
        cells[1]->AddCellProperty(pLabel2);



		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourAttackerAndNoTAttacker");
		simulator.SetDt(1.0/1200.0);
		simulator.SetSamplingTimestepMultiple(30u);

		//auto p_numerical_method = boost::make_shared<ForwardEulerNumericalMethodForCapsules<2,2>>();
		//simulator.SetNumericalMethod(p_numerical_method);

		/*
		 * We now create a force law and pass it to the simulation
		 * We use linear springs between cells up to a maximum of 1.5 ('relaxed' cell diameters) apart, and add this to the simulation class.
		 */

		auto p_capsule_force = boost::make_shared<CapsuleForce<2>>();
		simulator.AddForce(p_capsule_force);



		/* We then set an end time and run the simulation */
		simulator.SetEndTime(6.0);
		simulator.Solve();
		//PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
		 for (typename AbstractCellPopulation<2>::Iterator cell_iter = population.Begin();
		         cell_iter != population.End();
		         ++cell_iter)
		    {
			 PRINT_VARIABLE(pLabel1 -> GetCellTypeLabel());
			 TS_ASSERT_EQUALS(pLabel1 -> GetCellTypeLabel(), 1u);
			 PRINT_VARIABLE(pLabel2 -> GetCellTypeLabel());
			 TS_ASSERT_EQUALS(pLabel2 -> GetCellTypeLabel(), 0u);

			}
	}
	void NoTestCellLabelsToColourAttackerAndNotAttackerDivision()
	{
		EXIT_IF_PARALLEL;
		// Create some capsules
		std::vector<Node<2>*> nodes;
		nodes.push_back(new Node<2>(0u, Create_c_vector(0.0, 0.0)));
		nodes.push_back(new Node<2>(1u, Create_c_vector(1.5,1.5)));


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
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(TransitCellProliferativeType, p_type);
		MAKE_PTR(TypeSixMachineProperty, p_property);


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
			//MAKE_PTR(TypeSixMachineProperty, p_property2);
			p_property->rGetMachineData().emplace_back(std::pair<unsigned, std::vector<double>>(4, machine_coordinates));

			if (i==0)
			{

				p_property -> SetCellTypeLabel(0u);
			}
			else
			{
				p_property -> SetCellTypeLabel(1u);
			}
			p_cell->AddCellProperty(p_property);


			cells.push_back(p_cell);
		}

		// Create cell population
		NodeBasedCellPopulationWithCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CapsuleBasedDivisionRule<2,2>());
				 population.SetCentreBasedDivisionRule(p_division_rule);

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("TestCellLabelsToColourAttackerAndNotAttackerDivision");
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



		/* We then set an end time and run the simulation */
		simulator.SetEndTime(6.0);
		simulator.Solve();

		//PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
		 for (typename AbstractCellPopulation<2>::Iterator cell_iter = population.Begin();
		         cell_iter != population.End();
		         ++cell_iter)
		    {
			 PRINT_VARIABLE(p_property -> GetCellTypeLabel());
			 TS_ASSERT_EQUALS(p_property -> GetCellTypeLabel(), 0u);
			}
	}

}; 

#endif /*_TESTCAPSULELABELTOCOLOUR_HPP_*/

