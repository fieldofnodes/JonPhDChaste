
#ifndef _TESTCAPSULELABELTOCOLOURDIVISION_HPP_
#define _TESTCAPSULELABELTOCOLOURDIVISION_HPP_


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

class TestCapsuleLabelToColourDivision : public AbstractCellBasedTestSuite
{
public:
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Labelling cells --  division
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________

	void TestCellLabelsToColourDivision()
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

			if (i % 2 == 0)
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


			// if the node number is even then the label is 0u -- meaning attacker, if the node number is off
			// then the label is 1u -- meaning not attacker
			PRINT_VECTOR(mesh.GetAllNodeIndices());
			if (p_cell -> GetCellId() == 0)
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
		simulator.SetOutputDirectory("03TestCellLabelsToColourDivision");
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
		simulator.SetEndTime(1.0);
		simulator.Solve();
		//Determining the population size is the initial size, say p_0+floor[p_0*min_cell_cycle_duration] - where
		// the floor function is the greatest integer less than p_0*p_0*min_cell_cycle_duration
		PRINT_VARIABLE(simulator.rGetCellPopulation().GetNumRealCells());
		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),4u);
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


		        	//Currently this test passes at less than 4 cells. Based off of testing for inheritance --
		        	// this population is labelled as needed via the property labelling pointers, thought I
		        	// need to run a better test to match the original 2 nodes against the new populatin.
		        	// For some reason the Cell Id given in division does not go evenly to each cell or something
		        	// that causes my modulus 2 test to fail
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
}; 

#endif /*_TESTCAPSULELABELTOCOLOURDIVISION_HPP_*/

