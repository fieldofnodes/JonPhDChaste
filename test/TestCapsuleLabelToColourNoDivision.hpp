
#ifndef _TESTCAPSULELABELTOCOLOURNODIVISION_HPP_
#define _TESTCAPSULELABELTOCOLOURNODIVISION_HPP_


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
#include "NodeBasedCellPopulationCapsules.hpp"
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

class TestCapsuleLabelToColourNoDivision : public AbstractCellBasedTestSuite
{
public:
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
// Labelling cells -- no division
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________________

	void TestCellLabelsToColourAttackerNoDivision()
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
		NodeBasedCellPopulationCapsules<2> population(mesh, cells);

		population.AddCellWriter<CellIdWriter>();
		population.AddCellWriter<CapsuleOrientationWriter>();
		population.AddCellWriter<CapsuleScalingWriter>();
		population.AddCellWriter<CapsuleTypeLabelWriter>();

		// Create simulation
		OffLatticeSimulation<2> simulator(population);
		simulator.SetOutputDirectory("02TestCellLabelsToColourNoDivision");
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

}; 

#endif /*_TESTCAPSULELABELTOCOLOURNODIVISION_HPP_*/

