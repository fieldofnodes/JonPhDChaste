
#ifndef _TESTALILGNMENTCAPSULEFORCE_HPP_
#define _TESTALIGNMENTCAPSULEFORCE_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "CapsuleForce.hpp"
#include "AttractiveEndsCapsuleForce.hpp"
#include "AlignmentCapsuleForce.hpp"
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

class TestAlignmentCapsuleForce : public AbstractCellBasedTestSuite
{
public:

    void TestGenericEquailty() // throw(Exception)
    {
    	c_vector<double, 2u> VectorTest1;
    	c_vector<double, 2u> VectorTest2;
    	double dot;
    	VectorTest1(0) = -1;
    	VectorTest1(1) = 1;
    	VectorTest2(0) = 1;
    	VectorTest2(1) = 1;
    	dot = VectorTest1(0)*VectorTest2(0) + VectorTest1(1)*VectorTest2(1);
    	TS_ASSERT_DELTA(dot,0, 1e-9);
    }
}; // This bracket will close the class: TestAttractiveEndsCapsuleForce : public AbstractCellBasedTestSuite

#endif /*_TESTALIGNMENTSCAPSULEFORCE_HPP_*/

