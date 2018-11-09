
#include "CapsuleForce.hpp"
#include "AlignmentCapsuleForce.hpp"
#include "TypeSixSecretionEnumerations.hpp"
#include "NodeBasedCellPopulation.hpp"

#ifdef CHASTE_VTK
#include <vtkLine.h>
#endif // CHASTE_VTK

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <cmath>

#include "Debug.hpp"



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AlignmentCapsuleForce<ELEMENT_DIM, SPACE_DIM>::SetGamma(double Gamma)

{
	mGamma=Gamma;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AlignmentCapsuleForce<ELEMENT_DIM, SPACE_DIM>::GetGamma()

{
	return mGamma;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AlignmentCapsuleForce<ELEMENT_DIM, SPACE_DIM>::AlignmentCapsuleForce()
        : CapsuleForce<ELEMENT_DIM, SPACE_DIM>(),
          mGamma(100.0)

{
    // Has to be either element and space dimensions are both 2 or both 3.
    assert((ELEMENT_DIM == 2u && SPACE_DIM == 2u) || (ELEMENT_DIM == 3u && SPACE_DIM == 3u));
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AlignmentCapsuleForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    auto p_cell_population = dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);

    if (p_cell_population == nullptr)
    {
        EXCEPTION("Capsule force only works with AbstractCentreBasedCellPopulation");
    }


    // Calculate force and applied angle contributions from each pair
    for (auto& node_pair : p_cell_population->rGetNodePairs())
    {
        Node<SPACE_DIM>& r_node_a = *(node_pair.first);
        Node<SPACE_DIM>& r_node_b = *(node_pair.second);
        const double angle_theta_a = r_node_a.rGetNodeAttributes()[NA_THETA];
        const double angle_theta_b = r_node_b.rGetNodeAttributes()[NA_THETA];
        //const double angle_phi_a = r_node_a.rGetNodeAttributes()[NA_PHI];
        //const double angle_phi_b = r_node_b.rGetNodeAttributes()[NA_PHI];   
		double angle_between_a_b;
		angle_between_a_b = angle_theta_a-angle_theta_b;
		//PRINT_VARIABLE(fabs(angle_between_a_b));
		//if (SPACE_DIM==2u && angle_between_a_b < M_PI/2.0)
			
			//{
				r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += mGamma*sin(2*(angle_between_a_b));
				r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += (-1.0)*mGamma*sin(2*(angle_between_a_b));
			//}
			//else if (SPACE_DIM==2u && angle_between_a_b > M_PI/2.0)
			{
				/*
				r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product_3d(torque_vec_a, force_b_a)[2];
				r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += cross_product_3d(torque_vec_b, force_a_b)[2];
				r_node_a.rGetNodeAttributes()[NA_APPLIED_PHI] -= cross_product_3d(torque_vec_a, force_b_a)[0];
				r_node_b.rGetNodeAttributes()[NA_APPLIED_PHI] -= cross_product_3d(torque_vec_b, force_a_b)[0];
				*/
			}
		
	}
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AlignmentCapsuleForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class AlignmentCapsuleForce<1,1>;
template class AlignmentCapsuleForce<1,2>;
template class AlignmentCapsuleForce<2,2>;
template class AlignmentCapsuleForce<1,3>;
template class AlignmentCapsuleForce<2,3>;
template class AlignmentCapsuleForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AlignmentCapsuleForce)
