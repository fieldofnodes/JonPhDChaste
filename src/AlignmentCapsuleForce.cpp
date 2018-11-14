
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
    	auto location_a = r_node_a.rGetLocation();
    	auto location_b = r_node_b.rGetLocation();
        const double length_a = r_node_a.rGetNodeAttributes()[NA_LENGTH];
        const double length_b = r_node_b.rGetNodeAttributes()[NA_LENGTH];        
        const double angle_theta_a = r_node_a.rGetNodeAttributes()[NA_THETA];
        const double angle_theta_b = r_node_b.rGetNodeAttributes()[NA_THETA];
        const double angle_phi_a = r_node_a.rGetNodeAttributes()[NA_PHI];
        const double angle_phi_b = r_node_b.rGetNodeAttributes()[NA_PHI];   
		double angle_between_a_b;
		
		c_vector<double,SPACE_DIM> segment_a_point_1;
		c_vector<double,SPACE_DIM> segment_a_point_2;
		c_vector<double,SPACE_DIM> segment_b_point_1;
		c_vector<double,SPACE_DIM> segment_b_point_2;
		
		if (SPACE_DIM==3u)
		{
			segment_a_point_1[0] = location_a[0] + 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
			segment_a_point_1[1] = location_a[1] + 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
			segment_a_point_1[2] = location_a[2] + 0.5 * length_a * cos(angle_phi_a);
		//        PRINT_VECTOR(segment_a_point_1);

			segment_a_point_2[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a) * sin(angle_phi_a);
			segment_a_point_2[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a) * sin(angle_phi_a);
			segment_a_point_2[2] = location_a[2] - 0.5 * length_a * cos(angle_phi_a);
		//        PRINT_VECTOR(segment_a_point_2);

			segment_b_point_1[0] = location_b[0] + 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
			segment_b_point_1[1] = location_b[1] + 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
			segment_b_point_1[2] = location_b[2] + 0.5 * length_b * cos(angle_phi_b);
		//        PRINT_VECTOR(segment_b_point_1);

			segment_b_point_2[0] = location_b[0] - 0.5 * length_b * cos(angle_theta_b) * sin(angle_phi_b);
			segment_b_point_2[1] = location_b[1] - 0.5 * length_b * sin(angle_theta_b) * sin(angle_phi_b);
			segment_b_point_2[2] = location_b[2] - 0.5 * length_b * cos(angle_phi_b);
		}
		else
		{
			segment_a_point_1[0] = location_a[0] + 0.5 * length_a * cos(angle_theta_a);
			segment_a_point_1[1] = location_a[1] + 0.5 * length_a * sin(angle_theta_a);
		//        PRINT_VECTOR(segment_a_point_1);

			segment_a_point_2[0] = location_a[0] - 0.5 * length_a * cos(angle_theta_a);
			segment_a_point_2[1] = location_a[1] - 0.5 * length_a * sin(angle_theta_a);
		//        PRINT_VECTOR(segment_a_point_2);

			segment_b_point_1[0] = location_b[0] + 0.5 * length_b * cos(angle_theta_b);
			segment_b_point_1[1] = location_b[1] + 0.5 * length_b * sin(angle_theta_b);
		//        PRINT_VECTOR(segment_b_point_1);

			segment_b_point_2[0] = location_b[0] - 0.5 * length_b * cos(angle_theta_b);
			segment_b_point_2[1] = location_b[1] - 0.5 * length_b * sin(angle_theta_b);
		}
		
		
		c_vector<double,SPACE_DIM> vector_a;
		c_vector<double,SPACE_DIM> vector_b;
		
		vector_a = segment_a_point_1-segment_a_point_2;
		vector_b = segment_b_point_1-segment_b_point_2;
		
		angle_between_a_b = acos(inner_prod(vector_a,vector_b)/(norm_2(vector_a)*norm_2(vector_b)));
		//PRINT_VARIABLE(angle_between_a_b);
		
		
		
		if (SPACE_DIM==2u && angle_between_a_b < M_PI/2.0)
			
			{
				r_node_a.rGetNodeAttributes()[NA_APPLIED_THETA] += mGamma*sin(2*(angle_theta_a-angle_theta_b));
				r_node_b.rGetNodeAttributes()[NA_APPLIED_THETA] += (-1.0)*mGamma*sin(2*(angle_theta_a-angle_theta_b));
			}
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
