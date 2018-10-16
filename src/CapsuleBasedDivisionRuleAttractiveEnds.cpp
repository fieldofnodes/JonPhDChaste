#include "CapsuleBasedDivisionRuleAttractiveEnds.hpp"
#include "TypeSixSecretionEnumerations.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CapsuleBasedDivisionRuleAttractiveEnds<ELEMENT_DIM, SPACE_DIM>::CapsuleBasedDivisionRuleAttractiveEnds(c_vector<double, SPACE_DIM>& rDaughterLocation)
{
    mDaughterLocation = rDaughterLocation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& CapsuleBasedDivisionRuleAttractiveEnds<ELEMENT_DIM,SPACE_DIM>::rGetDaughterLocation() const
{
    return mDaughterLocation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CapsuleBasedDivisionRuleAttractiveEnds<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    c_vector<double, SPACE_DIM> axis_vector;
    //For the below DIM ==2 and DIM ==3, I changed the length to be a reference to the length of the cell.
    //Rather then the Distance = 1.5 -- what does that mean?
    switch (SPACE_DIM)
    {

        case 1:
        {
            EXCEPTION("CapsuleBasedDivisionRule is not implemented for SPACE_DIM==1");
        }
        case 2:
        {
        	Node<SPACE_DIM>* p_node = rCellPopulation.GetNodeCorrespondingToCell(pParentCell);

   	        const double orientation_theta = p_node->rGetNodeAttributes()[NA_THETA];
   	        const double length_capsule = p_node->rGetNodeAttributes()[NA_LENGTH];
   	        const double radius_capsule = p_node->rGetNodeAttributes()[NA_RADIUS];

        	const double distance=0.5*(length_capsule+2.0*radius_capsule);
            axis_vector(0) = distance*cos(orientation_theta);
            axis_vector(1) = distance*sin(orientation_theta);
            break;
        }
        case 3:
        {
        	Node<SPACE_DIM>* p_node = rCellPopulation.GetNodeCorrespondingToCell(pParentCell);

   	        const double orientation_theta = p_node->rGetNodeAttributes()[NA_THETA];
   	        const double orientation_phi = p_node->rGetNodeAttributes()[NA_PHI];
   	        const double length_capsule = p_node->rGetNodeAttributes()[NA_LENGTH];
   	        const double radius_capsule = p_node->rGetNodeAttributes()[NA_RADIUS];

        	const double distance=0.5*(length_capsule+2.0*radius_capsule);
            axis_vector(0) = distance*cos(orientation_theta)*sin(orientation_phi);
            axis_vector(1) = distance*sin(orientation_theta)*sin(orientation_phi);
            axis_vector(2) = distance*cos(orientation_phi);
            break;
        }
        default:
            NEVER_REACHED;
    }

    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - axis_vector;
    c_vector<double, SPACE_DIM> daughter_position = parent_position + 2.0*axis_vector;

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > positions(parent_position, daughter_position);
    return positions;
}

// Explicit instantiation
template class CapsuleBasedDivisionRuleAttractiveEnds<1,1>;
template class CapsuleBasedDivisionRuleAttractiveEnds<1,2>;
template class CapsuleBasedDivisionRuleAttractiveEnds<2,2>;
template class CapsuleBasedDivisionRuleAttractiveEnds<1,3>;
template class CapsuleBasedDivisionRuleAttractiveEnds<2,3>;
template class CapsuleBasedDivisionRuleAttractiveEnds<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleBasedDivisionRuleAttractiveEnds)
