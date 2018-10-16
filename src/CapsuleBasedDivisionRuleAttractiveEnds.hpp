#ifndef CAPSULEBASEDDIVISIONRULEATTRACTIVEENDS_HPP_
#define CAPSULEBASEDDIVISIONRULEATTRACTIVEENDS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class CapsuleBasedDivisionRuleAttractiveEnds : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
private:
	   /**
	     * The specified location of the new daughter cell.
	     * Initialized in the constructor.
	     */
	c_vector<double, SPACE_DIM> mDaughterLocation;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    /**
     * Default constructor.
     *
     * @param rDaughterLocation the specified location of the daughter cell
     */
    CapsuleBasedDivisionRuleAttractiveEnds(c_vector<double, SPACE_DIM>& rDaughterLocation);

    /**
     * Empty destructor.
     */
    virtual ~CapsuleBasedDivisionRuleAttractiveEnds()
    {
    }

    /**
     * @return mDaughterLocation.
     */
    const c_vector<double, SPACE_DIM>& rGetDaughterLocation() const;
    /**
     * Default constructor.
     */
    CapsuleBasedDivisionRuleAttractiveEnds()
    {
    }

    /**
     * Empty destructor.
     */


    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    virtual std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CapsuleBasedDivisionRuleAttractiveEnds)

#endif // CAPSULEBASEDDIVISIONRULEATTRACTIVEENDS_HPP_
