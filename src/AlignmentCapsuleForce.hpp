
#ifndef ALIGNMENTCAPSULEFORCE_HPP_
#define ALIGNMENTCAPSULEFORCE_HPP_


#include "CapsuleForce.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law between two capsules (cylinder with hemispherical caps), defined in Farrell et al:
 * J R Soc Interface. 2017 Jun;14(131). pii: 20170073. doi: 10.1098/rsif.2017.0073
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AlignmentCapsuleForce : public CapsuleForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestAlignmentCapsuleForce;

private:

    /** The gamme coefficient for the bending stiffness energy */
    double mGamma;


    /** Needed for serialization. */
    friend class boost::serialization::access;



    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mGamma;

    }

public:

    /**
     * Constructor.
     */
    AlignmentCapsuleForce();

    /**
     * Destructor.
     */
    virtual ~AlignmentCapsuleForce() = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

    void SetGamma(double Gamma);
    double GetGamma();

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AlignmentCapsuleForce)

#endif /*ALIGNMENTCAPSULEFORCE_HPP_*/
