#ifndef EXTERNALFORCE_HH
#define EXTERNALFORCE_HH


#include "../../StrandSim/Core/Definitions.hh"

namespace bogus
{

/*!
  \brief Object-indexed vector.

  Consider a dense vector as the concatenation of smaller ones,
  allowing easy access of the part of the whole vector that corresponds to
  a given object.
  */
struct AggregateVec
{
    typedef Eigen::VectorXd::SegmentReturnType VecT ;
    typedef Eigen::VectorXd::ConstSegmentReturnType ConstVecT ;

    /*! Constructor
        \param vec the whole dense vector
        \param indices the indices of the first element of each object in \p vec
        */
    AggregateVec( Eigen::VectorXd& vec, const std::vector< unsigned > &indices )
        :m_vec( vec ), m_indices( indices )
    {
    }

    //! \returns the number of elements of the object \p i
    unsigned nDofs( unsigned i ) const
    {
        assert( i < m_indices.size() ) ;
        return ( i+1 == m_indices.size() ? m_vec.rows() : m_indices[ i+1 ] ) - m_indices[i] ;
    }

    //! \returns the segment corresponding to object \p i
    ConstVecT operator() ( const unsigned oId ) const
    {
        return const_cast< const Eigen::VectorXd& >( m_vec ).segment( m_indices[ oId ], nDofs( oId ) ) ;
    }

    //! \returns the segment corresponding to object \p i
    VecT operator() ( const unsigned oId )
    {
        return m_vec.segment( m_indices[ oId ], nDofs( oId ) ) ;
    }

    Eigen::VectorXd &vec() { return m_vec ; }
    const std::vector< unsigned > &indices() const { return m_indices ; }

    //! \returns the number of objects that compose the whole vector
    size_t size() const { return m_indices.size() ; }

    private:
        Eigen::VectorXd &m_vec ;
        const std::vector< unsigned > &m_indices ;
} ;


/*!
  \ingroup bogus
  \brief Allows to plug-in external forces than can be updated during the contact solving
*/
class ExternalForce
{
    public:

    ExternalForce( unsigned objId )
    : m_objectID( objId )
    {
    }

    virtual ~ExternalForce()
    {
    }

    //! Computes the current external force
    /*! \param velocities Current generalized velocities at this point of the contact solving
        \param solverForces Sum of all the forces currently added by the solver
        All vectors can be indexed by an object id to get the relevant segment
        */
    virtual bool compute ( const Eigen::VectorXd& velocities, const Eigen::VectorXd& solverForces ) = 0;

    unsigned m_objectID ;

} ;

}

#ifdef BOOST_CLASS_EXPORT
BOOST_CLASS_EXPORT(yacfs::ExternalForce)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(yacfs::ExternalForce)
#endif

#endif // EXTERNALFORCE_HH
