#ifndef ELEMENTPROXY_HH_
#define ELEMENTPROXY_HH_

#include "../Mesh/TriMesh.h"
#include "../Mesh/TriMeshController.h"
#include "../Strand/StrandDynamics.h"
#include "../Strand/ElasticStrand.h"
#include "../Simulation/ImplicitStepper.h"
#include "CollisionUtils/BoundingBox.hh"
#include "TwistEdgeHandler.h"
#include <stack>

typedef BoundingBox<float> BBoxType;
class TwistEdgeHandler;

/* [H]
    If 'statique' then Bounding Box only involves start of time positions
        (assuming no motion)
    Else, if !statique, then Bounding Box expanded to include end of time
        positions in volume as well, in a cheap attempt to create swept volume BBox
*/

class ElementProxy
{
public:
    ElementProxy()
    {}

    virtual ~ElementProxy()
    {}

    virtual void computeBoundingBox( BBoxType& boundingBox, bool statique = false, TwistEdgeHandler* teh = NULL ) const = 0;

    void updateBoundingBox( bool statique = false, TwistEdgeHandler* teh = NULL )
    {
        computeBoundingBox( m_boundingBox, statique, teh );
    }

    const BBoxType& getBoundingBox() const
    { return m_boundingBox; }

    void resetBoundingBox()
    { m_boundingBox.reset(); }

    friend std::ostream& operator<<( std::ostream& os, const ElementProxy& elem );

protected:
    virtual void print( std::ostream& os ) const = 0;

    BBoxType m_boundingBox;
};


class EdgeProxy: public ElementProxy
{
public:
    EdgeProxy( ElasticStrand& strand, int vertexIndex, ImplicitStepper* stepper ):
        m_implicit( stepper ), 
        m_strand( strand ), 
        m_vertexIndex( vertexIndex )
    {
        std::cout << "ideally we dont need ImplicitStepper here anymore, just grab all dynamics from strand state and strand dynamics which we can get from strand" << std::endl;
    }

    void computeBoundingBox( BBoxType& boundingBox, bool statique ) const
    { 
        boundingBox.reset();
        boundingBox.insert( m_strand.getFutureVertex( m_vertexIndex ).cast< float >() );
        boundingBox.insert( m_strand.getFutureVertex( m_vertexIndex + 1).cast< float >() );

        if( !statique )
        {
            const Vec3f& vtx0 = m_strand.getVertex( m_vertexIndex ).cast< float >();
            const Vec3f& vtx1 = m_strand.getVertex( m_vertexIndex + 1 ).cast< float >();

            boundingBox.insert( vtx0 );
            boundingBox.insert( vtx1 );
        }
    }

    bool operator <( EdgeProxy* other ) const
    {
        if( m_strand.getGlobalIndex() == other->getStrandPointer()->getGlobalIndex() ){
            return m_vertexIndex < other->getVertexIndex();
        }
        else{
            return m_strand.getGlobalIndex() < other->getStrandPointer()->getGlobalIndex();        
        }
    }

    const ElasticStrand* getStrandPointer() const
    { return &m_strand; }

    ElasticStrand* getStrandPointer()
    { return &m_strand; }

    int getVertexIndex() const
    { return m_vertexIndex; }

    ElasticStrand& getStrand()
    { return m_strand; }

    ImplicitStepper* getImplicitStepper() const
    { return m_implicit; }

protected:
    virtual void print( std::ostream& os ) const;
    ImplicitStepper* m_implicit;
    ElasticStrand& m_strand;
    int m_vertexIndex;
};

class CylinderProxy: public EdgeProxy
{
public:
    CylinderProxy( ElasticStrand& strand, int vertexIndex, ImplicitStepper* stepper ) :
        EdgeProxy( strand, vertexIndex, stepper )
    {
        m_radius = m_strand.collisionParameters().getLargerCollisionsRadius( m_vertexIndex );
    }

    void computeBoundingBox( BBoxType& boundingBox, bool statique, TwistEdgeHandler* teh = NULL ) const
    {
        Vec3f min, max;
        static const Vec3f unit = Vec3f::Ones();
        EdgeProxy::computeBoundingBox( boundingBox, statique );

        min = boundingBox.min - m_radius * unit;
        max = boundingBox.max + m_radius * unit;

        boundingBox.insert( min );
        boundingBox.insert( max );
    }

    Scalar m_radius;
};


struct TwistIntersection
{
    TwistIntersection( Scalar refAngle )
    {
        originalTwistAngle = refAngle;
        coplanarTwists = 0;
        currAngle = 365;
    }
    ~TwistIntersection(){}

    Scalar originalTwistAngle;
    Scalar currAngle; // do we need to store this? its always computed on the fly...?
    int coplanarTwists;
};


class TwistEdge: public CylinderProxy
{
public:
    TwistEdge( ElasticStrand& strand, int vertexIndex, ImplicitStepper* stepper );
    TwistEdge( TwistEdge* first, TwistEdge* second );
    TwistEdge( TwistEdge* first, TwistEdge* second, int& uID );

    Scalar currAngle() const;
    int coplanarTwists() const;
    int intersectionTwists() const;
    Scalar originalTwistAngle() const;

    void computeBoundingBox( BBoxType& boundingBox, bool statique, TwistEdgeHandler* teh ) const;

    static void printIntersections( std::stack< TwistIntersection* >& intersections, std::ofstream& out );

    bool isTwistedBand;
    int uniqueID;
    int frozenID;

    std::stack< TwistIntersection* > intersections;

    bool flagged;
    bool traversed;

    std::vector< std::pair< int, TwistEdge* > > children; // partner(uniqueID) - child(TwistEdge*)
    std::pair< TwistEdge*, TwistEdge* > parents;

    TwistEdge* prev;
    TwistEdge* next;
};







//%//////////////////////////^



class FaceProxy: public ElementProxy
{
public:
    FaceProxy( const unsigned faceIndex, TriMeshController* controller ):
        m_faceIndex( faceIndex ),
        m_face( controller->getMesh()->getFace( faceIndex ) ),
        m_controller( controller )
    {}

    void computeBoundingBox( BBoxType& boundingBox, bool statique, TwistEdgeHandler* teh = NULL ) const
    {
        boundingBox.reset();

        const Vec3f vtx0 = getVertex( 0 ).cast< float >();
        const Vec3f vtx1 = getVertex( 1 ).cast< float >();
        const Vec3f vtx2 = getVertex( 2 ).cast< float >();

        boundingBox.insert( vtx0 - getDisplacement( 0 ).cast< float >() );
        boundingBox.insert( vtx1 - getDisplacement( 1 ).cast< float >() );
        boundingBox.insert( vtx2 - getDisplacement( 2 ).cast< float >() );

        if ( !statique )
        {
            boundingBox.insert( vtx0 );
            boundingBox.insert( vtx1 );
            boundingBox.insert( vtx2 );
        }
    }

    Vec3 getVertex( short apex ) const
    {
        return m_controller->getMesh()->getVertex( m_face.idx[apex] );
    }

    Vec3 getDisplacement( short apex ) const
    {
        return m_controller->getMesh()->getDisplacement( m_face.idx[apex] );
    }

    int getVertexIdx( short apex ) const
    {
        return m_face.idx[apex];
    }

    bool hasEnabledApex( short apex ) const
    {
        const std::vector<bool>& enabledVertices = m_controller->getEnabledVertices();

        if( enabledVertices.empty() ){
            return true;
        }
        else{
            return enabledVertices[m_face.idx[apex]];
        }
    }

    bool allApicesEnabled() const
    {
        const std::vector<bool>& enabledVertices = m_controller->getEnabledVertices();

        if( enabledVertices.empty() ){
            return true;
        }
        else{
            return enabledVertices[m_face.idx[0]] && enabledVertices[m_face.idx[1]]
                    && enabledVertices[m_face.idx[2]];
        }
    }

    double getFrictionCoefficient( double b0, double b1, double b2 ) const
    {
        const std::vector<double>& frictionCoefficients = m_controller->getFrictionCoefficients();

        if ( frictionCoefficients.empty() )
        {
            return m_controller->getDefaultFrictionCoefficient();
        }
        else
        {
            return b0 * frictionCoefficients[m_face.idx[0]]
                    + b1 * frictionCoefficients[m_face.idx[1]]
                    + b2 * frictionCoefficients[m_face.idx[2]];
        }
    }

    const TriMesh* getMesh() const
    {
        return m_controller->getMesh();
    }

    Vec3 getNormal() const
    {
        const Vec3 &q0 = getVertex( 0 );
        const Vec3 &q1 = getVertex( 1 );
        const Vec3 &q2 = getVertex( 2 );

        return ( q1 - q0 ).cross( q2 - q0 ).normalized();
    }

    bool collideOnBothSides() const
    {
        return m_controller->collideOnBothSides() ;
    }

    short knowsNormalSign( bool atPreviousStep, const ElasticStrand& strand, unsigned vertex ) const
    {
        return m_controller->knowsNormalSign( atPreviousStep, m_faceIndex, strand.getGlobalIndex(), vertex ) ;
    }

    void setNormalSign( short sign, float offset, const ElasticStrand& strand, unsigned vertex ) const
    {
        m_controller->setNormalSign( sign, offset, m_faceIndex, strand.getGlobalIndex(), vertex ) ;
    }

    const TriangularFace& getFace() const
    {
        return m_face;
    }

    //! Guarranteed to start with a zero and end with two zeros
    uint32_t uniqueId() const
    {
        //Not the real Id, but will do as long as it is unique
        // &m_face is aligned on 4 bytes, so
        // we can afford to lose the 4 least significant bits of the address
        // Since we will zero-out the last two, we can still shift by two bits

        return ( ( ( (size_t) &m_face ) >> 2 ) & 0x7FFFFFFCul );
    }

protected:
    virtual void print( std::ostream& os ) const;

    const unsigned m_faceIndex ;
    const TriangularFace& m_face;
    TriMeshController* m_controller;
};

class ElementProxyBBoxFunctor
{
public:
    ElementProxyBBoxFunctor( std::vector<ElementProxy*>& elements ) :
            m_elements( elements )
    {
        for( unsigned i = 0 ; i < m_elements.size() ; ++i )
        {
            if( m_elements[i]->getBoundingBox().isValid() )
            {
                m_valid.push_back( i );
            }
        }
        for( unsigned i = 0 ; i < m_valid.size() ; ++i )
        {
            swap( i, m_valid[i] ) ;
        }
    }

    unsigned int size() const
    {
        return ( unsigned int ) m_valid.size();
    }

    const BBoxType &operator[]( const unsigned int i ) const
    {
        return m_elements[i]->getBoundingBox();
    }

    void swap( unsigned int i, unsigned int j )
    {
        std::swap( m_elements[i], m_elements[j] );
    }

private:
    std::vector<ElementProxy*>& m_elements;
    std::vector<unsigned> m_valid;

};

class ElementProxySortedAABBFunctor
{
public:
   ElementProxySortedAABBFunctor( const std::vector< ElementProxy* > &proxies )
   {
        for( unsigned i = 0 ; i < proxies.size() ; ++ i )
        {
            if( dynamic_cast< EdgeProxy* >( proxies[ i ] ) ){
                m_edgeProxies.push_back( proxies[ i ] ) ;
            } 
            else{
                m_otherProxies.push_back( proxies[ i ] ) ;
            }
        }
   }

   ElementProxy* elementProxy( unsigned elementID ) const
   {
       const unsigned n = numEdgeProxies() ;
       return elementID < n
               ? m_edgeProxies[ elementID ]
               : m_otherProxies[ elementID - n ] ;
   }

   void getAABB( unsigned elementID, Vec3 &min, Vec3 &max ) const
   {
       const BBoxType& bbox = elementProxy( elementID )->getBoundingBox() ;
       min = bbox.min.cast< Scalar >() ;
       max = bbox.max.cast< Scalar >() ;
   }

   unsigned numEdgeProxies() const
   { return m_edgeProxies.size(); }

private:
   std::vector< ElementProxy* > m_edgeProxies ;
   std::vector< ElementProxy* > m_otherProxies ;

};

#endif /* ELEMENTPROXY_HH_ */
