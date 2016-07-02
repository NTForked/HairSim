#ifndef COLLISIONDETECTOR_HH_
#define COLLISIONDETECTOR_HH_

#include <vector>
#include <list>
#include <map>

#include "CollisionUtils/BVH.hh"

#include "CollisionUtils/SpatialHashMapFwd.hh"
#include "TwistEdgeHandler.h"

class ElementProxy;
class EdgeProxy;
class FaceProxy;
class Collision;
class ElementProxySortedAABBFunctor;

class CollisionDetector
{
public:
    CollisionDetector( std::vector< ElementProxy* >& elements );
    virtual ~CollisionDetector();

    void buildBVH( bool statique = false );
    void updateBoundingBoxes();
    void updateBoundingBox( BVHNodeType& node );
    void findCollisions( bool ignoreStrandStrand = false );
    void clear();

    bool empty()
    { return m_collisions.empty(); }

    const std::list< Collision* >& getCollisions() const
    { return m_collisions; }

    std::list< Collision* >& getCollisions()
    { return m_collisions; }

    static void setMaxSizeForElementBBox( double s )
    { s_maxSizeForElementBBox = s; }

    TwistEdgeHandler* m_proxyHistory;
    std::vector<ElementProxy*> m_elementProxies;

protected:

    void computeCollisions( const BVHNodeType& node_a, const BVHNodeType& node_b );
    bool appendCollision( ElementProxy* elem_a, ElementProxy* elem_b );
    bool appendCollision( EdgeProxy* edge_a, EdgeProxy* edge_b );
    bool appendCollision( EdgeProxy* edge_a, const FaceProxy* triangle_b );
    void filterWithSpatialHashMap( const Scalar largestBBox );

    BVH m_bvh;
    std::list< Collision* > m_collisions;
    bool m_ignoreStrandStrand;

    static Scalar s_maxSizeForElementBBox;
    ElementProxySortedAABBFunctor* m_sortedAABBFunctor;
    typedef SpatialHashMap<const ElementProxySortedAABBFunctor, unsigned, true> SpatialHashMapT;
    SpatialHashMapT* m_hashMap;

};

#endif /* COLLISIONDETECTOR_HH_ */
