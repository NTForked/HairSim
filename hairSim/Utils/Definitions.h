#ifndef S_DEFINITIONS_H
#define S_DEFINITIONS_H

#include <iostream>

#ifndef STRANDSIM_INCLUDE_VANILLA_EIGEN
#define EIGEN_VECTOR_IO_FORMAT Eigen::IOFormat(8, Eigen::DontAlignCols, ", ", ", ", "", "", "{ ", " }")
#define EIGEN_MATRIX_IO_FORMAT Eigen::IOFormat(8, 0, ", ", "\n", "{ ", " }", "{ ", " }")
#undef EIGEN_DEFAULT_IO_FORMAT // < To silence some warnings about redefining
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_VECTOR_IO_FORMAT
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#undef EIGEN_INITIALIZE_MATRICES_BY_ZERO // < To silence some warnings about redefining
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif // STRANDSIM_INCLUDE_VANILLA_EIGEN

#undef Success // Conflicts with Eigen
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Dense>

namespace Eigen
{
    template<typename _Scalar, int _Options, typename _Index>
    class SparseMatrix;
}

typedef double Scalar; ///< the scalar type
typedef uint16_t IndexType; ///< large unsigned int for IDs

typedef Eigen::Matrix<Scalar, 2, 1> Vec2x; ///< 2d scalar vector
typedef Eigen::Matrix<Scalar, 3, 1> Vec3x; ///< 3d scalar vector
typedef Eigen::Matrix<Scalar, 3, 1> Vec3; ///< 3d scalar vector
typedef Eigen::Matrix<Scalar, 4, 1> Vec4x; ///< 4d scalar vector
typedef Eigen::Matrix<Scalar, 11, 1> Vec11x; ///< 11d scalar vector (stencil for local forces)
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VecXx; ///< arbitrary dimension scalar vector
typedef Eigen::Matrix<unsigned, Eigen::Dynamic, 1> VecXu; ///< arbitrary dimension unsigned vector

typedef std::vector<Vec2x, Eigen::aligned_allocator<Vec2x> > Vec2xArray; ///< an array of 2d scalar vectors
typedef std::vector<Vec3x> Vec3xArray;
typedef std::vector<Vec11x, Eigen::aligned_allocator<Vec11x> > Vec11xArray; ///< an array of 11d scalar vectors

typedef Eigen::Matrix<float, 3, 1> Vec3f;
typedef Eigen::Matrix<float, 4, 1> Vec4f;
typedef Eigen::Matrix<float, 6, 1> Vec6f;
typedef std::vector<Vec3f, Eigen::aligned_allocator<Vec3f> > Vec3fArray; ///< an array of 3d float vectors
typedef std::vector<Vec4f, Eigen::aligned_allocator<Vec4f> > Vec4fArray; ///< an array of 4d float vectors

typedef Eigen::Matrix<double, 2, 1> Vec2d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > Vec2dArray; ///< an array of 2d double vectors
typedef std::vector<Vec3d > Vec3dArray; ///< an array of 3d double vectors
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VecXd; ///< arbitrary dimension scalar vector

typedef Eigen::Matrix<Scalar, 2, 2> Mat2x; ///< 2x2 scalar matrix
typedef Eigen::Matrix<Scalar, 3, 3> Mat3x; ///< 3x3 scalar matrix
typedef Eigen::Matrix<Scalar, 4, 4> Mat4x; ///< 4x4 scalar matrix
typedef Eigen::Matrix<Scalar, 6, 6> Mat6x; ///< 4x4 scalar matrix
typedef Eigen::Matrix<Scalar, 11, 11> Mat11x; ///< 11x11 scalar matrix (stencil for local forces)
typedef std::vector<Mat11x, Eigen::aligned_allocator<Mat11x> > Mat11xArray; ///< an array of 11d scalar matrices
typedef std::pair<Mat11x, Mat11x> Mat11xPair;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatXx; ///< arbitrary dimension scalar matrix

typedef Eigen::Quaternion<Scalar> Quaternion;

typedef Eigen::SparseMatrix< Scalar, Eigen::ColMajor, int > SparseMatx ;
typedef Eigen::SparseMatrix< Scalar, Eigen::RowMajor, int > SparseRowMatx ;

typedef std::vector<int> StrandSet;

template<typename ScalarT>
ScalarT EIGEN_STRONG_INLINE SMALL_NUMBER()
{
    return std::numeric_limits< ScalarT >::epsilon() ;
}

template<>
float EIGEN_STRONG_INLINE SMALL_NUMBER<float>()
{
    return 1e-6;
}

template<>
double EIGEN_STRONG_INLINE SMALL_NUMBER<double>()
{
    return 1e-12;
}

EIGEN_STRONG_INLINE Scalar square( const Scalar x )
{
    return x * x;
}

EIGEN_STRONG_INLINE Scalar cube( const Scalar x )
{
    return x * x * x;
}

template<typename ComparableT>
EIGEN_STRONG_INLINE ComparableT clamp( const ComparableT x, const ComparableT l, const ComparableT u )
{
    return ( x > u ) ? u : ( ( x > l ) ? x : l );
}

template<typename ScalarT>
EIGEN_STRONG_INLINE bool isSmall( ScalarT x )
{
  return fabs( x ) < SMALL_NUMBER<ScalarT>();
}

template<typename NormableT>
EIGEN_STRONG_INLINE bool isClose( const NormableT& x1, const NormableT& x2 )
{
    return isSmall( ( x1 - x2 ).norm() );
}

template<typename NormableT>
EIGEN_STRONG_INLINE bool isApproxUnit( const NormableT& x )
{
    return isSmall( x.squaredNorm() - 1 );
}

namespace std
{
    template < typename Derived >
    inline void swap ( Eigen::DenseBase< Derived >& a, Eigen::DenseBase< Derived >& b )
    {
        a.swap( b ) ;
    }

    template < typename Derived >
    inline void swap ( pair< Eigen::DenseBase< Derived >, Eigen::DenseBase< Derived > >& a,
                       pair< Eigen::DenseBase< Derived >, Eigen::DenseBase< Derived > >& b )
    {
        a.first.swap( b.first ) ;
        a.second.swap( b.second ) ;
    }
}


#endif 
