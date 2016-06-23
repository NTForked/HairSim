/*
 * \file BendingProducts.cc
 *
 *  Created on: 22/10/2012
 *     \author Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "BendingProducts.hh"
#include "../Core/ElasticStrandUtils.hh"

namespace strandsim
{

void BendingProducts::compute()
{
    m_value.resize( m_size );
    const Mat2x& bendingMatrix = m_bendingMatrixBase.get();
    const GradKArrayType& gradKappas = m_gradKappas.get();

    for ( IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx )
    {
        symBProduct<11>( m_value[vtx], bendingMatrix, gradKappas[vtx] );
    }

    setDependentsDirty();
}

}

