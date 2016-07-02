#ifndef S_COLOR_H
#define S_COLOR_H

#include "../Utils/Definitions.h"

class Color
{
public:

EIGEN_MAKE_ALIGNED_OPERATOR_NEW    ;

    typedef double Channel;

    explicit Color( const Channel& r = 1, const Channel& g = 1, const Channel& b = 1,
            const Channel& alpha = 1 )
    : m_color(r, g, b, alpha)
    {}

    Color(int r, int g, int b, int alpha = 255)
    : m_color( r / 255.0, g / 255.0, b / 255.0, alpha / 255.0 )
    {}

    Color(const Color& color)
    : m_color(color.m_color)
    {}

    Color& operator= (const Color& color)
    {
        m_color = color.m_color;
        return *this;
    }

    const Channel* data() const
    { return m_color.data(); }

    Channel* data()
    { return m_color.data(); }

    int size() const
    { return (int) m_color.size(); }

    Color inverse() const
    {
        return Color( 1.0 - m_color(0), 1.0 - m_color(1), 1.0 - m_color(2), m_color(3) );
    }

protected:

    Eigen::Matrix<Channel, 4, 1> m_color;
};

#endif // COLOR_H
