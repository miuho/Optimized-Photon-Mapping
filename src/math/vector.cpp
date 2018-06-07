/**
 * @file vector.cpp
 * @brief Vector classes.
 *
 * @author Eric Butler (edbutler)
 */

#include "math/vector.hpp"

namespace _462 {

std::ostream& operator<<( std::ostream& os, const Vector2& v )
{
    return os << '(' << v.x << ',' << v.y << ')';
}

std::ostream& operator<<( std::ostream& os, const Vector3& v )
{
    return os << '(' << v.x << ',' << v.y << ',' << v.z << ')';
}

std::ostream& operator<<( std::ostream& os, const Vector4& v )
{
    return os << '(' << v.x << ',' << v.y << ',' << v.z << ',' << v.w << ')';
}

} /* _462 */

