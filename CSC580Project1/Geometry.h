// Geometry.h
// NKU CSC 480/580 - Kirby 
// -------------------------------------------------------------------------- 
// Basic classes for inhomogeneous points and vectors. Spring 2006 version.
//
// These are in a "csc480" namespace (in case you find yourself using a 
// code with types that have the same name.)
//
// Maintains a clear point-vector distinction (though asVector() and asPoint()
// methods are available to treat one type as another if absolutely necessary).
// -------------------------------------------------------------------------- 

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <cmath>
#include <cassert>
#include <iostream>

namespace csc480 { 

const double PI= 3.14159265359 ;
inline double toRadians( double degrees ) { return PI * degrees / 180.0 ; }
inline double toDegrees( double radians ) { return 180.0 * radians / PI ; } 


// ---------------------------------------------------------------------------
// Vector3D
// --------
//
// A Vector3 object has three components, known as v[0],v[1],v[2] or,
// equivalently, as v.dx, v.dy, v.dz.
// 
// Vector3 objects support scaling, addition, norm, and dot and cross products.
// An asArray() member is provided for use with gl*v() functions.
// Stream I/O overloads of >> and << are available too. 
// ---------------------------------------------------------------------------

struct Point3 ; // forward declaration

struct Vector3
{
    double dx,dy,dz ;
    
    explicit Vector3( double dx=0, double dy=0, double dz=0 ) : dx(dx), dy(dy), dz(dz) {}

    void operator+=( const Vector3& v )     { dx+= v.dx ;  dy+= v.dy ;  dz+=v.dz ; }
    void operator-=( const Vector3& v )     { dx-= v.dx ;  dy-= v.dy ;  dz-=v.dz ; }

    void operator*=( double s )             { dx*= s ;  dy*= s ;  dz*= s ; }
    void operator/=( double s )             { dx/= s ;  dy/= s ;  dz/= s ; }

    double norm() const                     { return sqrt( dx*dx + dy*dy + dz*dz ) ; } 

    const double* asArray() const           { return reinterpret_cast<const double*>(  this ) ; }   
    const Point3& asPoint() const           { return reinterpret_cast<const Point3&>( *this ) ; }

    double* asArray()                       { return reinterpret_cast<double*>(  this ) ; }   
    Point3& asPoint()                       { return reinterpret_cast<Point3&>( *this ) ; }

    const double& operator[]( int i ) const { return asArray()[i] ; } 
    double& operator[]( int i )             { return asArray()[i] ; } 

} ;

inline Vector3 operator-( const Vector3& v )                     // - vector
    { return Vector3( -v.dx, -v.dy, -v.dz ) ; }

inline Vector3 operator+( const Vector3& v, const Vector3& w )   // vector + vector
    { return Vector3( v.dx+w.dx, v.dy+w.dy,  v.dz+w.dz ) ; }

inline Vector3 operator-( const Vector3& v, const Vector3& w )   // vector - vector
    { return Vector3( v.dx-w.dx, v.dy-w.dy,  v.dz-w.dz ) ; }

inline Vector3 operator*( double s, const Vector3& v )           // number * vector
    { return Vector3( v.dx*s, v.dy*s, v.dz*s ) ; }

inline Vector3 operator*( const Vector3& v, double s )           // vector * number
    { return Vector3( v.dx*s, v.dy*s, v.dz*s) ; }

inline Vector3 operator/( const Vector3& v, double s )           // vector / number
    { return Vector3( v.dx/s, v.dy/s, v.dz/s) ; }

inline double operator*( const Vector3& v, const Vector3& w )    // dot product 
    { return v.dx*w.dx + v.dy*w.dy + v.dz*w.dz ; }

inline Vector3 operator^( const Vector3& v, const Vector3& w )   // cross product 
    { return Vector3( v.dy*w.dz - v.dz*w.dy,  
                      v.dz*w.dx - v.dx*w.dz, 
                      v.dx*w.dy - v.dy*w.dx ) ; }

inline Vector3 normalize( const Vector3& v )                     // normalized vector
    { return v / v.norm() ; } 

inline std::ostream& operator<<( std::ostream& out, const Vector3& v ) // out << vector
    { out << "[ " << v.dx << " " << v.dy << " " << v.dz << " ]" ; return out ; }

inline std::istream& operator>>( std::istream& in, Vector3& v )        // in >> vector
    { in >> v.dx >> v.dy >> v.dz ; return in ; }




// ---------------------------------------------------------------------------
// Point3
// ------
//
// A Point3 object has three (Cartesian) coordinates: x,y,z.
// The usual point-vector arithmetic is provided.
// An asArray() member is provided for use with gl*v() functions.
//
// ---------------------------------------------------------------------------


struct Point3
{
    double x,y,z ;

    explicit Point3( double x=0, double y=0, double z=0 ) : x(x), y(y), z(z) {}

    void operator+=( const Vector3& v ) { x+= v.dx ;  y+= v.dy ;  z+=v.dz ; }  
    void operator-=( const Vector3& v ) { x-= v.dx ;  y-= v.dy ;  z-=v.dz ; }  

    const double* asArray() const       { return reinterpret_cast<const double*>( this ) ; } 
    const Vector3& asVector() const     { return reinterpret_cast<const Vector3&>( *this ) ; }

    double*  asArray()                  { return reinterpret_cast<double*>(  this ) ; }   
    Vector3& asVector()                 { return reinterpret_cast<Vector3&>( *this ) ; }
} ;


inline Point3 operator+( const Point3& p, const Vector3& v )  // point + vector
    { return Point3( p.x+v.dx, p.y+v.dy,  p.z+v.dz ) ; }

inline Point3 operator+( const Vector3& v, const Point3& p )  // vector + point
    { return p + v ; }

inline Vector3 operator-( const Point3& a, const Point3& b )  // point - point
    { return Vector3( a.x-b.x, a.y-b.y,  a.z-b.z ) ; }

inline double distance( const Point3& a, const Point3& b )    // distance( point, point )
    { return (a - b).norm() ; }

inline std::ostream& operator<<( std::ostream& out,  const Point3& v )  // out << point
    { out << "( " << v.x << ", " << v.y << ", " << v.z << " )" ; return out ; }

inline std::istream& operator>>( std::istream& in, Point3& v )          // in >> point
    { in >> v.x >> v.y >> v.z ; return in ; }


} // namespace
#endif