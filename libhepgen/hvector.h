/*!
 *  \file hvector.h
 *  \date Created on: Jan 16, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HVECTOR_H_
#define HVECTOR_H_
#include <cmath>
#include <vector>
#include <iostream>

/*!
 * \brief a simple 3-value-double-precision vector.
 *
 * This is somewhat similar to the TVector3 of root but it has only some basic features needed for the hepgen-project.
 * Not taking ROOT's TVector3 is a choice of because of the design-guideline: NO DEPENDENCIES
 *
 */
class HVector3
{
public:
  /*! \brief raw constructor from 3 doubles */
  HVector3();
  HVector3(double x, double y, double z);
  //HVector3(double _x=0.0, double _y=0.0, double _z=0.0);
  //todo: maybe add a constructor for reading std::vector or other stuff.
  /*! \brief copy constructor */
  HVector3(const HVector3& _copy);
  /*! \brief constructor from std::vector */
  HVector3(const std::vector< double >& _in) {
    fromStdVector(_in);
  };
    
    /*! \brief breaks down double precision to float precision - for test use only! */
    void roundToFloat();

    /*! \brief get the component at position _index (has to be 0,1,2 or else it returns nan */
    double at(unsigned int _index);


    /*! \brief prints the vector to cout */
    void print();

    /*! \brief destructor */
    ~HVector3() {};



    //static getters
    /*! \brief Read-Only access to X */
    double X() const {
        return x;
    };

    /*! \brief Read-Only access to Y */
    double Y() const {
        return y;
    };

    /*! \brief Read-Only access to Z */
    double Z() const {
        return z;
    };
    //setters


    /*! \brief Sets a new X */
    void setX(double _x) {
        x=_x;
    };

    /*! \brief Sets a new Y */
    void setY(double _y) {
        y=_y;
    };

    /*! \brief Sets a new Z */
    void setZ(double _z) {
        z=_z;
    };

    /*! \brief Sets a new X,Y,Z */
    void setXYZ(double _x, double _y, double _z);  /*! gives write access to components */




    /*! \brief Linearly extrapolates along this vector from start vector to end vector */
    void extrapNeutral(const HVector3& _startPoint, HVector3& _endPoint, double distanceZ);


    /*! \brief rotates vector around axis with angle phi in radians */
    void rotateAxisAngle(const HVector3& _axis, const double _phi);


    //just add some more math stuff here as needed
    /*! \brief returns a simple dot-product between the 2 vectors*/
    double dotProduct(const HVector3& _in) const;
    /*! \brief returns the vectorial-product (AxB) between the 2 vectors*/
    HVector3 crossProduct(const HVector3& _in) const;
    /*! \brief returns the length */
    const double length() const;

    /*! \brief normalizes the vector to a new length */
    void normalize(double _newlength);

    /*! \brief rotates the X to the Z axis (COMPASS compat) */
    void rotXToZ();


    /*! \brief rotates the x/y plane with angle _angle */
    void rotPhi(double _angle);


    /*! \brief returns a std::vector with the entries */
    const std::vector<double> toStdVector() const;

    /*! \brief initializes from an std-vector */
    void fromStdVector(const std::vector<double>& _in);


    //overloaded operators

    bool operator== (const HVector3& _in);
    HVector3& operator= (const HVector3& _in);
    HVector3 operator- (const HVector3& _in)const;
    HVector3 operator+ (const HVector3& _in)const;
    HVector3 operator* (double _in);
    const double& operator[] (int const& _index)const;
    double& operator[] (int const& _index);

    /* double operator () (int) const; */
    /* inline double operator [] (int) const; */
    
    /* double & operator () (int); */
    /* inline double & operator [] (int); */
    // Set components by index.
    
    /* inline double x()  const; */
    /* inline double y()  const; */
    /* inline double z()  const; */
    /* inline double X()  const; */
    /* inline double Y()  const; */
    /* inline double Z()  const; */
    inline double Px() const;
    inline double Py() const;
    inline double Pz() const;
    inline void SetX(double);
    inline void SetY(double);
    inline void SetZ(double);

    // The angle w.r.t. another 3-vector.
    double Angle(const HVector3 &) const;

   // Cross product.
    inline HVector3 Cross(const HVector3 &) const;
    // Scalar product.
    inline double Dot(const HVector3 &) const;
    void RotateX(double);
    // Rotates the Hep3Vector around the x-axis.
    
    void RotateY(double);
    // Rotates the Hep3Vector around the y-axis.
    
    void RotateZ(double);
    // Rotates the Hep3Vector around the z-axis.
 
    inline HVector3 operator - () const;
    // Unary minus.

    // Unit vector parallel to this.
    HVector3 Unit() const;

   inline double Mag2() const;
   // The magnitude squared (rho^2 in spherical coordinate system).
   double Mag() const { return sqrt(Mag2()); }
   // The magnitude (rho in spherical coordinate system).
 
private:
    double x,y,z;

    friend class HLorentzVector;

};
////////////////////////////////////////////////////////////////////////////////
/// Constructors
inline HVector3::HVector3()
:  x(0.0),  y(0.0),  z(0.0) {}

inline HVector3::HVector3(double xx, double yy, double zz)
:  x(xx),  y(yy),  z(zz) {}

inline HVector3 HVector3::Cross(const HVector3 & p) const {
   return HVector3(y*p.z-p.y*z, z*p.x-p.z*x, x*p.y-p.x*y);
}

inline double HVector3::Dot(const HVector3 & p) const {
   return x*p.x + y*p.y + z*p.z;
}

inline HVector3 HVector3::operator - () const {
   return HVector3(-x, -y, -z);
}

inline double HVector3::Mag2() const { return x*x + y*y + z*z; }

/* inline double & HVector3::operator[] (int i)       { return operator()(i); } */
/* inline double   HVector3::operator[] (int i) const { return operator()(i); } */

/* inline double HVector3::x()  const { return x; } */
/* inline double HVector3::y()  const { return y; } */
/* inline double HVector3::z()  const { return z; } */
/* inline double HVector3::X()  const { return x; } */
/* inline double HVector3::Y()  const { return y; } */
/* inline double HVector3::Z()  const { return z; } */
inline double HVector3::Px() const { return x; }
inline double HVector3::Py() const { return y; }
inline double HVector3::Pz() const { return z; }

inline void HVector3::SetX(double xx) { x  = xx; }
inline void HVector3::SetY(double yy) { y  = yy; }
inline void HVector3::SetZ(double zz) { z  = zz; }

#endif /* HVECTOR_H_ */


