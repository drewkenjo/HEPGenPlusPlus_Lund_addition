/*!
 *  \file hlorentzvector.h
 *  \date Created on: Jan 16, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HLVECTOR_H_
#define HLVECTOR_H_


#include "hvector.h"
#include <iostream>
#include <cstdio>


/*!
 * \brief a simple Lorentz-4-Vector
 *
 * This is somewhat similar to the TLorentzVector of root but it has only some basic features needed for the hepgen-project.
 * Not taking ROOT's TLorentzVector is a choice of because of the design-guideline: NO DEPENDENCIES
 *
 */
class HLorentzVector
{
public:
    /*! \brief raw constructor from 4 doubles, note that e could be e or t, depending on the interpretation of the other 3 doubles */
    HLorentzVector(double _x=0.0, double _y=0.0, double _z=0.0, double _e=0.0);
    /*! \brief  copy constructor */
    HLorentzVector(const HLorentzVector& _copy);
    /*! \brief  constructor from 3-vector of momentum and energy */
    HLorentzVector(HVector3 _vec, double _e);
    ~HLorentzVector() {};

    enum { kX=0, kY=1, kZ=2, kT=3, kNUM_COORDINATES=4, kSIZE=kNUM_COORDINATES };

    /*! \brief returns 3-vector of momentum CONST */
    const HVector3 getVector() const;

    inline HVector3 Vect() const ;
   // Get spatial component.

    HVector3& getVectorRef() {
        return threeVector;
    };

    /*! \brief returns energy */
    double getEnergy() const {
        return energy;
    };
    
    /*! \brief returns the invariant mass */
    double getMass() const{
        return sqrt(getQuare());
    }
    
    /*! \brief rounds this vector to float precision only */
    void roundToFloat();

    /*! \brief gets the square of the vector */
    double getQuare() const;

    /*! \brief rotates the vector from x to z axis */
    void rotXToZ() {
        threeVector.rotXToZ();
    };
    
    
    /*! \brief makes a dot product of two lorentz vectors */
    double dotProduct(const HLorentzVector& _in) const{
      return (-getEnergy()*_in.getEnergy() + getVector().dotProduct(_in.getVector()));
    }

    double operator () (int i) const;
    inline double operator [] (int i) const;
    // Get components by index.
    
    double & operator () (int i);
    inline double & operator [] (int i);
    // Set components by index.


    inline double X() const;
    inline double Y() const;
    inline double Z() const;
    inline double T() const;

    inline void SetX(double a);
    inline void SetY(double a);
    inline void SetZ(double a);
    inline void SetT(double a);
   // Set position and time.

    inline HVector3 BoostVector() const ;
    // Returns the spatial components divided by the time component.
    
    void Boost(double, double, double);
    inline void Boost(const HVector3 &);
    // Lorentz boost.
    
    /*! \brief prints vector to cout */
    void print();

    //setter
    /*! \brief  sets 3-vector */
    void setVector(const HVector3& _in);
    /*! \brief sets the energy */
    void setEnergy(double _e) {
        energy = _e;
    };

    /*! \brief  sets the vector-properties from spherical coordinates of the momentum and the energy */
    void setLVectorAngular(double _P, double _theta, double _phi,double _E);


    /*! \brief gets transverse momentum wrs z-axis*/
    double getPtrans();

    /*! \brief Boosts the vector */
    void boost(double _u, const HLorentzVector& _ps);
    /*! \brief lorenf reimpls from fortran */
    void lorenf(double _u, const HLorentzVector& _ps);


    //overloaded operators

    bool operator== (const HLorentzVector& _in);
    HLorentzVector& operator= (const HLorentzVector& _in);

    HLorentzVector operator- (const HLorentzVector & _in) const;
    HLorentzVector operator+ (const HLorentzVector & _in) const;
    






private:
    HVector3 threeVector;
    double energy;

};

inline HVector3 HLorentzVector::Vect() const { return threeVector; }

inline HVector3 HLorentzVector::BoostVector() const {
   return HVector3(X()/T(), Y()/T(), Z()/T());
}
inline void HLorentzVector::Boost(const HVector3 & b) {
   Boost(b.X(), b.Y(), b.Z());
}
inline double_t HLorentzVector::X() const { return threeVector.X(); }
inline double_t HLorentzVector::Y() const { return threeVector.Y(); }
inline double_t HLorentzVector::Z() const { return threeVector.Z(); }
inline double_t HLorentzVector::T() const { return energy; }

inline void HLorentzVector::SetX(double a) { threeVector.SetX(a); }
inline void HLorentzVector::SetY(double a) { threeVector.SetY(a); }
inline void HLorentzVector::SetZ(double a) { threeVector.SetZ(a); }
inline void HLorentzVector::SetT(double a) { energy = a; }


inline double & HLorentzVector::operator [] (int i)       { return (*this)(i); }
inline double   HLorentzVector::operator [] (int i) const { return (*this)(i); }

inline double HLorentzVector::operator () (int i) const
{
   //dereferencing operator const
   switch(i) {
      case kX:
	 return threeVector.X();
      case kY:
         return threeVector.Y();
      case kZ:
         return threeVector.Z();
      case kT:
         return energy;
      default:
         return energy;
	//         Error("operator()()", "bad index (%d) returning 0",i);
   }
   return 0.;
}

inline double & HLorentzVector::operator () (int i)
{
   //dereferencing operator
   switch(i) {
      case kX:
	return threeVector.x;
      case kY:
	return threeVector.y;
      case kZ:
	return threeVector.z;
      case kT:
         return energy;
      default:
         return energy;
	 //	Error("operator()()", "bad index (%d) returning &energy",i);
   }
   return energy;
}

#endif /* HLVECTOR_H_ */

