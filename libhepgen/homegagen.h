/*!
 *  \file homegagen.h
 *  \date Created on: July 27, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"
#include "hrotmat.h"

#ifndef HPhysicsGenOMEGA_H_
#define HPhysicsGenOMEGA_H_




/*!
 * \brief Implementation of a omega-Generator
 *
 */
class HPhysicsGenOMEGA : public HPhysicsGen
{


public:
    HPhysicsGenOMEGA() {};
    HPhysicsGenOMEGA(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenOMEGA();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    /*! \brief static version of calcweights for usage in the hepgen_in_phast or other features
     * Please note: This function calls hreweightKine::getFluxCompensator and therefore needs myInt to be set with:
     * Qsq,tprime, y,nu,slpin -!! it will not error if they are not set correctly!
     */
    static void calcWeights(hWeightInterface* myInt,double& WEIGHTRET);

    void generatePolarisation();

    void setUserVars();
    void setParameters();
    
    double phi_vm(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  LzVecOut3, HLorentzVector  pr_r, int i4);
    double PHI_VM(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  LzVecOut3, int i4);
    double PHI_PI(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  LzVecOut3, int i4);

    double sdme(double cos_t, double phi_d, double phi_m, int ii);
    double cos_theta_vm(HLorentzVector LzVecOut1, HLorentzVector LzVecOut2, HLorentzVector LzVecOut3, HLorentzVector pr_r,   int i4);

    int polarization;

    double weight;
    double weight_sdme;
    double thetapi;
    double costhetapi;
    double theta_n;
    double costheta_n;
    double phipi;
    double phi_n;
    double x,y,z;
    
    double xD,yD;




    int eventcount;
    int resetcount;
private:
    /*! \brief changes the momentum components cyclish style */
    void changeVector(HVector3& _in);
    HRotMat rotationMatrix;
    vector<HParticle*> decayParticles;

};















#endif





