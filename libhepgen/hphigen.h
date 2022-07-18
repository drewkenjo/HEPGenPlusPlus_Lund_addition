/*!
 *  \file hphigen.h
 *  \date Created on: Feb 26, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"


#ifndef HPHIGEN_H_
#define HPHIGEN_H_




/*!
 * \brief Implementation of a Phi-Generator
 *
 */
class HPhysicsGenPHI : public HPhysicsGen
{


public:
    HPhysicsGenPHI() {};
    HPhysicsGenPHI(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenPHI();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    static void calcWeights(hWeightInterface* myInt,double& WEIGHTRET);

    double phi_vm(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  pr_r, int i4);
    double PHI_VM(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, int i4);
    double PHI_PI(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, int i4);

    double sdme(double cos_t, double phi_d, double phi_m, int ii);
    double cos_theta_vm(HLorentzVector LzVecOut1, HLorentzVector LzVecOut2, HLorentzVector pr_r,   int i4);

     void generatePolarisation();

    void setUserVars();
    void setParameters();

    int polarization;

    double weight;
    double weight_sdme;
    double thetapi;
    double costhetapi;
    double phipi;
    double x,y,z;

    int eventcount;
    int resetcount;


};



























#endif





