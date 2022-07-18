#include "homegagen.h"
#include "reweightKine.h"
#define DEBUG 0
void HPhysicsGenOMEGA::generateEvent()
{

    bool eventOK = false;
    //HPhysicsGen::generateEvent();

    while (!eventOK)
    {
        event->getStruct()->listOfParticles.resize(11);
        event->getStruct()->recoil.setParticleAuxFlag(1);

        //generate Nu
        if(!generateNu())
            continue;

        //this will set the beam parameters in the future when beamfilereading is implemented
        if (!setBeam())
            continue;

        //generate Q^2
        if (!generateQSQ())
            continue;

        //meson mass to rho and a bit randomness
        if (!generateMesonMass())
            continue;


        if (!generatePhiGamma())
            continue;
        //maybe this will be implemented in the futute
        if (!generateElastic())
            continue;


        //generate the smearing
        if (!generateSmearing())
            continue;
        //generate the mandelstam t
        if (!generatet())
            continue;
        //generate the gamma kinematics
        if (!generateOutgoingParticle())
            continue;
        //calc the pt and transform outgoing vectors to lab system
        event->calcPT();


        generatePolarisation();

        generateDecay();

        if (!calculatePhir())
            continue;

        //here we need to get
        calcWeights();

        if (doDD)
            ddGotoLab();
        eventOK = true;

    }
    event->rotXToZ();
    if (enabledBeamFile)
        event->rotateEventToBeam(beamEntry.momentum.getVector());

    setParameters();
     //AB    setUserVars();
    if (!doDD){
	weightCounter += weight;
// 	cout << "added weight: " << weight <<" totalling to " << weightCounter << endl;
    }
    generateListOfParticles();

////////////////////////////////////////////////////////////////////////////////////////////////////////
//   SDME FROM WITOLD
////////////////////////////////////////////////////////////////////////////////////////////////////////

    HLorentzVector ebeamlepton = event->getBeam().getVector();
    HLorentzVector escatlepton = event->getScat().getVector();
    HLorentzVector outPart1    = event->getOutPart2_Lab().getVector();
    HLorentzVector outPart2    = event->getOutPart3_Lab().getVector();
    HLorentzVector outPart3    = event->getOutPart4_Lab().getVector();
    HLorentzVector recoilP     = event->getRecoil().getVector();
    // FOR TEST:
    // HLorentzVector ebeamlepton( 0.000e+00,	 0.000e+00,	 1.600e+02,	 1.600e+02);
    // HLorentzVector escatlepton( 2.345e+00,	 1.422e-01,	 9.295e+01,	 9.298e+01);
    // HLorentzVector outPart1(   -1.851e-01,	-2.685e-01,	 9.290e+00,	 9.297e+00);
    // HLorentzVector outPart2(   -2.396e+00,	-2.119e-02,	 5.762e+01,	 5.767e+01);
    // HLorentzVector recoilP(     2.357e-01,	 1.475e-01,	 1.362e-01,	 9.880e-01);
  
    double big_phi = PHI_VM(ebeamlepton,escatlepton,outPart1,outPart2,outPart3, 4);   
    double sml_phi = phi_vm(ebeamlepton,escatlepton,outPart1,outPart2,outPart3,recoilP, 4);
    double cos_the = cos_theta_vm(outPart1,outPart2,outPart3,recoilP, 4);

    costheta_n = cos_the;
    phi_n      = sml_phi;

    int beam_charge=paramMan->getclept(); 
    weight_sdme = 1;
    if( paramMan->getSDME() == 1 ) weight_sdme    =  sdme(cos_the, big_phi, sml_phi, beam_charge); 

    setUserVars();

    if (paramMan->getKeyContents( "ENABLE_DEBUG" ).at ( 1 ) == "1" || paramMan->getDEBUG()==1 )
    {
   
      cout << endl; 
      cout << endl; 
      cout.setf(ios::fixed);
      cout <<"SDME flag         : " << paramMan->getSDME() << "\n" << endl; 
      cout <<"SDME weights         : " << endl ;
      cout <<" type     "   << setprecision (0) << setw( 8)<<event->getRecoil().getParticleType()
	   <<" recop x  "   << setprecision (5) << setw( 8)<<recoilP.X()
	   <<" recop y  "   << setprecision (5) << setw( 8)<<recoilP.Y()
	   <<" recop z  "   << setprecision (5) << setw( 8)<<recoilP.Z()
	   <<" recop E  "   << setprecision (5) << setw( 8)<<recoilP.T()
	   <<" recop M  "   << setprecision (5) << setw( 8)<<recoilP.getMass()
	   << endl;
      cout <<" type 1   "   << setprecision (0) << setw( 8)<<event->getOutPart2_Lab().getParticleType()
	   <<" out 1 x  "   << setprecision (5) << setw( 8)<<outPart1.X()
	   <<" out 1 y  "   << setprecision (5) << setw( 8)<<outPart1.Y()
	   <<" out 1 z  "   << setprecision (5) << setw( 8)<<outPart1.Z()
	   <<" out 1 E  "   << setprecision (5) << setw( 8)<<outPart1.T()
	   <<" out 1 M  "   << setprecision (5) << setw( 8)<<outPart1.getMass()
	   << endl;
      cout <<" type 2   "   << setprecision (0) << setw( 8)<<event->getOutPart3_Lab().getParticleType()
	   <<" out 2 x  "   << setprecision (5) << setw( 8)<<outPart2.X()
	   <<" out 2 y  "   << setprecision (5) << setw( 8)<<outPart2.Y()
	   <<" out 2 z  "   << setprecision (5) << setw( 8)<<outPart2.Z()
	   <<" out 2 E  "   << setprecision (5) << setw( 8)<<outPart2.T()
	   <<" out 2 M  "   << setprecision (5) << setw( 8)<<outPart2.getMass()
	   << endl;
      cout <<" type 3   "   << setprecision (0) << setw( 8)<<event->getOutPart4_Lab().getParticleType()
	   <<" out 3 x  "   << setprecision (5) << setw( 8)<<outPart3.X()
	   <<" out 3 y  "   << setprecision (5) << setw( 8)<<outPart3.Y()
	   <<" out 3 z  "   << setprecision (5) << setw( 8)<<outPart3.Z()
	   <<" out 3 E  "   << setprecision (5) << setw( 8)<<outPart3.T()
	   <<" out 3 M  "   << setprecision (5) << setw( 8)<<outPart3.getMass()
	   << endl;
      
      cout << endl; 
      cout <<" COS(th)  "   << setprecision (5) << setw( 8)<<cos_the
	   <<" PHI      "   << setprecision (5) << setw( 8)<<big_phi    
	   <<" phi      "   << setprecision (5) << setw( 8)<<sml_phi         
	   <<" SDME w   "   << setprecision (5) << setw( 8)<<weight_sdme
	   <<" weight   "   << setprecision (5) << setw( 8)<<weight      
	   << endl;
      event->getStruct()->listOfParticles.at(0)->printDebugHeader();
      for (int i =0; i < event->getStruct()->listOfParticles.size(); i ++ )
        {
	  cout << i+1 << " ";
	  event->getStruct()->listOfParticles.at(i)->printDebug();
        }
    }
}

HPhysicsGenOMEGA::HPhysicsGenOMEGA(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "w-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    //bookMan->addHBook1D(&phir,100,0,2*M_PI,"Phi_r");

    event->getStruct()->type = hepconst::OMEGA;



    eventcount = 0;
    resetcount = 0;


    //build the shuffle list for the three-body-decay
    decayParticles.clear();
    decayParticles.push_back(&event->getStruct()->outPart2);//pi+
    decayParticles.push_back(&event->getStruct()->outPart3);//pi-
    decayParticles.push_back(&event->getStruct()->outPart4);//pi0

    decayParticles.at(0)->setParticleType(hepconst::typePiPlus);
    decayParticles.at(1)->setParticleType(hepconst::typePiMinus);
    decayParticles.at(2)->setParticleType(hepconst::typePi0);


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle); //0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle); //1
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt); //2
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle); //3
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab); //4 omega
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab); //5 pi+
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab); //6 pi-
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart4_Lab); //7 pi0
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart5_Lab); //8 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart6_Lab); //9 gamma
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil); //10





    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(6)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(7)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(8)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(9)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(10)->setSharpMass(true);


    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(6)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(7)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(8)->setParticleOrigin(8);
    event->getStruct()->listOfParticles.at(9)->setParticleOrigin(8);
    event->getStruct()->listOfParticles.at(10)->setParticleOrigin(2);



    //set the aux-flags of the particles, 3 are dead already 1 is decayed rest is aLivE!
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(6)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(7)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(8)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(9)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(10)->setParticleAuxFlag(1);

    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(11);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(9);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(10)->setParticleDaughter1(0);


    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(8);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(10);
    event->getStruct()->listOfParticles.at(8)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(9)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(10)->setParticleDaughter2(0);




}


void HPhysicsGenOMEGA::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,1,"weight(#omega)");
    bookMan->addHBook1D(&weight_sdme,&event->getStruct()->dummyW,100,-0.01,0.15,"weight SDME (#omega)");

    bookMan->addHBook1D(&phi_out     , &weight, 63     , 0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&theta_out   , &weight,100     , 0,  M_PI,"Theta_out");
    bookMan->addHBook1D(&phir        , &weight, 63     , 0,2*M_PI,"Phi_r");
    bookMan->addHBook1D(&thetapi     , &weight, 32     , 0,  M_PI,"Theta_pi");
    bookMan->addHBook1D(&thetapi     , &weight_sdme, 32, 0,  M_PI,"Theta_pi (SDME)");
    bookMan->addHBook1D(&costhetapi  , &weight,100     ,-1,     1,"Cos(Theta_pi) hepgen w");
    bookMan->addHBook1D(&costhetapi  , &weight_sdme,100,-1,     1,"Cos(Theta_pi) (SDME)");
    bookMan->addHBook1D(&costheta_n  , &weight,100     ,-1,     1,"Cos(Theta_n) hepgen w");
    bookMan->addHBook1D(&costheta_n  , &weight_sdme,100,-1,     1,"Cos(Theta_n) (SDME)");
    bookMan->addHBook1D(&phipi       , &weight, 63     , 0,2*M_PI,"Phi_pi");
    bookMan->addHBook1D(&phipi       , &weight_sdme,63 , 0,2*M_PI,"Phi_pi (SDME)");
    bookMan->addHBook1D(&phi_n       , &weight, 63     , 0,2*M_PI,"Phi_n");
    bookMan->addHBook1D(&phi_n       , &weight_sdme, 63, 0,2*M_PI,"Phi_n (SDME)");

    bookMan->addHBook1D(&m_meson,&weight,200,0.5,1.0,"Mass of #omega meson");
 
    bookMan->addHBook1D(&beta_proton,&weight,3000,0,2,"Beta_Proton_Recoil");
    bookMan->addHBook1D(&theta_proton,&weight,3000,-6,6,"theta proton");
    bookMan->addHBook1D(&x,&weight,3000,-6,6,"X Comp Momentum");
    bookMan->addHBook1D(&y,&weight,3000,-6,6,"Y Comp Momentum");
    bookMan->addHBook1D(&z,&weight,3000,-6,6,"Z Comp Momentum");
    
    bookMan->addHBook2D(&xD,&yD,&event->getStruct()->dummyW,1000,1000,-0.35,-0.05,0.35,0.65,"DalitzPlot");
}





void HPhysicsGenOMEGA::generateDecay()
{

    //generate the particle decay
    //this is a bit problematic here, as w->pi+ pi- pi0
    //three-body-decay that cannot be simply calculated

    //we start like with the rho generator

    /*----- standard vektor meson part ----*/
    // first we generate the angles of the first pion

    double r = myRandom->flat();
    double costh;
    double angle;

    if( paramMan->getSDME() == 0 ) {
      
      if(polarization == 0)
	{
	  if(r < 0.5)
            costh=-1.0*pow(abs(2*r-1.),(1./3.));
	  else if(r > 0.5)
            costh=+1.0*pow(abs(2*r-1.),(1./3.));
	  else
            costh=0.0;
	}
      else
	{
	  angle=acos(abs(2*r-1.));
	  if(r<0.5)
            costh=-2.0*cos((M_PI+angle)/3.);
	  else
            costh=+2.0*cos((M_PI+angle)/3.);
	}
    } else {
      costh= -1.0 + 2* r;
    }
    
    //generate phi now
    r = myRandom->flat();
    phipi=r*2*M_PI;
    // now we have all the angles
    costhetapi = costh;
    thetapi = acos (costh);

    /* ----- special omega part ----- */

    //this is in to not break the numerical compat with hepgen

    r = myRandom->flat();
    HParticle *part1,*part2,*part3;
    


    if (r < 0.33)
    {
        part2= &event->getStruct()->outPart2;//+
        part3= &event->getStruct()->outPart3;//-
        part1= &event->getStruct()->outPart4;//0
    }
    else if (r >= 0.33 &&  r < 0.66)
    {
        part3= &event->getStruct()->outPart2;//+
        part1= &event->getStruct()->outPart3;//-
        part2= &event->getStruct()->outPart4;//0
    }
    else if (r >= 0.66)
    {

        part1= &event->getStruct()->outPart2;//+
        part2= &event->getStruct()->outPart3;//-
        part3= &event->getStruct()->outPart4;//0
    }


    event->getStruct()->outPart2.setSharpMass(true);
    event->getStruct()->outPart3.setSharpMass(true);
    event->getStruct()->outPart4.setSharpMass(true);


    //this leads to segfaults - i really dont know why... this is so fucked up...

    //whatever we just do it the fortran way
//     decayParticles.clear();
//     decayParticles.push_back(part1);
//     decayParticles.push_back(part2);
//     decayParticles.push_back(part3);
//




    //shuffle the outparticles
//     random_shuffle(decayParticles.begin(),decayParticles.end());


    //nota bene: this yields 6 cases whereas hepgen only yields 3 cases
    //as it assumes pi+,pi- indistinguishable! Therefore one should respect this
    double massDiPionMin = part2->getMass()+part3->getMass();

    /*  for (int i=0; i < 3;i++)
        decayParticles.at(i)->setSharpMass(true);
    */
    //maybe this would be more beautiful with using the square-function of the lorentzvectors instead of mass-getters? i dont know!
    double maxEnergy = (pow(hepconst::mOmega,2.) + part1->getMassSq() - pow(massDiPionMin,2.))/(2*hepconst::mOmega);
    double maxMomentum = sqrt(pow(maxEnergy,2.)-part1->getMassSq());


    //now roll for the momentum of particle A
    bool momOkay = false;
    double thetapq,qDiPion,massDiPion,energyDiPion,momA ;
    while (!momOkay)
    {
        //choose momentum of pion A
        double r = myRandom->flat();
        momA = r*maxMomentum;

        //now we can calculate the kinematics of the virtual di-pion
        energyDiPion = hepconst::mOmega - sqrt(part1->getMassSq()+pow(momA,2.));
        massDiPion = sqrt(pow(energyDiPion,2.)-pow(momA,2.));
        //this is the momentum, but we call it q for historical reasons
        qDiPion = sqrt((pow(massDiPion,2.)-pow(massDiPionMin,2.))*(pow(massDiPion,2.)-pow((part2->getMass()-part3->getMass()),2.)))/(2*massDiPion);

        //now we choose cosine angle between pion and di-pion
        r = myRandom->flat();
        double cospq = 2*r-1;
        thetapq = acos(cospq);


        double prob = 9 * pow(momA,2.) * pow(qDiPion,2.) * ( 1 - pow(cospq,2.));


        r = myRandom->flat();
	if ( r * hepconst::probMax < prob)
            momOkay = true;
    }


    HVector3 momentumB = HVector3(qDiPion*cos(thetapq),qDiPion*sin(thetapq),0.);
    part2->getVector().setVector(momentumB);
    part2->getVector().setEnergy(sqrt(part2->getMassSq()+pow(qDiPion,2.)));
    part3->getVector().setVector(HVector3(0.,0.,0.)-momentumB);
    part3->getVector().setEnergy(sqrt(part3->getMassSq()+pow(qDiPion,2.)));



    //here this is already in cm-system
    HVector3 momentumA = HVector3(momA,0.,0.);
    part1->getVector().setVector(momentumA);
    part1->getVector().setEnergy(sqrt(pow(momA,2.)+part1->getMassSq()));
    HParticle diPion;
    diPion.setSharpMass(false);
    HLorentzVector diPionVec(-momA,0,0,energyDiPion);

    //now we rot the diPion-products into the cm system as well by lorentz-boosting along the diPionVec lorentzvector
    part2->getVector().boost(massDiPion,diPionVec);
    part3->getVector().boost(massDiPion,diPionVec);


    //choose random phi angle and rotate all the vectors in it
    r = myRandom->flat();
    double phiToRot = 2 * M_PI * r;
    //make the rotation matrix ready
    rotationMatrix.setFromThetaPhiVector(thetapi,phipi,event->getStruct()->outPart1_Lab.getTVector());


    //rotate the outgoing particles

    part1->getVector().getVectorRef().rotPhi(phiToRot);
    part2->getVector().getVectorRef().rotPhi(phiToRot);
    part3->getVector().getVectorRef().rotPhi(phiToRot);

    changeVector(part1->getVector().getVectorRef());
    changeVector(part2->getVector().getVectorRef());
    changeVector(part3->getVector().getVectorRef());
    
    double momsum = hepconst::mOmega - 2 * hepconst::mpic - hepconst::mpi;
    double tplus = event->getStruct()->outPart2.getVector().getEnergy() - hepconst::mpic;
    double tminus= event->getStruct()->outPart3.getVector().getEnergy() - hepconst::mpic;
    double tzero = event->getStruct()->outPart4.getVector().getEnergy() - hepconst::mpi;
    xD = (tminus - tplus)/(sqrt(3.)*momsum);
    yD = tzero/momsum;
    

//        part2= &event->getStruct()->outPart2;//+
//         part3= &event->getStruct()->outPart3;//-
//         part1= &event->getStruct()->outPart4;//0


    part1->getVector().getVectorRef() = rotationMatrix.rotateVector(part1->getVector().getVectorRef());
    part2->getVector().getVectorRef() = rotationMatrix.rotateVector(part2->getVector().getVectorRef());
    part3->getVector().getVectorRef() = rotationMatrix.rotateVector(part3->getVector().getVectorRef());


    //okay we are almost there, now we go to the lab system
    event->getStruct()->outPart2_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart2.getVector(),event->getStruct()->outPart1_Lab.getVector()));
    event->getStruct()->outPart3_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart3.getVector(),event->getStruct()->outPart1_Lab.getVector()));
    event->getStruct()->outPart4_Lab.setVector(event->goToLabSystem(event->getStruct()->outPart4.getVector(),event->getStruct()->outPart1_Lab.getVector()));

    //finally we still have the pi0 decay to do
    double mom = sqrt(hepconst::w2pi)/2;
    double theta_gamma = acos(myRandom->flat()*2-1);
    double phi_gamma = 2*myRandom->flat()*M_PI;
    HLorentzVector gammaVect;
    gammaVect.setLVectorAngular(mom,theta_gamma,phi_gamma,mom);
    HVector3 myVect = gammaVect.getVector();

    //invert momentum for second gamma
    HLorentzVector gammaVect2;
    gammaVect2.setVector(HVector3(0,0,0)-myVect);
    gammaVect2.setEnergy(mom);

    //bring both to the lab system
    HLorentzVector gammaLab1 = event->goToLabSystem(gammaVect,event->getStruct()->outPart4_Lab.getVector());
    HLorentzVector gammaLab2 = event->goToLabSystem(gammaVect2,event->getStruct()->outPart4_Lab.getVector());

    //set it correct to the particles that get written to the event
    event->getStruct()->outPart5_Lab.setParticleType(22);
    event->getStruct()->outPart6_Lab.setParticleType(22);

    event->getStruct()->outPart5_Lab.setVector(gammaLab1);
    event->getStruct()->outPart6_Lab.setVector(gammaLab2);

}








HPhysicsGenOMEGA::~HPhysicsGenOMEGA()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenOMEGA::changeVector(HVector3& _in)
{
    _in.setZ(_in.Y());
    _in.setY(_in.X());
    _in.setX(0.0);
}


void HPhysicsGenOMEGA::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    event->getStruct()->USERVAR.at(3)=weight_sdme;
   if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;

}



void HPhysicsGenOMEGA::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeOmega;
    paramMan->getStruct()->PARL.at(7)=hepconst::typePiPlus;
    paramMan->getStruct()->PARL.at(8)=hepconst::typePiMinus;
    paramMan->getStruct()->PARL.at(9)=hepconst::typePi0;
}



bool HPhysicsGenOMEGA::generateMesonMass()
{
    eventcount++;

    double r  = myRandom->flat();

    //pretty much a copy of the rho mass generation with different constants
    //maybe put this into the hphysicsgen class?
    m_meson = (0.84-0.73) * r + 0.73;
    double qcms=sqrt((pow(m_meson,2.0)-pow((2*hepconst::mpic+hepconst::mpi),2.0))*(pow(m_meson,2.0)-pow((2*hepconst::mpic-hepconst::mpi),2.0)))/(2*m_meson);
    double grho=hepconst::gOmega * pow((qcms/hepconst::qOmega),3) * (hepconst::mOmega/m_meson);
    //double grhop=hepconst::gRho0 * pow((qcms/hepconst::qCMS0),3) * (2*pow(hepconst::qCMS0,2)/(pow(hepconst::qCMS0,2)+pow(qcms,2)));
    double prob=m_meson*hepconst::mOmega*grho*hepconst::gOmega / ( pow((pow(m_meson,2)-pow(hepconst::mOmega,2)),2)
                +pow(hepconst::mOmega,2)*pow(grho,2));
    r = myRandom->flat();
    //the thing is a relativistic breit wigner distribution ala HERA hep-ex/9507011v2
    //for theory see the old ancient nuovo cimenta paper: http://www-theory.lbl.gov/jdj/NuoCim.pdf
    //check if this mass is good, if not we just do it again!
    if(r > prob)
    {
        resetcount++;
        this->generateMesonMass();
    }
    //now we set the meson mass
    event->getStruct()->m_meson = m_meson;

    //and check for inconsistencies with the exclusive generation
    //if so, we throw it away and start all over
    if (sqrt(event->getWsq()) < hepconst::w2prma + m_meson)
        return false;
    else
    {
        HPhysicsGen::generateMesonMass();
        return true;
    }//this sets the parl correctly

}



void HPhysicsGenOMEGA::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeOmega);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typePiPlus);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typePiMinus);
    event->getStruct()->outPart4_Lab.setParticleType(hepconst::typePi0);
    event->getStruct()->outPart5_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart6_Lab.setParticleType(hepconst::typeGamma);
}





void HPhysicsGenOMEGA::generatePolarisation()
{

    //ratio longitudinal to transversal polarized mesons:
    //this seriously is evil!
    double ratiolt = pow((1.+ event->getQsq() / paramMan->getStruct()->amt2_rho ),2) * paramMan->getStruct()->aksi2 *
                     pow((M_PI/2.*paramMan->getStruct()->aml2_rho/event->getQsq() - paramMan->getStruct()->aml2_rho*sqrt(paramMan->getStruct()->aml2_rho)/(sqrt(event->getQsq())*(event->getQsq()+paramMan->getStruct()->aml2_rho))
                          -paramMan->getStruct()->aml2_rho/event->getQsq()*atan(sqrt(paramMan->getStruct()->aml2_rho)/sqrt(event->getQsq()))),2);


    double fractionl = (event->getStruct()->epsilon + event->getStruct()->delta) * ratiolt / (1+(event->getStruct()->epsilon + event->getStruct()->delta)*ratiolt);

    //roll the dice
    double r = myRandom->flat();

    if (r < fractionl)
        polarization = 0;
    else
        polarization = 1;
    event->getStruct()->epp_polarized_longitudinal=polarization;


}


void HPhysicsGenOMEGA::calcWeights ( hWeightInterface* myInt, double& WEIGHTRET ) {
    //omega cross section is pretty much identical to cross section of Rho0
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    double beta = 1.96;
    double siggam=SIG0RO*pow((QSQ0/myInt->qsq),beta)*exp(-myInt->slpin*myInt->tprim)*myInt->slpin;
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;
    
    // this needs to be multiplied by flux and phasefactor!
    
    WEIGHTRET=siggam;
    WEIGHTRET*= hreweightKine::getFluxCompensator(*myInt);
    

    if (WEIGHTRET < 0)
        WEIGHTRET = 0;

    //theres just a factor of 10 to divide by for omega now:
    WEIGHTRET /= 10.;
    //thats it :)

}


void HPhysicsGenOMEGA::calcWeights()
{
    //double ALF = 0.007299;
    //double AMP = 0.938256;
    double SIG0RO =27.4;
    double QSQ0 = 6.;
    double beta = 1.96;
    double siggam=SIG0RO*pow((QSQ0/event->getQsq()),beta)*exp(-paramMan->getSlopeIncoherent()*event->getStruct()->tprim)*paramMan->getSlopeIncoherent();
// c ---- renormalised to GK prediction at Q2=4 GeV2 and W=10 GeV
// c ---- sigma_L = 97.5 nb, sigma_T = 51.8 nb
    siggam *=2.33;

    double funrho0=event->getStruct()->flux*siggam;
    weight = funrho0 * event->getTotalPhaseFactor();
    // cout << " weights " << funrho0 << " " << event->getTotalPhaseFactor() << " cut " << paramMan->getStruct()->CUT.at(0) <<endl;
    // cout << " weights " << funrho0 << " t " << event->getStruct()->PFt << " nu " <<  event->getStruct()->PFnu << " Q2 " <<  event->getStruct()->PFqsq <<endl;
    
    if (weight < 0)
        weight = 0;

    weight /= 10.;

    //add some histogram-infos here

    theta_proton = acos(event->getRecoil().getTVector().X()/event->getRecoil().getTVector().length());
    z = event->getRecoil().getTVector().Z();
    y = event->getRecoil().getTVector().Y();
    x = event->getRecoil().getTVector().X();


    beta_proton = (event->getRecoil().getTVector().length()/hepconst::w2prma);

}

/////////////////////////////////////////////////////////////////////////////////////////
// phi small  \phi .....  
// mu_in[ ] incaming miuon mu_out outcoming miuon 
// pi_1 , pi_2, pi_3 pro_rec as above ...................  
///////////////////////////////////////////////////////////////////////////////////////
double HPhysicsGenOMEGA::phi_vm( HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  LzVecOut3, HLorentzVector  pr_r, int i4)
{
  const double M_Pr = 0.93827203; // proton mass
  HLorentzVector q,LzVecOut,HCMS ;
  HLorentzVector p(0.0,0.0,0.0,M_Pr);

  q = LzVecMu0 - LzVecMu1;

  LzVecOut= LzVecOut1+LzVecOut2+LzVecOut3;
  HCMS= q + p;
  
  //  HVector3  b = -HCMS.BoostVector();          // boost vector TO HCMS system
  HVector3  b = -LzVecOut.BoostVector();          // boost vector TO VM system
  HLorentzVector lv1p_vm(LzVecOut1); lv1p_vm.Boost(b); // transform p1 4vector to HCMS system
  HLorentzVector lv2p_vm(LzVecOut2); lv2p_vm.Boost(b); // transform p2 4vector to HCMS system
  HLorentzVector q_vm(q);               q_vm.Boost(b); // transform q  4vector to HCMS system
  HLorentzVector pr_rvm(LzVecOut);    pr_rvm.Boost(b); // transform VM 4vector to HCMS
  HLorentzVector pr_rpr(pr_r    );    pr_rpr.Boost(b); // transform p' 4vector to HCMS
  HVector3 llv1_vm =  lv1p_vm.Vect();
  HVector3 llv2_vm =  lv2p_vm.Vect();
  HVector3 enne    =  llv1_vm.Cross(llv2_vm);
  HVector3  qq3    =     q_vm.Vect();
  HVector3  VV  =  pr_rvm.Vect(); // for phi and rho use VM vector
  HVector3  PP  =  pr_rpr.Vect(); // for omega use scattered p vector 
  double ccos=(( qq3.Cross( PP)).Unit()).Dot((PP.Cross(enne)).Unit());            //      for omega change sign of sin   
  double ssin=(((qq3.Cross( PP)).Cross(PP)).Unit()).Dot((enne.Cross(PP)).Unit());
  double phimale=atan2(ssin,ccos);
  return  phimale ;  //poprawione 12*XI*2016
}

////////////////////////////////////////////////////////////////////////////////////////
//  \Phi big    -----  
//  the same def. as mentioned above
/////////////////////////////////////////////////////////////////////////////
double HPhysicsGenOMEGA::PHI_VM(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2,  HLorentzVector  LzVecOut3, int i4)
{
  
  const double M_Pr = 0.93827203; // proton mass
  HLorentzVector q,LzVecOut ;
  HLorentzVector p(0.0,0.0,0.0,M_Pr);
  q= LzVecMu0- LzVecMu1;
  LzVecOut= LzVecOut1+LzVecOut2+LzVecOut3;
  HVector3 lv = LzVecOut. Vect(); 
  HVector3 lv1= LzVecOut1.Vect();
  HVector3 lv2= LzVecOut2.Vect();
  HVector3 lv3= LzVecOut3.Vect();
  HVector3 qq= q.Vect();  
  HVector3 beta, k1xk2,qxv,k1xk2u,qxvu, ksin,qrp3u  ;
  double thet_1,phi_1,cme,k1z;
  
  thet_1=acos(qq[2]/qq.Mag());
  phi_1= atan2(qq[1],qq[0]);
  
  HVector3 km1r = LzVecMu0.Vect(); // lepton
  HVector3 km2r = LzVecMu1.Vect(); // lepron rozpr
  HVector3 qrx = q.Vect();         // q
  HVector3 vmr = lv;               // VM

 //     cout<<" vmr 2 "<< vmr[2] << endl;
  
  km1r.RotateZ(-phi_1); km1r.RotateY(-thet_1);
  km2r.RotateZ(-phi_1); km2r.RotateY(-thet_1);
  qrx.RotateZ(-phi_1); qrx.RotateY(-thet_1);
  vmr.RotateZ(-phi_1); vmr.RotateY(-thet_1);
  
  /*      cout << " km1r "<< km1r[0]<< " " << km1r[1] << endl;
	  cout << " km2r "<< km2r[0]<< " " << km2r[1] << endl;
	  cout << " qrx "<< qrx[0]<< " " << qrx[1]<< " "<< qrx[2] << endl;
	  cout << " vmr "<< vmr[0]<< " " << vmr[1] << endl;
  */     
  HLorentzVector k1r, k2r, qr, vr ;
  k1r[0]=km1r[0]; k1r[1]=km1r[1]; k1r[2]= km1r[2]; k1r[3]= LzVecMu0[3]; // mion beam rot
  k2r[0]=km2r[0]; k2r[1]=km2r[1]; k2r[2]=km2r[2];  k2r[3]=LzVecMu1[3]; // mion scattered
  qr[0]= qrx[0]; qr[1]=qrx[1] ; qr[2]= qrx[2]; qr[3]=q[3];         // rotated 
  vr[0]= vmr[0]; vr[1]=vmr[1]; vr[2]=vmr[2]; vr[3]=LzVecOut[3];        //
  //     cout << " proton "<< p[3] <<endl; 
  HLorentzVector b_cmr =  qr + p; // cm frame
  /*      cout << " k1r "<< k1r[0]<< " " << k1r[1]<<" "<< k1r[2]<<" " << k1r[3] << endl;
	  cout << " k2r "<< k2r[0]<< " " << k2r[1]<<" " << k2r[2]<<" "<< k2r[3] << endl;
	  cout << " qr "<< qr[0]<< " " << qr[1]<<" "<< qr[2]<< " "<< qr[3]  << endl;
	  cout << " vr "<< vr[0]<< " " << vr[1]<<" "<< vr[2]<<" "<< vr[3] << endl;
	  
  */ 
  HVector3 b_cm= -b_cmr.BoostVector();          // df of boost
  HLorentzVector k1rp(k1r); k1rp.Boost(b_cm);  // u0 -beam
  HLorentzVector k2rp(k2r); k2rp.Boost(b_cm);  // u1 scatt.
  HLorentzVector qrp(qr); qrp.Boost(b_cm);     // q_cm
  HLorentzVector vrp(vr); vrp.Boost(b_cm);
  
  //   cout<< " b_cm "<< b_cm[0] << " " << b_cm[1] << " "<< b_cm[2] << endl;   
  
  HVector3 k11=k1rp.Vect();
  HVector3 k12=k2rp.Vect();
  HVector3 qrp3=qrp.Vect();
  HVector3 vrp3=vrp.Vect();
  k1xk2=k11.Cross(k12);
  qxv= qrp3.Cross(vrp3);
  k1xk2u=k1xk2.Unit();
  qxvu=qxv.Unit();
  qrp3u= qrp3.Unit();
  
  double CosPHI,SinPHI,bbb,PHID;
  CosPHI=k1xk2u.Dot(qxvu);
  ksin= qxvu.Cross(k1xk2u);
  SinPHI=ksin.Dot(qrp3u); 
  PHID=atan2(SinPHI,CosPHI);
  
  return PHID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// azym for pion from omega
// TRENTO CONVENTION .....................
/////////////////////////////
double HPhysicsGenOMEGA::PHI_PI(HLorentzVector  LzVecMu0, HLorentzVector LzVecMu1, HLorentzVector LzVecOut1, HLorentzVector  LzVecOut2, HLorentzVector  LzVecOut3, int i4)
{
  //       LzVecOutX --- pion  
  const double M_Pr = 0.93827203; // proton mass
  HLorentzVector  LzVecOut, q,HCMS,LzVecOutX  ;
  HLorentzVector p(0.0,0.0,0.0,M_Pr);
  
  q= LzVecMu0- LzVecMu1;  
  LzVecOut= LzVecOut1+LzVecOut2+ LzVecOut3;
  LzVecOutX= LzVecOut1;
  HCMS= q + p;
  HVector3  b = -HCMS.BoostVector();
  HLorentzVector q_cm(q); q_cm.Boost(b); // transform q  4vector to HCMS system
  HLorentzVector V_cm(LzVecOut); V_cm.Boost(b); // transform VM 
  HLorentzVector PI_cm(LzVecOutX); PI_cm.Boost(b); // transform PIONU 
  HLorentzVector K1_cm(LzVecMu1); K1_cm.Boost(b); // MU1 to HCMS system
  HLorentzVector K_cm(LzVecMu0); K_cm.Boost(b); // mu   4vector to HCMS system
  HVector3 QQ= q_cm.Vect();       
  HVector3 VV= V_cm.Vect();
  HVector3 PIP= PI_cm.Vect(); // pion cm 
  HVector3 KK= K_cm.Vect();
  HVector3 KP= K1_cm.Vect();
  //        double  COSPHI=((QQ.Cross(VV)).Unit()).Dot((KK.Cross(KP)).Unit());
  double  COSPHI=((QQ.Cross(PIP)).Unit()).Dot((KK.Cross(KP)).Unit());
  //         TVector3   SIN1 = (QQ.Cross(VV)).Unit();
  HVector3   SIN1 = (QQ.Cross(PIP)).Unit();
  HVector3   SIN2 = (KK.Cross(KP)).Unit();
  double  SINPHI=(SIN1.Cross(SIN2)).Dot(QQ.Unit()); 
  double  PHID_PI=atan2(SINPHI,COSPHI);
  
  return -PHID_PI; // TRENTO CONVENTION
}
///////////////////////////////////////////////////////////////////////////
//   arguments for sdme  : cos_t cos(\theta); phi_d - big PHI;  phi_m - small phi;
//     epsl is \epsilon  pol_beam pol. of muons   ii - mu beam mu+ =1 mu- =-1
//    helicity def as pol_beam * (-1 ) * ii last line in sdme function

///////////////////////////////////////////////////////////////////////////
// SEQENCE of SDME's :
// USED SET of 
//   NU  SDME        Value         ERROR 
// ==============================================================
  // 1   r^04_00     0.49926       0.50619E-02   0.12679E-04   0.49926    
  // 2   r^04_10     0.56259E-01   0.30369E-02   0.66754E-05   0.56259E-01
  // 3   r^04_1-1   -0.15964E-01   0.37399E-02   0.44353E-04  -0.15964E-01
  // 4   r^1_11     -0.20392E-01   0.42056E-02   0.46044E-04  -0.20392E-01
  // 5   r^1_00     -0.36005E-01   0.87320E-02   0.99864E-04  -0.36005E-01
  // 6   r^1_10     -0.48245E-01   0.39147E-02   0.49147E-04  -0.48245E-01
  // 7   r^1_1-1     0.20074       0.48648E-02   0.63559E-04   0.20074    
  // 8   r^2_10      0.49471E-01   0.38722E-02   0.49080E-04   0.49471E-01
  // 9   r^2_1-1    -0.20010       0.49862E-02   0.65039E-04  -0.20010    
  //10   r^5_11      0.36627E-02   0.21844E-02   0.23683E-04   0.36627E-02
  //11   r^5_00      0.13666       0.46404E-02   0.49318E-04   0.13666    
  //12   r^5_10      0.15592       0.20262E-02   0.25188E-04   0.15592    
  //13   r^5_1-1    -0.60333E-03   0.26790E-02   0.66687E-05  -0.60333E-03
  //14   r^6_10     -0.13933       0.19729E-02   0.47813E-05  -0.13933    
  //15   r^6_1-1     0.86775E-03   0.25712E-02   0.12674E-05   0.86775E-03
  //16   r^3_10      0.48685E-02   0.12843E-01   0.30190E-04   0.48685E-02
  //17   r^3_1-1     0.35502E-01   0.15227E-01   0.20047E-03   0.35502E-01
  //18   r^7_10      0.16234E-01   0.16809E-01   0.22204E-03   0.16234E-01
  //19   r^7_1-1     0.71746E-02   0.23442E-01   0.30024E-03   0.71746E-02
  //20   r^8_11      0.31996E-02   0.17392E-01   0.99061E-03   0.31996E-02
  //21   r^8_00      0.27686E-01   0.35493E-01   0.42193E-03   0.27686E-01
  //22   r^8_10      0.36608E-01   0.15852E-01   0.21436E-03   0.36608E-01
  //23   r^8_1-1     0.13810E-03   0.21497E-01   0.55996E-04   0.13810E-03
///////////////////////////////////////////////////////////////////////////
// Please  note that in presentation SDME's are ordered in CLASSES
/////////////////////////////////////////////////////////////////////////
double   HPhysicsGenOMEGA::sdme(double cos_t, double phi_d, double phi_m, int ii){
//double cos_t,phi_d,phi_m;
//   ii sign of beam 
  double wl1,wl2,wl3;
  double sd[23];
  double wu=0.0;
  double wl=0.0; 
  double ww=0.0; 
  //const double epsl=0.6112; for W 12-17 GeV  
  const double epsl=0.9545;  
  const double pol_beam=0.8;
  double wu1=0.0;
  double wu2=0.0; 
  double wu3=0.0; 
  double wu4=0.0;
  double sin_2t,pi_1,pi_2;
  // sdme zb sdme == rho-uncamis-Q0P0W0.cutd
  // ~/COMPASS/DATA-XI/UNBINNED/
  // average value of epsilon //
  sd[ 0]= 0.346;  sd[ 1]=  0.018;  sd[ 2]=  0.022;  sd[ 3]= -0.025;
  sd[ 4]=-0.078;  sd[ 5]= -0.081;  sd[ 6]= -0.041;  sd[ 7]=  0.061;
  sd[ 8]= 0.031;  sd[ 9]= -0.016;  sd[10]=  0.132;  sd[11]=  0.103;
  sd[12]=-0.017;  sd[13]= -0.089;  sd[14]=  0.023;  sd[15]=  0.057;
  sd[16]= 0.095;  sd[17]=  0.005;  sd[18]=  0.150;  sd[19]= -0.106;
  sd[20]= 0.125;  sd[21]=  0.093;  sd[22]= -0.009;
  
  sin_2t =2.0*cos_t*sqrt(1.0-cos_t*cos_t);
  
  wu =     0.5*(1.0-sd[0])+0.5*(3.0*sd[0]-1.0)*cos_t*cos_t;
  wu = wu - sqrt(2.0)*sd[1]*sin_2t*cos(phi_m)-sd[2]*(1.0-cos_t*cos_t)*cos(2.0*phi_m) ;
  
  wu1= sd[3]*(1.0-cos_t*cos_t) + sd[4]*cos_t*cos_t -sqrt(2.0)*sd[5]*cos(phi_m)*sin_2t;
  
  wu1=wu1 -sd[6]*(1.0-cos_t*cos_t)*cos(2.0*phi_m);
  wu = wu -epsl*cos(2.0*phi_d)*wu1;
  
  wu2=sqrt(2.0)*sd[7]*sin_2t*sin(phi_m) + sd[8]*sqrt(1.0-cos_t*cos_t)*sin(2.0*phi_m) ;
  wu= wu - epsl*sin(2.0*phi_d)*wu2;
  
  wu3=sd[9]*(1.0-cos_t*cos_t) + sd[10]*cos_t*cos_t -sqrt(2.0)*sd[11]*sin_2t*cos(phi_m);
  wu3=wu3-sd[12]*(1.0-cos_t*cos_t)*cos(2.0*phi_m);
  // wu3=wu3+sd[12]*(1.0-cos_t*cos_t)*cos(2.0*phi_m);
  wu= wu + sqrt(2.0*epsl*(1.0+ epsl))*cos(phi_d)*wu3;
  
  wu4= sqrt(2.0)*sd[13]*sin_2t*sin(phi_m) + sd[14]*(1.0-cos_t*cos_t)*sin(2.0*phi_m);
  wu4= wu4*sqrt(2.0*epsl*(1.0+epsl))*sin(phi_d);
  wu= wu+wu4;
  pi_1= M_PI     ;
  pi_2= pi_1*pi_1;
  wu= wu*(3.0/(8.0*pi_2)); 
  //////////////////////////////////////////////////////// 
  wl1= sqrt(2.0)*sd[15]*sin_2t*sin(phi_m) + sd[16]*(1.0-cos_t*cos_t)*sin(2.0*phi_m);
  wl1= sqrt(1.0-epsl*epsl)*wl1;
  
  wl2= sd[17]*sqrt(2.0)*sin_2t*sin(phi_m)+ sd[18]*(1.0-cos_t*cos_t)*sin(2.0*phi_m);
  wl2= wl2*sqrt(2.0*epsl*(1.0-epsl))*cos(phi_d);
  
  wl3= sd[19]*(1.0-cos_t*cos_t)+sd[20]*cos_t*cos_t - sd[21]*sqrt(2.0)*sin_2t*cos(phi_m);
  wl3=wl3 -sd[22]*(1.0-cos_t*cos_t)*cos(2.0*phi_m);
  wl3=wl3*sqrt(2.0*epsl*(1.0-epsl))*sin(phi_d);
  wl=(wl1+wl3+wl3)*(3.0*(-1)*ii*pol_beam)/(8.0*pi_2);
  // -1*ii*pol_beam helicity of beam
  return wu+wl;
}
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
double  HPhysicsGenOMEGA::cos_theta_vm(HLorentzVector LzVecOut1, HLorentzVector LzVecOut2,  HLorentzVector LzVecOut3, HLorentzVector pr_r,   int i4)
{
  double ss=0. ;
  HVector3 lv1,lv2,b;
  HLorentzVector LzVecOut ;

  LzVecOut= LzVecOut1+LzVecOut2+LzVecOut3;
  lv1= LzVecOut1.Vect(); 
  lv2= LzVecOut2.Vect(); 
  b = - LzVecOut.BoostVector();                  // boost vector of            p1+p2+p3 system
  HLorentzVector lv1p(LzVecOut1); lv1p.Boost(b); // transform p1    4vector to p1+p2+p3 system
  HLorentzVector lv2p(LzVecOut2); lv2p.Boost(b); // transform p2    4vector to p1+p2+p3 system
  HLorentzVector lvr(pr_r);        lvr.Boost(b); // transform proton vector to p1+p2+p3 system
  HVector3 llv1 =  lv1p.Vect();
  HVector3 llv2 =  lv2p.Vect();
  HVector3 enne =  llv1.Cross(llv2).Unit();
  HVector3 llpr =  -lvr.Vect(); // change sign
  ss = cos(llpr.Angle(enne));
  return ss ; 
} 
/////////////////////////////////////////////////////////////////////////////
