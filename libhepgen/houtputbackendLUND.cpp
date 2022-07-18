#include "houtputbackendLUND.h"


void HOutPutBackEndLUND::closeFile()
{
    outFile.flush();
    outFile.close();
}


// outputs event in CLAS12 LUND format: https://gemc.jlab.org/gemc/html/documentation/generator/lund.html

void HOutPutBackEndLUND::dumpEvent(HEvent* event)
{
    int nactive = 0;
    for (unsigned int ipart=0; ipart < eventPointer->getStruct()->listOfParticles.size(); ipart++)
      if(event->getStruct()->listOfParticles.at(ipart)->getParticleAuxFlag() == 1)
        nactive++;

    outFile << nactive << " " << 0 << " ";
    outFile << event->getStruct()->PARL.at(1) << " " << event->getStruct()->PARL.at(0) << " ";
    outFile << 0 << " " << event->getStruct()->epp_polarized_longitudinal << " ";
    outFile << event->getStruct()->incBeamParticle.getParticleType() << " ";
    outFile << event->getStruct()->PARL.at(20) << " ";
    outFile << event->getStruct()->targetParticle.getParticleType() << " ";
    outFile << 0 << " " << 0 << std::endl;

    int iactive = 0;
    for (unsigned int ipart=0; ipart < eventPointer->getStruct()->listOfParticles.size(); ipart++) {
      if(event->getStruct()->listOfParticles.at(ipart)->getParticleAuxFlag() == 1) {
        iactive++;

        outFile << "  " << iactive << " " << 0 << " " << 1 << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getParticleType() << " " << 0 << " " << 0 << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getVector().getVector().X() << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getVector().getVector().Y() << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getVector().getVector().Z() << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getEnergy() << " ";
        outFile << event->getStruct()->listOfParticles.at(ipart)->getMass() << " ";
        outFile << event->getStruct()->USERVAR.at(0) << " " << event->getStruct()->USERVAR.at(1) << " ";
        outFile << 0 << std::endl;

      }
    }

    outFile.flush();
}




void HOutPutBackEndLUND::initFile(string _fileName, HEvent* _eventPointer)
{
    fileName = _fileName+".lund";
    eventPointer = _eventPointer;
    outFile.open(fileName.c_str());
}


