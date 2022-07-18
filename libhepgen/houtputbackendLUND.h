/*!
 *  \file houtputbackendLUND.h
 *
 *  \date Created on: July 18, 2022
 *  \author Andrey Kim<kenjo@jlab.org>
 *  Copyright (c) 2022 All Right Reserved
 */


#ifndef HOUTPUTBACKENDLUND_H_
#define HOUTPUTBACKENDLUND_H_


#include "houtputbackend.h"
#include "hparammanager.h"
#include "hevent.h"
#include <iostream>
#include <fstream>


/*! \brief Output-Backend for events in CLAS12 LUND-format */
class HOutPutBackEndLUND: public HOutPutBackEnd
{
public:
    void dumpEvent(HEvent* event);
    void closeFile();
    void initFile(std::string _fileName, HEvent* _eventPointer);


private:
    std::string fileName;
    ofstream outFile;
};


#endif

