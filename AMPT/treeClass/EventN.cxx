#include "EventN.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

//_____________________________________________________________________________
ClassImp(EventN)

EventN::EventN():
TObject(),
multChEta1(-1),
nMPI(-1)
{
    // default constructor
}

//_________________________________________________________
EventN::~EventN(){}
