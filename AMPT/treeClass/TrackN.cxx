#include "TrackN.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

ClassImp(TrackN);

//
// Track class
//

//_________________________________________________________
TrackN::TrackN():
TObject(),
pdgCode(-999999999),
status(-999),
px(-999),
py(-999),
pz(-999),
e(-999),
phi(-999),
eta(-999),
q(-999)
{
  // default constructor
}


//_________________________________________________________
TrackN::~TrackN(){}
