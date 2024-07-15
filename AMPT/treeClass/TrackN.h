#ifndef TRACKN_H
#define TRACKN_H

#include <TObject.h>

class TrackN : public TObject
{
 public:
    Int_t     pdgCode; //pdg code
    Int_t     status;  //status
    Float_t   px;      //px
    Float_t   py;      //py
    Float_t   pz;      //pz
    Float_t   e;       //energy
    Float_t   phi;     //phi
    Float_t   eta;     //pseudorapidity
    Float_t   q;       //charge

    TrackN();
    virtual ~TrackN(); // default destructor
  
    ClassDef(TrackN, 1);    // Help class

};

#endif
