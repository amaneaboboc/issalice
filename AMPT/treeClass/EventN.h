#ifndef EVENTN_H
#define EVENTN_H

#include <TObject.h>

class EventN : public TObject
{
public:
    Int_t     multChEta1;    // multiplicity |eta|<1
    Int_t     nMPI;         // number of MPIs
    
    EventN(); // default constructor
    virtual ~EventN(); // default destructor
  
    ClassDef(EventN, 1)  // Event information
    
};

#endif
