//#ifndef"Filter.h"

#include "Filter.h"
//#include "Particle.h"

#include "TrackN.h"

class ParticleFilter: public Filter<TrackN>
{
public:

ParticleFilter()
:
Filter<TrackN>()
{

}

bool accept(const TrackN * particle)
{

unsigned int nConditions = getNConditions();
//cout<<nConditions<<endl;
if (nConditions<1) return true;
double value;
bool accepting;

for (unsigned int k = 0; k<nConditions; k++)
{
Condition & condition = *(conditions[k]);
unsigned int filterType = condition.filterType;
unsigned int filterSubType = condition.filterSubtype;

//cout << "filterType: " << filterType <<endl;
//cout<< "filterSubtype: "<< filterSubType<<endl;

  switch (filterType)
  {
  
  	case 0:
  	{
  	switch (filterSubType)
  	  {
  	  case 0: accepting = 0; break; //decayed particles
  	  case 1: accepting = 1; break; //undecayed particles
  	  }
  	  //cout<< " case 0: accepting "<<accepting;
  	}
  	
  	break;
  	
  	case 1://charge separation
  	{
  	  int charge = 1;
  	  switch (filterSubType)
  	  {
  	  case 0: accepting = (particle->q==0); break;
  	  case 1: accepting = (particle->q!=0); break;
  	  case 2: accepting = (particle->q > 0.);  break;
  	  case 3: accepting = (particle->q < 0.);  break;
  	  }
  	  //cout<< " case 1: accepting "<<accepting;
  	}
	break;
	
	case 2: //PDG code
	{
	//accepting = condition.accept(particle->pdgCode);
	
	
	}
	
	//cout<< " case 2: accepting "<<accepting;
	
	break;

	}
	if(!accepting) return false;
    }
return true;
}




};



//#endif
