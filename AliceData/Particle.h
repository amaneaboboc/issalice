class Particle
{


protected:

  Particle * type;      //!< type of this particle
  //vector<Particle*> parents;  //!< array containing the parents of this particle.
  //vector<Particle*> children; //!< array containing the children of this particle.
  Particle * truth;  //!< pointer to the truth particle corresponding to this particle.
  //bool live; //!< whether this particle is live or dead (measurable or not)
  long pid;  //!< used defined identified used in some applications
  int  sourceIndex;  //!<  source index  used in some applications


public:

Particle()
:
momentum (),
position (),
type     (nullptr),
//parents  (),
//children (),
//truth    (nullptr),
//live     (false),
pid      (-1)
//sourceIndex(-1),
//ixEtaPhi (0),
//ixYPhi   (0)
{

}
/*
virtual ~Particle(){}

  TString getName()
  {
  if (type)
    return type->getName();
  else
    return "UnknownType";
  }

void set(Particle* _type,
                   double p_x, double p_y, double p_z, double p_e,
                   double _x,  double _y,  double _z,  double _t,
                   )
{
  //clear();
  type = _type;
  momentum.SetPxPyPzE (p_x,p_y,p_z,p_e);
  position.SetXYZT    (_x,_y,_z,_t);
  //live       = _live;
  //parents.clear();
  //children.clear();
  //truth = nullptr;
}
*/
};
