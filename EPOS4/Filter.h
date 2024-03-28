//#ifndef

class Condition
{
public:
  int filterType;
  int filterSubtype;
  double minimum;
  double maximum;

  Condition (int _filterType, int _filterSubtype, double _minimum,
             double _maximum)
      : filterType (_filterType), filterSubtype (_filterSubtype),
        minimum (_minimum), maximum (_maximum)
  {
  }

  virtual ~Condition () {}

  virtual bool
  accept (double value)
  {
    switch (filterType)
      {
      default:
        return (value >= minimum && value < maximum);
      case 2:
      case 3:
        return std::fabs (value - filterSubtype) < 0.1;
      }
  }

  virtual bool
  accept (int value)
  {
    switch (filterType)
      {
      default:
        return (value >= minimum && value < maximum);
      case 2:
      case 3:
        return filterSubtype == value;
      }
  }
};

template <class T> class Filter
{
protected:
  TString name;
  TString title;
  std::vector<Condition *> conditions;

public:
  Filter () : name (""), title (""), conditions ()
  {
    //---
  }

  virtual ~Filter ()
  {
    for (unsigned int k = 0; k < conditions.size (); k++)
      {
        delete conditions[k];
      }
    conditions.clear ();
  }

  virtual bool
  accept (const T &object __attribute__ ((unused)))
  {
    return true;
  }

  void
  addCondition (unsigned int type, unsigned int subtype, double minimum,
                double maximum)
  {
    Condition *condition = new Condition (type, subtype, minimum, maximum);

    // cout<"Filter:"<<name<<endl;

    cout << "type:" << type << endl;

    cout << "subtype:" << subtype << endl;

    cout << "min:" << minimum << endl;

    cout << "max:" << maximum << endl;

    conditions.push_back (condition);
  }

  inline unsigned int
  getNConditions () const
  {
    // cout<<"Nr. Conditii: "<< conditions.size()<<endl;
    return conditions.size ();
  }

  void
  setTitle (const TString &newTitle)
  {
    title = newTitle;
  }

  TString
  getTitle ()
  {
    return title;
  }

  void
  setName (const TString &newName)
  {
    name = newName;
  }

  TString
  getName ()
  {
    return name;
  }
};
