/*
  This is a class that is build to unique-ify the values of the parameters 
*/
#include <iostream>
#include <map>
#include <utility>
#include <algorithm>
#include <vector>

#include <string.h>
#include <math.h>

using namespace std;

class BondedParams {
private:
  int natoms;  
  int atom[4]; // indices of the atoms involved in this bonded term

  static const int maxparms = 8;
  int nparms; // number parameters used to describe this bonded term in the force field
  double parm[maxparms]; // parameters used to describe this bonded term in the force field

public:
  BondedParams() { natoms = nparms = 0;}
  BondedParams(int n) { natoms = 0; nparms = n;}
  // copy constructor
  BondedParams(const BondedParams& B) { 
    natoms = B.natoms;
    memcpy((void *)atom, (void *)B.atom, 4*sizeof(int));
    nparms = B.nparms;
    memcpy((void *)parm, (void *)B.parm, maxparms*sizeof(double));
  }

  void Set_Atoms(int a, int b) { 
    natoms=2; 
    if(a<b) {atom[0]=a; atom[1]=b;}
    else    {atom[0]=b; atom[1]=a;}
  }
  void Set_Atoms(int a, int b, int c) { 
    natoms=3; 
    if(a < c) { atom[0]=a; atom[1]=b; atom[2]=c;}
    else {      atom[0]=c; atom[1]=b; atom[2]=a;}
  }
  void Set_Atoms(int a, int b, int c, int d) { 
    natoms=4; 
    if(a<d){ atom[0]=a; atom[1]=b; atom[2]=c; atom[3]=d;}
    else {   atom[0]=d; atom[1]=c; atom[2]=b; atom[3]=a;}    
  }
  void Set_Atoms_AsIs(int a, int b, int c, int d) { 
    natoms=4; 
    atom[0]=a; atom[1]=b; atom[2]=c; atom[3]=d;
  }

  void Set_Atoms_Inversion(int a, int b, int c, int d)
  {
    vector<int> v;
    natoms=4;
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    sort(v.begin(),v.end());
    atom[0] = a;
    atom[1] = v[0];
    atom[2] = v[1];
    atom[3] = v[2];
  }

  void Set_Parm(int id, double val) { if(id < nparms) parm[id] = val; }

  void Set_NParms(int n) { nparms = n; }
  void Set_NAtoms(int n) { natoms = n; }

  int Get_Atom(int id) const { return (id < natoms ? atom[id] : -1); }
  double Get_Parm(int id) const { if(id < nparms) return parm[id]; else return -1.0;}

  void print() const { 
    int i;
    cout << "*** ";
    for(i=0;i<natoms;i++) cout << atom[i] << " ";
    cout << "params: ";
    for(i=0;i<nparms;i++) cout << parm[i] << " ";
    cout << " ***" << endl;
  }
  
  bool operator< (const BondedParams &b2) const
  {
    int i;

    if(this->natoms < b2.natoms) return true;
    if(this->nparms < b2.nparms) return true;

    // get the first different value
    i=0;
    while( (fabs(this->Get_Parm(i) - b2.Get_Parm(i)) < 1.0e-5) && i < this->nparms) {
      i++;
    }

    // if it is the last, then they're equal
    if(i == this->nparms) {
      return false;
    } else {
      // return the comparison
      return(this->Get_Parm(i) < b2.Get_Parm(i));
    }

  }

  bool operator== (const BondedParams &b2) const
  {
    int i;
    double diff;

    if(this->natoms != b2.natoms) return false;
    if(this->nparms != b2.nparms) return false;
            
    for(i=0; i < this->nparms; i++) {
      diff = fabs( this->Get_Parm(i)-b2.Get_Parm(i) );
      if(diff > 1.0e-4 ) 
	{
	  return false;
	}
    }
    
    return true;
  }
    
};

void print_multimap(multimap<BondedParams, int> &m)
{
  multimap<BondedParams, int>::iterator it;
  vector<BondedParams> ub;
  BondedParams lastkey;

  for( it=m.begin(); it != m.end(); it++) {
    (*it).first.print();
    cout << (*it).second << endl;
  }
}

// Returns a vector with the unique parameters found
vector<BondedParams> unique_params(multimap<BondedParams, int> &m)
{
  multimap<BondedParams, int>::iterator it;
  vector<BondedParams> ub;
  BondedParams lastkey;

  for( it=m.begin(); it != m.end(); it++) {
    if(!(lastkey == (*it).first)) {
      lastkey = (*it).first;
      ub.push_back(lastkey);
    }
  }

  return ub;
}

int bond_type(BondedParams &bp, vector<BondedParams> &ub)
{
  vector<BondedParams>::iterator it;
  int ret = 0;

  for( it=ub.begin(); it != ub.end(); it++, ret++) {
    if( *it == bp ) return(ret);
  }

  return -1; // error
}
