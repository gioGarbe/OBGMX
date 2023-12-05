/**********************************************************************
obgmx.cpp - generates GROMACS topologies using the OpenBabel Library

Copyright (C) 2012 Giovanni Garberoglio - garberoglio _at_ fbk _dot_ eu 
Interdisciplinary Laboratory for Computational Science (LISC)
FBK-CMM and University of Trento (more the former than the latter)
http://lisc.fbk.eu

Please cite: 
G. Garberoglio. OBGMX: A web-based generator of GROMACS topologies for
molecular and periodic systems using the Universal Force Field.
J. Comput. Chem. 33 (2012), 2204-2208. doi: 10.1002/jcc.23049

Based on obenergy.cpp - Copyright (C) 2006 Tim Vandermeersch

this program needs a -gg openbabel distribution, which knows about periodic
molecules and exports the forcefield data.
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/plugin.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/mol.h>
#include <openbabel/chargemodel.h>
#include <openbabel/obutil.h>


#ifndef _MSC_VER
  #include <unistd.h>
#endif

#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bondedparams.h"

#define USE_BOND_LENGTHS (1<<0)
#define USE_ANGLES       (1<<1)
#define USE_DIHEDRALS    (1<<2)
#define USE_INVERSIONS   (1<<3)

using namespace std;
using namespace OpenBabel;

void write_credits(ofstream *s)
{
  char line[512];

  sprintf(line,"; Topology generated using OBGMX.\n;\n");
  s->write(line,strlen(line));
  sprintf(line,"; Please cite:\n");
  s->write(line,strlen(line));
  sprintf(line,"; G. Garberoglio. OBGMX: A web-based generator of GROMACS topologies for\n");
  s->write(line,strlen(line));
  sprintf(line,"; molecular and periodic systems using the Universal Force Field.\n");
  s->write(line,strlen(line));
  sprintf(line,"; J. Comput. Chem. 33 (2012), 2204-2208. doi: 10.1002/jcc.23049\n\n");
  s->write(line,strlen(line));
}

int main(int argc,char **argv)
{
  int i, j;
  char *program_name= argv[0];
  int c;
  unsigned int natoms = 0;
  bool hydrogens = false;
  string basename, filename = "", option, option2, ff = "UFF";
  char logfilename[256] = "obgmx.ffout";
  OBMol mol, tmpmol;

  char line[512], foo[512];

  char residuename[16] = "UFF";
  char molname[256] = "obgmx";
  char topfilename[256] = "obgmx.top";
  char itpfilename[256] = "obgmx.itp";

  // combination rules for non-bonded interaction
#define LORENTZ_BERTHELOT 1
#define GEOMETRIC_RULE    2
  int comb_rule = GEOMETRIC_RULE;

  // 1-4 multiplication factor.
  // no indication in the UFF paper, apart from the neglect of the 1-2
  // and 1-3 interactions
  double mult14 = 1.0;            

  // functional form for the angular potential
#define GROMACS_HARMONIC_ANGLE 1
#define GROMACS_COSINE_ANGLE   2
  int angletype = GROMACS_COSINE_ANGLE;

  OBConversion conv;
  OBFormat *format_in = NULL;

  // in linux, I have to make these declarations, otherwise nothing will work
  //OBChargeModel foo1("gasteiger",true);
  //OBChargeModel foo2("mmff94",false);
  //OBChargeModel foo3("qeq",false);
  //OBChargeModel foo4("qtpie",false);

  //OBChargeModel *chg_model = NULL;
  char *chg_method;
  std::vector<double> partialCharges;
  std::vector<double> formalCharges;
  // atomic partial charges. not recommended in UFF.
  double *charge = NULL;
  double total_charge = 0.0, formal_charge = 0.0;
  double dipole[3];
  int use_geometry = 0;

  bool debug = false;

  chg_method = (char *)calloc(256,sizeof(char));

  if (argc < 2) {
    cout << "Usage: " << argv[0] << " [options] <filename>" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -h          add hydrogens before calculating energy" << endl;
    cout << "  -d          include debug statements" << endl;
    cout << "  -LB         use Lorentz-Berthelot mixing rule [def. geometric]" << endl;
    cout << "  -14 mult    1-4 multiplication parameter [def. 1.0]" << endl;
    cout << "  -H          use harmonic form for angle potential [def. cosine]" << endl;
    cout << "  -G idx      use geometrical data. An OR of:" << endl;
    cout << "              " << USE_BOND_LENGTHS << ": bond lengths" << endl;
    cout << "              " << USE_ANGLES       << ": angles" << endl;
    cout << "              " << USE_DIHEDRALS    << ": dihedral" << endl;
    //cout << "              " << USE_INVERSIONS   << ": inversions" << endl;

    cout << endl;
#ifndef OBGMX_RELEASE
    cout << "  -ff ffid    select a forcefield" << endl;
    cout << "              available forcefields:" << endl;
    OBPlugin::List("forcefields", "verbose");
    cout << endl;
    cout << "  -cm model    select a charge model" << endl;
    cout << "              available charge models:" << endl;
    OBPlugin::List("charges", NULL);
    cout << endl;
#endif
    exit(-1);
  } else {
    int ifile = 1;
    for (int i = 1; i < argc; i++) {
      option = argv[i];
      
      if (option == "-h") {
        hydrogens = true;
        ifile++;
      }

      if (option == "-LB") {
	comb_rule = LORENTZ_BERTHELOT;
        ifile++;
      }

      if (option == "-d") {
	debug = true;
        ifile++;
      }

      if (option == "-H") {
	angletype = GROMACS_HARMONIC_ANGLE;
        ifile++;
      }
     
      if ((option == "-14") && (argc > (i+1))) {
        mult14 = atof(argv[i+1]);
        ifile += 2;
      }


      if ((option == "-G") && (argc > (i+1))) {
        use_geometry = atoi(argv[i+1]);
        ifile += 2;
      }
      
      if ((option == "-ff") && (argc > (i+1))) {
        ff = argv[i+1];
        ifile += 2;
      }

      if ((option == "-cm") && (argc > (i+1))) {
        strcpy(chg_method, argv[i+1]);
        ifile += 2;
      }

    }
    
    basename = filename = argv[ifile];
    size_t extPos = filename.rfind('.');

    if (extPos!= string::npos) {
      basename = filename.substr(0, extPos);
    }

  }

  cout << "\n" << endl;
  cout << "############################################################" << endl;
  cout << "# OBGMX: a generator of GROMACS topologies using OpenBabel #" << endl;
  cout << "############################################################" << endl;
  cout << "\n" << endl;


  // Find Input filetype
  format_in = conv.FormatFromExt(filename.c_str());    
  if (!format_in || !conv.SetInFormat(format_in)) {
    cerr << program_name << ": cannot read input format!" << endl;
    exit (-1);
  }

  // Fake-output filetype is MOL2, so we have the charges
  conv.SetOutFormat("MOL2");

  ifstream ifs;
  ofstream ofs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBForceField* pFF = OBForceField::FindForceField(ff);
  if (!pFF) {
    cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
    exit (-1);
  }


  ofstream *outstream = new ofstream;
  outstream->open(logfilename);
  cout << "Setting output to " << logfilename << endl;

  //pFF->SetLogFile(&std::cout);
#ifdef CRIPPLE
  pFF->SetLogLevel(OBFF_LOGLVL_LOW);
#else
  pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
#endif
  pFF->SetLogFile(outstream);

  // end of command line parsing
  if(strlen(chg_method)) {
    conv.AddOption("partialcharges",OBConversion::INOPTIONS,chg_method);
  }


  double energy;
  for (c=1;;c++) {
    tmpmol.Clear();

    if (!conv.Read(&tmpmol, &ifs)) {
      conv.Convert();
      break;
    }
    if (tmpmol.Empty()) break;
    mol = tmpmol;

    if (hydrogens)
      mol.AddHydrogens();

    natoms = mol.NumAtoms();

#ifdef CRIPPLE
#define MAX_ATOMS 1000
    if(natoms > MAX_ATOMS) 
      {
	// create obgmx.top and obgmx.itp with warning
	ofstream *s = new ofstream;
	s->open("obgmx.top");
	write_credits(s);
	sprintf(line,"Cannot create topology.\n");
	s->write(line,strlen(line));
	sprintf(line,"This version is limited to %d atoms and you requested %d.\n",MAX_ATOMS,natoms);
	s->write(line,strlen(line));
	s->close();
	s->open("obgmx.itp");
	write_credits(s);
	sprintf(line,"Cannot create topology.\n");
	s->write(line,strlen(line));
	sprintf(line,"This version is limited to %d atoms and you requested %d.\n",MAX_ATOMS,natoms);
	s->write(line,strlen(line));
	s->close();
	exit(10);
      }
#endif
    
    if (!pFF->Setup(mol)) {
      cerr << program_name << ": could not setup force field." << endl;
      exit (-1);
    }

    /* ok. here we have the molecule. let's rock! */
    cout << "Molecule name: '" << mol.GetTitle() << "'" << endl;

    // default to zero partial charges
    charge = new double[natoms]; 
    for(i=0; i<natoms;i++) charge[i] = 0.0;

    // see if we have to get charges
#ifndef OBGMX_RELEASE
    if(strlen(chg_method)) 
      {
        
        cout << "Trying to use the charge method: " << chg_method << endl;
        
        OBChargeModel *chg_model = OBChargeModel::FindType(chg_method);
        
        if(chg_model) {
          
          mol.UnsetPartialChargesPerceived();
          total_charge = 0.0;
          dipole[0] = dipole[1] = dipole[2] = 0.0;
          
          
          if( chg_model->ComputeCharges(mol) ) 
            {
              partialCharges = chg_model->GetPartialCharges();
              formalCharges  = chg_model->GetFormalCharges();
              
            } 
          else 
            {
              cout << "Cannot compute charges using: " << chg_method << endl;
            }
        }
        else 
          {
            cout << "Cannot setup partial charge calculation using " << chg_method << endl;
            cout << "Reverting to gasteiger which is the default" << endl;
            sprintf(chg_method,"%s","gasteiger");
            
            cout <<
              "----------------------------------------------------------------------"
                 << endl;
          }
      
        OBAtom *atom;
        i = 0; 
        FOR_ATOMS_OF_MOL(atom,mol) 
          {
            //cout << atom->GetType() << " " << atom->GetPartialCharge() << endl;
            charge[i] = atom->GetPartialCharge();
            total_charge += charge[i];
            dipole[0] += charge[i] * atom->x();
            dipole[1] += charge[i] * atom->y();
            dipole[2] += charge[i] * atom->z();
            i++;
          }      
        cout << "Total  charge: " << total_charge << endl;
        cout << "Dipole moment: " << dipole[0] << " " << dipole[1] << " " << dipole[2] << endl;
        cout << "----------------------------------------------------------------------" << endl;           
      }      
#endif
    // Calculate the energy and - most important - fill the logfile
    // with all the information which is needed :)
    energy = pFF->Energy(false);
    outstream->close();
    
  } // end for loop


  char names[natoms][5];
  ifstream infile;
  unsigned int nbonds, nangles, ntors, noop;

  // Here we should in principle do different things for different
  // force fields. Right now we are limited to UFF.

  if(ff != "UFF")
    {
      cout << "Topology generation limited to UFF, but you chose " << ff << endl;
      cout << "exiting..." << endl;
      exit(2);
    }

  cout << "Setting up UFF parameters for " << natoms << " atoms" << endl;

  // MAKING THE TOP FILE WITH THE DEFINITION OF THE FORCE FIELD
  // AND ATOMIC TYPES

  ofstream *topstream = new ofstream;

  sprintf(topfilename,"%s.top",molname);
  topstream->open(topfilename);
  cout << "Saving topology to: " << topfilename << endl;

  write_credits(topstream);
  sprintf(line,"[ defaults ]\n");
  topstream->write(line,strlen(line));
  sprintf(line,"; nbfunc      comb-rule      gen-pairs       fudgeLJ    fudgeQQ\n");
  topstream->write(line,strlen(line));
  sprintf(line,"1             %d              yes             %.1f        %.1f\n\n",
	  comb_rule,mult14,mult14);
  topstream->write(line,strlen(line));

  sprintf(line,"[ atomtypes ]\n");
  topstream->write(line,strlen(line));
  sprintf(line,"; name1 name2   mass     charge  ptype   sigma   epsilon\n");
  topstream->write(line,strlen(line));

#ifndef OBGMX_RELEASE
  if(strlen(chg_method)) 
    {
      sprintf(line,"; charges set using the %s model in the itp file\n",chg_method);
      topstream->write(line,strlen(line));
    }
#endif
  
  for(i=0; i<pFF->NAtomTypes(); i++) {
    const double cfact = 1000.0/(6.022*1.38);
    sprintf(line,"%-6s   %-6s   %9.4f   %9.4f   A   %9.4f    %9.4f\n",
	    pFF->TypeName(i),pFF->TypeName(i),
	    pFF->MassOfAtomType(i), 
	    0.0,
	    0.1*pFF->LJSig(pFF->TypeName(i)),
	    pFF->LJEps(pFF->TypeName(i)));
    topstream->write(line,strlen(line));
  }
    
  sprintf(line,"\n#include <%s.itp>\n\n",molname);
  topstream->write(line,strlen(line));

  sprintf(line,"[ system ]\nUFF MOLECULE\n\n");
  topstream->write(line,strlen(line));

  sprintf(line,"[ molecules ]\n%-6s   1\n",residuename);
  topstream->write(line,strlen(line));

  topstream->close();

  // NOW LET'S MAKE THE ITP FILE RELATIVE TO THIS SYSTEM

  ofstream *itpstream = new ofstream;

  vector<string> outlist;
  string line_as_string;
  vector<string>::iterator it;

  // the following variables are used on each of the four possible
  // bonded potentials (bond stretching, angle bending, proper and
  // improper torsions)

  // see the BONDS section below for an explanation

  // a multimap is used to find the number of equivalent parameters
  multimap<BondedParams, int> m;
  BondedParams B;
  vector<BondedParams> Bvec; // a vector with all the bonded parameters found
  vector<BondedParams> ub;   // a vector with the UNIQUE set of parameters
  vector<BondedParams>::iterator bit;
    
  sprintf(itpfilename,"%s.itp",molname);
  itpstream->open(itpfilename);
  cout << "Saving molecular topology to: " << itpfilename << endl;

  write_credits(itpstream);
  sprintf(line,"[ moleculetype ]\n; Name       nrexcl\n");
  itpstream->write(line,strlen(line));
  sprintf(line,"%s    3\n",residuename);
  itpstream->write(line,strlen(line));

  sprintf(line,"\n[ atoms ]\n");
  itpstream->write(line,strlen(line));
  sprintf(line,"; nr type  resnr    residue    atom     cgnr    charge       mass \n");
  itpstream->write(line,strlen(line));

#ifndef OBGMX_RELEASE
  if(strlen(chg_method)) 
    {
      sprintf(line,"; charges set using the %s model\n",chg_method);
      itpstream->write(line,strlen(line));
      sprintf(line,"; total charge: %f\n",total_charge);
      itpstream->write(line,strlen(line));
      sprintf(line,"; dipole moment: %f %f %f\n",dipole[0],dipole[1],dipole[2]);
      itpstream->write(line,strlen(line));
    }
#endif
  if(debug)
    cout << "#\n# CONNECTIVITY (indices start from 1)" << endl;
    
  for(i=0; i<natoms; i++) {
    int n = pFF->AtomTypeAsNumber(i);
    sprintf(line,"%-6d     %-6s   1   %-6s  %-6s   %d  %9.4f  %9.4f\n",
	    i+1,pFF->TypeName(n),residuename,pFF->TypeName(n),
	    i+1,charge[i],pFF->MassOfAtomType(n));
    //cout << line;
    itpstream->write(line,strlen(line));

    // debug statement: print the neighbors for each atom
    if(debug)
      {
	unsigned int nneigh = 0;
	OBAtom *atom, *neigh;
	atom = mol.GetAtom(i+1);
	cout << "# Neighbors of " << pFF->TypeName(n) << "(" << i << "): ";
	
	FOR_NBORS_OF_ATOM(neigh,atom)
	  {
	    unsigned int idx = neigh->GetIdx();
	    cout << pFF->TypeName(pFF->AtomTypeAsNumber(neigh->GetIdx()-1));
	    cout << "(" << idx << ") ";
	    nneigh++;
	  }
	cout << "--- " << nneigh << endl;
      }

  }
  if(debug) cout << endl << "#" << endl ;
	
  // BONDS
  sprintf(line,"\n[ bonds ]\n");
  itpstream->write(line,strlen(line));

  B.Set_NAtoms(2); // a bond-stretching terms involves two atoms
  B.Set_NParms(2); // and is described by two parameters in UFF

  for(i=0; i<pFF->NBonds(); i++) {
    double r0 = 0.1   * pFF->BondParam(i,2);
    double kb = 200.0 * pFF->BondParam(i,1);

    int idx1 = pFF->BondAtomIdx(i,0);
    int idx2 = pFF->BondAtomIdx(i,1);

    if(use_geometry & USE_BOND_LENGTHS) {
      OBAtom *a1, *a2;
      a1 = mol.GetAtom(idx1);
      a2 = mol.GetAtom(idx2);
      r0 = 0.1 * a1->GetDistance(a2); // distance in nm
    }

    sprintf(line,"%-6d %-6d  %d  %9.4f  %9.4f ; %-6s %-6s ",
	    idx1, idx2,
	    1,  // harmonic bond in GROMACS
	    r0, // r0 in nm
	    kb, // kb in kj/mol/nm-2
	    pFF->AtomType(pFF->BondAtomIdx(i,0)-1),
	    pFF->AtomType(pFF->BondAtomIdx(i,1)-1));

    line_as_string.assign(line);
    outlist.push_back(line_as_string);
    //itpstream->write(line,strlen(line));
      
    // setup the multimap
    B.Set_Atoms(pFF->AtomTypeAsNumber(pFF->BondAtomIdx(i,0)-1),
		pFF->AtomTypeAsNumber(pFF->BondAtomIdx(i,1)-1));
    //B.Set_Atoms(pFF->BondAtomIdx(i,0),
    //pFF->BondAtomIdx(i,1));
    B.Set_Parm(0,r0);
    B.Set_Parm(1,kb);

    // save the parameters for this particular bond into a a vector
    Bvec.push_back(B);
    m.insert(pair<BondedParams, int>(B,i));
  }

  if(use_geometry & USE_BOND_LENGTHS)
    cout << "Using bond lengths from geometry" << endl;

  ub = unique_params(m); // find the unique terms present in the system
  cout << "FOUND " << ub.size() << " UNIQUE BOND TERMS ";
  cout << "OUT OF " << pFF->NBonds() << endl;

  sprintf(line,"; found %d unique bond terms\n",(int)ub.size());    
  itpstream->write(line,strlen(line));

  // let's build the unique bonds and print them out
  for(i=0,bit=ub.begin(); bit != ub.end(); ++bit,i++) {
    sprintf(line,"; [%d] \t %-6s %-6s     %d     %9.4f  %9.4f\n",i,
	    pFF->TypeName((*bit).Get_Atom(0)),
	    pFF->TypeName((*bit).Get_Atom(1)),
	    1, // harmonic bond in GROMACS
	    (*bit).Get_Parm(0),  // r0 in nm
	    (*bit).Get_Parm(1)); // kb in kj/mol/nm-2      
    cout << line;
    itpstream->write(line,strlen(line));
  }


  for(i=0, it = outlist.begin(); it < outlist.end(); it++, i++) {
    itpstream->write((*it).c_str(), (*it).length() );
    sprintf(line,"type %d\n", bond_type(Bvec[i],ub));
    itpstream->write(line,strlen(line));
  }
    
  // we clean up the variables to be used for the next kind of bonded potential.
  Bvec.clear();
  ub.clear();
  m.clear();
  outlist.clear();

  // ANGLES
  // GROMACS does not have the UFF functional form
  // so we just have to approximate it.
  // we use either harmonic angles or cosine based (GROMOS96) angles
  // with a suitable approximation for the parameters

  B.Set_NAtoms(3);
  B.Set_NParms(2); // theta and kappa

  sprintf(line,"\n[ angles ]\n");
  itpstream->write(line,strlen(line));
  if(angletype == GROMACS_HARMONIC_ANGLE)
    sprintf(line,"; harmonic angles\n");
  if(angletype == GROMACS_COSINE_ANGLE)
    sprintf(line,"; G96 (cosine) angles\n");
  itpstream->write(line,strlen(line));

  for(i=0; i<pFF->NAngles(); i++) {
    int idx1, idx2, idx3;
    int coord = (int)pFF->AngleParam(i,0);
    double ka = pFF->AngleParam(i,1);
    double c0 = pFF->AngleParam(i,2);
    double c1 = pFF->AngleParam(i,3);
    double c2 = pFF->AngleParam(i,4);

    double thetamin, kappa;

    // harmonic approximation of the UFF functional form

    // don't check coord, it's already done in the UFF force field setup
    switch(coord) {
    case 1:
      thetamin = 180.0;
      kappa = ka;
      break;
    case 2:
      thetamin = 120.0;
      kappa = 4.0*ka/3.0;
      break;
    case 4:
    case 6:
      thetamin = 90.0;
      kappa = 2.0*ka;
      break;
    case 7:
      double alpha, c7;
      alpha = 2.0 * M_PI / 5.0;
      c7 = sin(alpha)*(cos(alpha) - cos(2*alpha));
      thetamin = 72.0;
      kappa = 2.0 * c7*c7 * ka * c1;
      break;
    default:
      thetamin = M_PI-acos(c1/(4.0*c2));
      thetamin = thetamin * 180.0 / M_PI;
      kappa = ka * (16.0*c2*c2 - c1*c1) / (4.0* c2);
      break;
    }	

    if(angletype == GROMACS_COSINE_ANGLE && coord != 1) {
      double s = sin(thetamin*M_PI/180.0);
      kappa = kappa / (s*s);
    }

    idx1 = pFF->AngleAtomIdx(i,0);
    idx2 = pFF->AngleAtomIdx(i,1);
    idx3 = pFF->AngleAtomIdx(i,2);

    if(use_geometry & USE_ANGLES) {
      OBAtom *a1, *a2, *a3;
      a1 = mol.GetAtom(idx1);
      a2 = mol.GetAtom(idx2);
      a3 = mol.GetAtom(idx3);
      thetamin = a1->GetAngle(a2,a3); // angle in degrees
    }
      
    sprintf(line,"%-6d %-6d %-6d  %d  %9.4f  %9.4f ; %-6s %-6s %-6s",
	    idx1, idx2, idx3,
	    (coord == 1 ? GROMACS_HARMONIC_ANGLE : angletype), 
	    thetamin,
	    kappa,
	    pFF->AtomType(pFF->AngleAtomIdx(i,0)-1),
	    pFF->AtomType(pFF->AngleAtomIdx(i,1)-1),
	    pFF->AtomType(pFF->AngleAtomIdx(i,2)-1));
      
    B.Set_Atoms(pFF->AtomTypeAsNumber(pFF->AngleAtomIdx(i,0)-1),
		pFF->AtomTypeAsNumber(pFF->AngleAtomIdx(i,1)-1),
		pFF->AtomTypeAsNumber(pFF->AngleAtomIdx(i,2)-1));
    B.Set_Parm(0,thetamin);
    B.Set_Parm(1,kappa);

    Bvec.push_back(B);
    m.insert(pair<BondedParams, int>(B,i));

    line_as_string.assign(line);
    outlist.push_back(line_as_string);
  }


  if(use_geometry & USE_ANGLES)
    cout << "Using angles from geometry" << endl;

  ub = unique_params(m);
  cout << endl;
  cout << "FOUND " << ub.size() << " UNIQUE ANGULAR TERMS ";
  cout << "OUT OF " << pFF->NAngles() << endl;

  sprintf(line,"; found %d unique angle terms\n",(int)ub.size());    
  itpstream->write(line,strlen(line));

  // print the unique angular terms that have been found
  for(i=0,bit=ub.begin(); bit != ub.end(); ++bit,i++) {
    sprintf(line,"; [%d] \t %-6s %-6s %-6s   %9.4f %9.4f\n",i,
	    pFF->TypeName((*bit).Get_Atom(0)),
	    pFF->TypeName((*bit).Get_Atom(1)),
	    pFF->TypeName((*bit).Get_Atom(2)),
	    (*bit).Get_Parm(0), 
	    (*bit).Get_Parm(1)); 

    cout << line;
    itpstream->write(line,strlen(line));
  }

  for(i=0, it = outlist.begin(); it < outlist.end(); it++, i++) {
    itpstream->write((*it).c_str(), (*it).length() );
    sprintf(line,"type %d\n", bond_type(Bvec[i],ub));
    itpstream->write(line,strlen(line));
  }

  Bvec.clear();
  ub.clear();
  m.clear();
  outlist.clear();

  // DIHEDRALS
  B.Set_NAtoms(4);
  B.Set_NParms(3);

  sprintf(line,"\n[ dihedrals ]\n; proper torsion terms\n");
  itpstream->write(line,strlen(line));

  for(i=0; i<pFF->NTorsions(); i++) {
    int idx1, idx2, idx3, idx4;
    double V = pFF->TorsionParam(i,0);
    double phi0 = pFF->TorsionParam(i,1);      
    int n = (int)pFF->TorsionParam(i,2);
    double nphi0, phi_s;

    idx1 = pFF->TorsionAtomIdx(i,0);
    idx2 = pFF->TorsionAtomIdx(i,1); 
    idx3 = pFF->TorsionAtomIdx(i,2);
    idx4 = pFF->TorsionAtomIdx(i,3);

    nphi0 = n*phi0;

    if(fabs(sin(nphi0*M_PI/180.0)) > 1.0e-3)
      cout << "WARNING!!! nphi0 = " << nphi0 << endl;

    // there is a - sign between the cos(n phi0) of GROMACS and that
    // of UFF. Therefore, we subtract 180 degrees and to make sure
    // phi0 is between -180 and 180
    phi_s = nphi0-180.0;     

    if(use_geometry & USE_DIHEDRALS) {
      OBAtom *a1, *a2, *a3, *a4;
      a1 = mol.GetAtom(idx1);
      a2 = mol.GetAtom(idx2);
      a3 = mol.GetAtom(idx3);
      a4 = mol.GetAtom(idx4);
      phi_s = n*mol.GetTorsion(a1,a2,a3,a4) + 180.0; // angle in degrees
      // make it between 0 and 360
      phi_s -= 360.0*rint(phi_s/360.0);
    }


    sprintf(line,"%-6d %-6d %-6d %-6d  %d  %9.4f  %9.4f  %d ; %-6s %-6s %-6s %-6s ",
	    idx1, idx2, idx3, idx4,
	    1, // proper dihedrals in GROMACS
	    phi_s, // phi_s in degrees
	    V, // kphi in kJ/mol
	    n, // multiplicity
	    pFF->AtomType(pFF->TorsionAtomIdx(i,0)-1),
	    pFF->AtomType(pFF->TorsionAtomIdx(i,1)-1),
	    pFF->AtomType(pFF->TorsionAtomIdx(i,2)-1),
	    pFF->AtomType(pFF->TorsionAtomIdx(i,3)-1));
      
      
    B.Set_Atoms(pFF->AtomTypeAsNumber(pFF->TorsionAtomIdx(i,0)-1),
		pFF->AtomTypeAsNumber(pFF->TorsionAtomIdx(i,1)-1),
		pFF->AtomTypeAsNumber(pFF->TorsionAtomIdx(i,2)-1),
		pFF->AtomTypeAsNumber(pFF->TorsionAtomIdx(i,3)-1));

    B.Set_Parm(0,(double)n);
    B.Set_Parm(1,phi_s);
    B.Set_Parm(2,V);

    Bvec.push_back(B);
    m.insert(pair<BondedParams, int>(B,i));

    line_as_string.assign(line);
    outlist.push_back(line_as_string);      
  }

  if(use_geometry & USE_DIHEDRALS)
    cout << "Using dihedrals from geometry" << endl;

  ub = unique_params(m);
  cout << endl;
  cout << "FOUND " << ub.size() << " UNIQUE DIHEDRAL TERMS ";
  cout << "OUT OF " << pFF->NTorsions() << endl;

  sprintf(line,"; found %d unique dihedral terms\n",(int)ub.size());    
  itpstream->write(line,strlen(line));

  // print the unique dihedral parameters
  for(i=0,bit=ub.begin(); bit != ub.end(); ++bit,i++) {
    sprintf(line,"; [%d] %-6s %-6s %-6s %-6s   %9.4f  %9.4f %9.4f\n",i,
	    pFF->TypeName((*bit).Get_Atom(0)),
	    pFF->TypeName((*bit).Get_Atom(1)),
	    pFF->TypeName((*bit).Get_Atom(2)),
	    pFF->TypeName((*bit).Get_Atom(3)),
	    (*bit).Get_Parm(0),
	    (*bit).Get_Parm(1),
	    (*bit).Get_Parm(2));
    cout << line;
    itpstream->write(line,strlen(line));
  }
    
  for(i=0, it = outlist.begin(); it < outlist.end(); it++, i++) {
    itpstream->write((*it).c_str(), (*it).length() );
    sprintf(line,"type %d\n", bond_type(Bvec[i],ub));
    itpstream->write(line,strlen(line));
  }

  Bvec.clear();
  ub.clear();
  m.clear();
  outlist.clear();

  // IMPROPER DIHEDRALS
  // the functional form of UFF is not available in GROMACS
  // where impropers are only quadratic
  // we will have to appoximate somehow...
  B.Set_NAtoms(4);
  B.Set_NParms(2);

  sprintf(line,"\n[ dihedrals ]\n; inversion terms (improper dihedrals)\n");
  itpstream->write(line,strlen(line));

  for(i=0; i<pFF->NInversions(); i++) {
    int idx1, idx2, idx3, idx4;
    double k = pFF->InversionParam(i,0);
    double c0 = pFF->InversionParam(i,1);
    double c1 = pFF->InversionParam(i,2);
    double c2 = pFF->InversionParam(i,3);

    double csi0;
    double kcsi; 

    if(fabs(c2) < 1.0e-5) {
      csi0 = 0.0;
      kcsi = k;
    } else {
      csi0 = acos(-c1/(4.0*c2))*M_PI/180.0;
      kcsi = (16.0*c2*c2-c1*c1)/(4.0*c2*c2);
      kcsi = k*kcsi;
    }

    /*
      This inversion is OK and has been checked with benzene
      (tends to oscillate on a plane)
     */
    idx1 = pFF->InversionAtomIdx(i,1);
    idx2 = pFF->InversionAtomIdx(i,0); 
    idx3 = pFF->InversionAtomIdx(i,2);
    idx4 = pFF->InversionAtomIdx(i,3);

    sprintf(line,"%-6d %-6d %-6d %-6d  %d  %9.4f  %9.4f ; %-6s %-6s %-6s %-6s ",
            // central atom goes first here
            idx1, idx2, idx3, idx4,
            2, // improper dihedrals in GROMACS
            csi0, // csi_0 in degrees
            kcsi, // kcsi in kJ/mol/rad^2

            pFF->AtomType(pFF->InversionAtomIdx(i,1)-1),
            pFF->AtomType(pFF->InversionAtomIdx(i,0)-1),
            pFF->AtomType(pFF->InversionAtomIdx(i,2)-1),
            pFF->AtomType(pFF->InversionAtomIdx(i,3)-1)); 
    
    line_as_string.assign(line);
    outlist.push_back(line_as_string);
      
    B.Set_Atoms_Inversion(pFF->AtomTypeAsNumber(pFF->InversionAtomIdx(i,1)-1),
			  pFF->AtomTypeAsNumber(pFF->InversionAtomIdx(i,0)-1),
			  pFF->AtomTypeAsNumber(pFF->InversionAtomIdx(i,2)-1),
			  pFF->AtomTypeAsNumber(pFF->InversionAtomIdx(i,3)-1));
    B.Set_Parm(0,csi0);
    B.Set_Parm(1,kcsi);

    Bvec.push_back(B);
    m.insert(pair<BondedParams, int>(B,i));
  }

  ub = unique_params(m);
  cout << endl;
  cout << "FOUND " << ub.size() << " UNIQUE INVERSION TERMS ";
  cout << "OUT OF " << pFF->NInversions() << endl;

  sprintf(line,"; found %d unique inversion terms\n",(int)ub.size());    
  itpstream->write(line,strlen(line));

  // let's build the unique bonds and print them out
  for(i=0,bit=ub.begin(); bit != ub.end(); ++bit,i++) {
    sprintf(line,"; [%d] %-6s %-6s %-6s %-6s     %d     %9.4f  %9.4f\n",i,
	    pFF->TypeName((*bit).Get_Atom(0)),
	    pFF->TypeName((*bit).Get_Atom(1)),
	    pFF->TypeName((*bit).Get_Atom(2)),
	    pFF->TypeName((*bit).Get_Atom(3)),
	    2, 
	    (*bit).Get_Parm(0), 
	    (*bit).Get_Parm(1));
    cout << line;
    itpstream->write(line,strlen(line));
  }

  for(i=0, it = outlist.begin(); it < outlist.end(); it++, i++) {
    itpstream->write((*it).c_str(), (*it).length() );
    sprintf(line,"type %d\n", bond_type(Bvec[i],ub));
    itpstream->write(line,strlen(line));
  }


  //for(i=0; i<pFF->NInversions(); i++)
  // cout << "Inversion " << i << " is of inversion type " << bond_type(Bvec[i],ub) << endl;

  Bvec.clear();
  ub.clear();
  m.clear();


  itpstream->close();

  ifs.close();
  ofs.close();

  cout << "\nPlease cite:" << endl;
  cout << "G. Garberoglio. OBGMX: A web-based generator of GROMACS topologies for" << endl;
  cout << "molecular and periodic systems using the Universal Force Field." << endl;
  cout << "J. Comput. Chem. 33 (2012), 2204-2208. doi: 10.1002/jcc.23049" << endl;
  cout << "\n" << endl;

  return(1);
}

/* obgmx man page*/
/** \page generate GROMACS topologies
*
* \n
* \par SYNOPSIS
*
* \b obgmx [options] \<filename\>
*
* \par DESCRIPTION
*
* The obgmx tool can be used to generate GROMACS topologies
*
* \par OPTIONS
*
* If no filename is given, obgmx will give all options including the
* available forcefields.
*
* \b -ff \<forcefield\>:
*     Select the forcefield \n\n
*
* \par EXAMPLES
*  - View the possible options, including available forcefields: 
*   obgmx
*  - Calculate the energy for the molecule(s) in file test.mol2:
*   obgmx test.mol2
*  - Calculate the energy for the molecule(s) in file test.mol2 using the UFF forcefield:
*   obgmx -ff UFF test.mol2 
*  - Calculate the energy for the molecule(s) in file test.mol2 and print out all individual energy interactions:
*    obgmx -v test.mol2
*
* \par AUTHORS
*
* The obgmx program was contributed by \b Giovanni \b Garberoglio.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 2012 by Giovanni Garberoglio. \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at: http://openbabel.org/ \n
*   The web pages for Open Babel Molecular Mechanics can be found at: 
*   http://openbabel.org/wiki/Molecular_mechanics \n
**/
