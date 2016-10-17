//////////////////////////////////////////////////////////////
// Read in binary preprocessed data and convert to ROOT files
/// G. Dedes, Department of Medical Physics, LMU
/// 25.08.2015
//////////////////////////////////////////////////////////////

#include "pctlluconverter_ggo.h"

// Standard libs
#include <iostream>
#include <fstream>
#include <vector>

// ROOT headers
#include "TApplication.h"
#include "TFile.h"
#include "TNtuple.h"

#include <rtkGgoFunctions.h>
#include <rtkMacro.h>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char **argv)
{
  GGO(pctlluconverter, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

  std::ifstream infile;
  infile.open(args_info.input_arg, std::ios::binary | std::ios::in);

  // reading header
  std::cout << "Reading header..." << std::endl << std::endl;

  char magicNumber[4];
  infile.read((char*)&magicNumber, sizeof(magicNumber));
  std::cout << "Magic number: " << magicNumber[0] << magicNumber[1] << magicNumber[2] << magicNumber[3] << std::endl;

  int formaVersionID;
  infile.read((char*)&formaVersionID, sizeof(formaVersionID));
  std::cout << "Format version ID: " << formaVersionID << std::endl;

  int numberOfEvents;
  infile.read((char*)&numberOfEvents, sizeof(numberOfEvents));
  std::cout << "Number of events: " << numberOfEvents << std::endl;

  float projectionAngle;
  infile.read((char*)&projectionAngle, sizeof(projectionAngle));
  std::cout << "Projection angle [deg]: " << projectionAngle << std::endl;

  float beamEnergy;
  infile.read((char*)&beamEnergy, sizeof(beamEnergy));
  std::cout << "Beam energy [MeV]: " << beamEnergy << std::endl;

  int acqGendate;
  infile.read((char*)&acqGendate, sizeof(acqGendate));
  std::cout << "Acq/Gen date: " << acqGendate << std::endl;

  int preProcessdate;
  infile.read((char*)&preProcessdate, sizeof(preProcessdate));
  std::cout << "Preprocess date: " << preProcessdate << std::endl;

  int phantomNameLength;
  infile.read((char*)&phantomNameLength, sizeof(phantomNameLength));
  std::cout << "Phantom name length: " << phantomNameLength << std::endl;

  char *phantomName = new char[phantomNameLength+1];
  infile.read(phantomName, sizeof(char)*phantomNameLength);
  std::cout << "Phantom name: " << phantomName << std::endl;
  delete phantomName;

  int dataSourceNameLength;
  infile.read((char*)&dataSourceNameLength, sizeof(dataSourceNameLength));
  std::cout << "Data source name length: " << dataSourceNameLength << std::endl;

  char *dataSourceName = new char[dataSourceNameLength+1];
  infile.read(dataSourceName, sizeof(char)*dataSourceNameLength);
  std::cout << "Data source name: " << dataSourceName << std::endl;
  delete dataSourceName;

  int personNameLength;
  infile.read((char*)&personNameLength, sizeof(personNameLength));
  std::cout << "Person name length: " << personNameLength << std::endl;

  char *personName = new char[personNameLength+1];
  infile.read(personName, sizeof(char)*personNameLength);
  std::cout << "Prepared by: " << personName << std::endl;
  delete personName;

  std::cout << "End of header" << std::endl;

  // t coordinate
  float t0, t1, t2, t3;
  std::vector <float > t0_vec, t1_vec, t2_vec, t3_vec;
  for(int t0_index=0; t0_index<numberOfEvents; t0_index++)
    {
    infile.read((char*)&t0, sizeof(t0));
    t0_vec.push_back(t0);
    }
  for(int t1_index=0; t1_index<numberOfEvents; t1_index++)
    {
    infile.read((char*)&t1, sizeof(t1));
    t1_vec.push_back(t1);
    }
  for(int t2_index=0; t2_index<numberOfEvents; t2_index++)
    {
    infile.read((char*)&t2, sizeof(t2));
    t2_vec.push_back(t2);
    }
  for(int t3_index=0; t3_index<numberOfEvents; t3_index++)
    {
    infile.read((char*)&t3, sizeof(t3));
    t3_vec.push_back(t3);
    }

  // v coordinate
  float v0, v1, v2, v3;
  std::vector <float > v0_vec, v1_vec, v2_vec, v3_vec;
  for(int v0_index=0; v0_index<numberOfEvents; v0_index++)
    {
    infile.read((char*)&v0, sizeof(v0));
    v0_vec.push_back(v0);
    }
  for(int v1_index=0; v1_index<numberOfEvents; v1_index++)
    {
    infile.read((char*)&v1, sizeof(v1));
    v1_vec.push_back(v1);
    }
  for(int v2_index=0; v2_index<numberOfEvents; v2_index++)
    {
    infile.read((char*)&v2, sizeof(v2));
    v2_vec.push_back(v2);
    }
  for(int v3_index=0; v3_index<numberOfEvents; v3_index++)
    {
    infile.read((char*)&v3, sizeof(v3));
    v3_vec.push_back(v3);
    }

  // u coordinate
  float u0, u1, u2, u3;
  std::vector <float > u0_vec, u1_vec, u2_vec, u3_vec;
  for(int u0_index=0; u0_index<numberOfEvents; u0_index++)
    {
    infile.read((char*)&u0, sizeof(u0));
    u0_vec.push_back(u0);
    }
  for(int u1_index=0; u1_index<numberOfEvents; u1_index++)
    {
    infile.read((char*)&u1, sizeof(u1));
    u1_vec.push_back(u1);
    }
  for(int u2_index=0; u2_index<numberOfEvents; u2_index++)
    {
    infile.read((char*)&u2, sizeof(u2));
    u2_vec.push_back(u2);
    }
  for(int u3_index=0; u3_index<numberOfEvents; u3_index++)
    {
    infile.read((char*)&u3, sizeof(u3));
    u3_vec.push_back(u3);
    }

  // wepl
  float wepl;
  std::vector <float > wepl_vec;
  for(int wepl_index=0; wepl_index<numberOfEvents; wepl_index++)
    {
    infile.read((char*)&wepl, sizeof(wepl));
    wepl_vec.push_back(wepl);
    }

  std::cout << "Writing into ROOT file: " << args_info.output_arg << std::endl;
  TFile *my_rootfile = new TFile(args_info.output_arg,"RECREATE");
  TNtuple *my_Ntuple = new TNtuple("recoENTRY", "recoENTRY","v_hit0:t_hit0:u_hit0:v_hit1:t_hit1:u_hit1:v_hit2:t_hit2:u_hit2:v_hit3:t_hit3:u_hit3:calculated_WEPL:projection_angle");
  float v_hit0, t_hit0, u_hit0, v_hit1, t_hit1, u_hit1, v_hit2, t_hit2, u_hit2, v_hit3, t_hit3, u_hit3, calculated_WEPL, projection_angle;
  for(int numberOfEvents_index=0; numberOfEvents_index<numberOfEvents; numberOfEvents_index++)
    {
    v_hit0 = v0_vec[numberOfEvents_index];
    t_hit0 = t0_vec[numberOfEvents_index];
    u_hit0 = u0_vec[numberOfEvents_index];
    v_hit1 = v1_vec[numberOfEvents_index];
    t_hit1 = t1_vec[numberOfEvents_index];
    u_hit1 = u1_vec[numberOfEvents_index];
    v_hit2 = v2_vec[numberOfEvents_index];
    t_hit2 = t2_vec[numberOfEvents_index];
    u_hit2 = u2_vec[numberOfEvents_index];
    v_hit3 = v3_vec[numberOfEvents_index];
    t_hit3 = t3_vec[numberOfEvents_index];
    u_hit3 = u3_vec[numberOfEvents_index];
    calculated_WEPL = wepl_vec[numberOfEvents_index];
    projection_angle = projectionAngle;
    if (calculated_WEPL>-40. && calculated_WEPL<300.)
      my_Ntuple->Fill(v_hit0, t_hit0, u_hit0, v_hit1, t_hit1, u_hit1, v_hit2, t_hit2, u_hit2, v_hit3, t_hit3, u_hit3, calculated_WEPL, projection_angle);
    }
  // Write data to ROOT file
  my_Ntuple->Write();
  my_rootfile->Close();

  std::cout << "Done" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  return 0;
}
