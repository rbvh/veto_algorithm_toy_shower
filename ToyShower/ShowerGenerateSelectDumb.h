#ifndef GenerateSelectDumb_h
#define GenerateSelectDumb_h

#include "ShowerBase.h"

struct ShowerGenerateSelectDumb : public Shower {

  ShowerGenerateSelectDumb(int seed) : Shower(seed) {}
  struct Antenna : public Shower::Antenna {
    Antenna(Shower* showerPtr, int indexi, int indexj, int indexSelf, int indexLeft, int indexRight) : 
      Shower::Antenna(showerPtr, indexi, indexj, indexSelf, indexLeft, indexRight) {
    }

    void print() {
      Shower::Antenna::print();
    }  
  };

  void run() {
    Shower::reset();

    // Clear antennae
    antennae.clear();

    // Add the first antenna
    antennae.push_back(Antenna(this, 0, 1, 0, -1, -1));

    // Run the shower
    double Q2now = s;
    while(true) {

      // Compute Iz by resetting every antenna
      double IzSum = 0;
      for (int i=0; i<antennae.size(); i++) {
        antennae[i] = Antenna(this, antennae[i].indexi, antennae[i].indexj, 
          i, antennae[i].indexLeft, antennae[i].indexRight);
        IzSum += antennae[i].Iz;
      }

      // Generate a new scale
      Q2now = Q2now*pow(rng(), 1./a/IzSum);
    
      // Stopping condition
      if (Q2now < Q2cut) {break;}

      // Select an element the dumb way
      double IzSelect = IzSum*rng();
      double IzCount = 0;
      int indexNow = -1;
      for (int i=0; i<antennae.size(); i++) {
        IzCount += antennae[i].Iz;
        if (IzCount > IzSelect) {
          indexNow = i;
          break;
        }
      }

      // Generate a z, sij, sjk
      double zNow = antennae[indexNow].zMin + rng()*(antennae[indexNow].zPls - antennae[indexNow].zMin);
      double sijNow = antennae[indexNow].m2*exp(zNow);
      double sjkNow = Q2now/exp(zNow);

      // Now veto in that channel
      if (antennae[indexNow].veto(sijNow, sjkNow)) {
        continue;
      }
      
      // Do the kinematics
      vector<Vec4> pThree;
      antennae[indexNow].kinematics(pThree, sijNow, sjkNow);

      // Overwrite two momenta 
      int indexi = antennae[indexNow].indexi;
      int indexk = antennae[indexNow].indexj;
      particles[indexi] = pThree[0];
      particles[indexk] = pThree[2];

      // Add the new one
      particles.push_back(pThree[1]);
      int indexj = particles.size()-1;

      // Update antennae
      int indexLeft  = antennae[indexNow].indexLeft;
      int indexRight = antennae[indexNow].indexRight;
      int indexNew = antennae.size();
      // Reuse the old one as the left-hand antenna
      antennae[indexNow] = Antenna(this, indexi, indexj, indexNow, indexLeft, indexNew);

      // And add a new one as the right-hand antenna
      antennae.push_back(Antenna(this, indexj, indexk, indexNew, indexNow, indexRight));

      // Update links of left and right ones
      if (indexLeft != -1) {
        antennae[indexLeft] = Antenna(this, antennae[indexLeft].indexi, antennae[indexLeft].indexj, 
          indexLeft, antennae[indexLeft].indexLeft, indexNow);
      }
      if (indexRight != -1) {
        antennae[indexRight] = Antenna(this, antennae[indexRight].indexi, antennae[indexRight].indexj,
          indexRight, indexNew, antennae[indexRight].indexRight);
      }
    }

    Shower::finalize();
  }

  int getMultiplicity() {
    return antennae.size()+1;
  }

  void print() {
    cout << "---------------------------" << endl;
    for (int i=0; i<antennae.size(); i++) {
      antennae[i].print();
      if (i!=antennae.size()-1) {cout << endl;}
    }
    cout << "---------------------------" << endl;
  }

  void checkAntennaIndices() {
    // Find the left-hand antenna
    cout << "Checking" << endl;
    int iNow = -1;
    for (int i=0; i<antennae.size(); i++) {
      if (antennae[i].indexLeft == -1) {
        iNow = i;
        break;
      }
    }
    if (iNow == -1) {
      cout << "Did not find left-hand index" << endl;
      exit(0);
    }

    // Follow the chain
    int iNowLeft = antennae[iNow].indexLeft;
    int iNowRight = antennae[iNow].indexRight;
    while (iNowRight != -1) {
      // Check internal index
      if (iNow != antennae[iNow].indexSelf) {
        cout << "Internal index wrong" << endl;
        exit(0);
      }

      // Check right side
      if (antennae[iNowRight].indexLeft != iNow) {
        cout << "Left side wrong" << endl;
        exit(0);
      }
      if (iNowLeft != -1) {
        if (antennae[iNowLeft].indexRight != iNow) {
          cout << "Right side wrong" << endl;
        }
      }

      // Shift
      iNow = iNowRight;
      iNowRight = antennae[iNow].indexRight;
      iNowLeft = antennae[iNow].indexLeft;
    }
    cout << "Passed" << endl;
  }

  vector<Antenna> antennae;
};

#endif