#ifndef ShowerVetoMaxDumb_h
#define ShowerVetoMaxDumb_h

#include "ShowerBase.h"

struct ShowerVetoMaxDumb : public Shower {

  ShowerVetoMaxDumb(int seed) : Shower(seed) {}

  struct Antenna : public Shower::Antenna {
    Antenna(Shower* showerPtr, int indexi, int indexj, int indexSelf, int indexLeft, int indexRight, double Q2start) : 
      Shower::Antenna(showerPtr, indexi, indexj, indexSelf, indexLeft, indexRight) {
      
      // Generate a scale
      generate_scale(Q2start);
    }


    void generate_scale(double Q2start) {
      if (Iz == 0.) {
        Q2trial = 0;
        sijTrial = 0;
        sjkTrial = 0;
        return;
      }

      Q2trial = Q2start;
      while(true) {
        Q2trial = Q2trial*pow(showerPtr->rng(), 1./a/Iz);
        double zTrial = zMin + showerPtr->rng()*(zPls - zMin);
        sijTrial = m2*exp(zTrial);
        sjkTrial = Q2trial/exp(zTrial);

        if (!veto(sijTrial, sjkTrial)) {break;}
      }
    }

    void print() {
      Shower::Antenna::print();
      cout << "Q2trial = " << Q2trial << endl;
    }

    // Trial variables
    double Q2trial;
    double sijTrial, sjkTrial;
  };

  void run() {
    Shower::reset();

    // Clear antennae
    antennae.clear();

    // Add the first antenna
    antennae.push_back(Antenna(this, 0, 1, 0, -1, -1, s));

    // Run the shower
    double Q2now = s;
    while(true) {

      // Loop over all trials and find the highest
      Q2now = 0;
      int indexNow = -1;
      for (int i=0; i<antennae.size(); i++) {
        if (antennae[i].Q2trial > Q2now) {
          Q2now = antennae[i].Q2trial;
          indexNow = i;
        }
      }

      // Stopping condition
      if (Q2now < Q2cut) {break;}

      // Do the kinematics
      vector<Vec4> pThree;
      double sij = antennae[indexNow].sijTrial;
      double sjk = antennae[indexNow].sjkTrial;
      antennae[indexNow].kinematics(pThree, sij, sjk);

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
      antennae[indexNow] = Antenna(this, indexi, indexj, indexNow, indexLeft, indexNew, Q2now);

      // And add a new one as the right-hand antenna
      antennae.push_back(Antenna(this, indexj, indexk, indexNew, indexNow, indexRight, Q2now));

      // Update links of left and right ones
      if (indexLeft != -1) {
        antennae[indexLeft] = Antenna(this, antennae[indexLeft].indexi, antennae[indexLeft].indexj, 
          indexLeft, antennae[indexLeft].indexLeft, indexNow, Q2now);
      }
      if (indexRight != -1) {
        antennae[indexRight] = Antenna(this, antennae[indexRight].indexi, antennae[indexRight].indexj,
          indexRight, indexNew, antennae[indexRight].indexRight, Q2now);
      }
 
      // Because we're being dumb, generate new scales for all others too
      for (int i=0; i<antennae.size(); i++) {
        if (i!=indexNow && i!=indexNew && i!=indexLeft && i!=indexRight) {
          antennae[i] = Antenna(this, antennae[i].indexi, antennae[i].indexj, 
            i, antennae[i].indexLeft, antennae[i].indexRight, Q2now);
        }
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