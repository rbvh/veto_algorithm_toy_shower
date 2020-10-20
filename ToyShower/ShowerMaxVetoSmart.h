#ifndef ShowerMaxVetoSmart_h
#define ShowerMaxVetoSmart_h

#include "ShowerBase.h"

struct ShowerMaxVetoSmart : public Shower {

  ShowerMaxVetoSmart(int seed) : Shower(seed) {}
  
  struct Antenna : public Shower::Antenna {
    Antenna(Shower* showerPtr, int indexi, int indexj, int indexSelf, int indexLeft, int indexRight, double Q2start) : 
      Shower::Antenna(showerPtr, indexi, indexj, indexSelf, indexLeft, indexRight), trialId(-1) {
      
      // Generate a scale
      generate_scale(Q2start);

      // Add it to the queue
      dynamic_cast<ShowerMaxVetoSmart*>(showerPtr)->queue.push(Trial(trialId, indexSelf, Q2trial));
    }


    void generate_scale(double Q2start) {
      double arg = -2*log(showerPtr->rng())/a + pow2(log(Q2start/m2));
      /*
      cout << log(showerPtr->rng()) << " " << pow2(log(Q2start/m2)) << endl;
      if (arg < 0) {
        Q2trial = 0;
        sijTrial = 0;
        sjkTrial = 0;
        cout << "here" << endl;
      }
      */
      Q2trial = m2*exp(-sqrt(arg));
      double zTrial = pow(Q2trial/m2, showerPtr->rng());
      sijTrial = Q2trial/zTrial;
      sjkTrial = m2*zTrial;

      // Set trial id
      trialId = showerPtr->nextCount();
    }

    void print() {
      Shower::Antenna::print();
      cout << "Q2trial = " << Q2trial << endl;
    }

    // Trial variables
    int trialId;
    double Q2trial;
    double sijTrial, sjkTrial;
  };

  struct Trial {
    Trial(int id, int index, double Q2trial): id(id), index(index), Q2trial(Q2trial) {}
    // Unique trial id to track if trial is still valid
    int id;
    // Antenna index
    int index;
    // Trial value
    double Q2trial;
  };

  struct CompareTrial {
    bool operator() (Trial &a1, Trial &a2) {return a1.Q2trial < a2.Q2trial;}
  };

  void run() {
    Shower::reset();

    // Clear antennae
    antennae.clear();
    // Clear queue
    queue = priority_queue<Trial, vector<Trial>, CompareTrial>();

    // Add the first antenna
    antennae.push_back(Antenna(this, 0, 1, 0, -1, -1, s));

    // Run the shower
    double Q2now = s;
    while(true) {
      // Loop until we find a valid trial
      int indexNow, idNow;
      do {
        indexNow = queue.top().index;
        idNow = queue.top().id;
        queue.pop();
      }
      while (antennae[indexNow].trialId != idNow);

      // Record current scale
      Q2now = antennae[indexNow].Q2trial;

      if (Q2now < Q2cut) {break;}

      // Veto
      if (antennae[indexNow].veto(antennae[indexNow].sijTrial, antennae[indexNow].sjkTrial)) {
        // Make a new trial
        antennae[indexNow] = Antenna(this, antennae[indexNow].indexi, antennae[indexNow].indexj, 
            indexNow, antennae[indexNow].indexLeft, antennae[indexNow].indexRight, Q2now);
        continue;
      }

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
  priority_queue<Trial, vector<Trial>, CompareTrial> queue;
};

#endif