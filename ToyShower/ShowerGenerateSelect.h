#ifndef GenerateSelect_h
#define GenerateSelect_h

#include "ShowerBase.h"

struct ShowerGenerateSelect : public Shower {

  ShowerGenerateSelect(int seed) : Shower(seed) {}
  struct Antenna : public Shower::Antenna {
    Antenna(Shower* showerPtr, int indexi, int indexj, int indexSelf, int indexLeft, int indexRight, double Q2start) : 
      Shower::Antenna(showerPtr, indexi, indexj, indexSelf, indexLeft, indexRight), weightId(-1) {
      
      // Make a new weight id
      weightId = showerPtr->nextCount();

      // Add it to the queue
      dynamic_cast<ShowerGenerateSelect*>(showerPtr)->queue.push(Weight(weightId, indexSelf, Iz));
    }

    void print() {
      Shower::Antenna::print();
    }

    // Weight id for the queue
    int weightId;
  };

  struct Weight {
    Weight(int id, int index, double weight): id(id), index(index), weight(weight) {}
    // Unique trial id to track if trial is still valid
    int id;
    // Antenna index
    int index;
    // Trial value
    double weight;
  };

  struct CompareWeight {
    bool operator() (Weight &a1, Weight &a2) {return a1.weight < a2.weight;}
  };

  void run() {
    Shower::reset();

    // Clear antennae
    antennae.clear();
    // Clear queue
    queue = priority_queue<Weight, vector<Weight>, CompareWeight>();

    // Add the first antenna
    antennae.push_back(Antenna(this, 0, 1, 0, -1, -1, s));

    // Track the total weight
    double IzSum = antennae[0].Iz;

    // Run the shower
    double Q2now = s;
    while(true) {
      // Generate a new scale
      Q2now = Q2now*pow(rng(), 1./a/IzSum);

      // Stopping condition
      if (Q2now < Q2cut) {break;}

      // Select an element - first pop the queue until we find an Iz that still exists
      int indexMax, idMax;
      while(true) {
        indexMax = queue.top().index;
        idMax = queue.top().id;
        if (antennae[indexMax].weightId == idMax) {break;}
        queue.pop();
      }
     
      double IzMax = antennae[indexMax].Iz;

      // Roulette wheel to select an element
      int indexNow;
      while(true) {
        // Select one at random
        indexNow = rng()*antennae.size();

        // Accept or reject
        if (rng() < antennae[indexNow].Iz/IzMax) {break;}
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

      // Subtract Iz that's getting replaced
      IzSum -= antennae[indexNow].Iz;
      // Reuse the old one as the left-hand antenna
      antennae[indexNow] = Antenna(this, indexi, indexj, indexNow, indexLeft, indexNew, Q2now);
      IzSum += antennae[indexNow].Iz;

      // And add a new one as the right-hand antenna
      antennae.push_back(Antenna(this, indexj, indexk, indexNew, indexNow, indexRight, Q2now));
      IzSum += antennae[indexNew].Iz;

      // Update links of left and right ones
      if (indexLeft != -1) {
        IzSum -= antennae[indexLeft].Iz;
        antennae[indexLeft] = Antenna(this, antennae[indexLeft].indexi, antennae[indexLeft].indexj, 
          indexLeft, antennae[indexLeft].indexLeft, indexNow, Q2now);
        IzSum += antennae[indexLeft].Iz;
      }
      if (indexRight != -1) {
        IzSum -= antennae[indexRight].Iz;
        antennae[indexRight] = Antenna(this, antennae[indexRight].indexi, antennae[indexRight].indexj,
          indexRight, indexNew, antennae[indexRight].indexRight, Q2now);
        IzSum += antennae[indexRight].Iz;
      }
    }

    Shower:finalize();
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
  priority_queue<Weight, vector<Weight>, CompareWeight> queue;
};

#endif