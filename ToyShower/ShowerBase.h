#ifndef Showerbase_h
#define Showerbase_h

#include "Basics.h"
#include <queue>
#include <random>
#include <chrono>

using namespace std;

const double s = 1e6;
const double Q2cut = 1;
const double a = 0.1;

//---------------------------------------------------------
//---------------------------------------------------------
//---------------------------------------------------------
struct Shower {
  Shower(int seed): urd(0,1) {
      seed_seq seq{seed, seed, seed, seed, seed};
      gen = mt19937(seq);
  }

  struct Antenna;

  // Reset the event
  void reset() {
    particles.clear();
    counter = 0;
    functionCalls = 0;

    double E = sqrt(s)/2;
    particles.push_back(Vec4(0,0,E,E));
    particles.push_back(Vec4(0,0,-E,E));

    // Record start time
    startTime = chrono::high_resolution_clock::now();
  }

  void finalize() {
    stopTime = chrono::high_resolution_clock::now();
  }

  double rng() {return urd(gen);}

  int nextCount() {counter++; return counter;}

  virtual void run() = 0;

  virtual int getMultiplicity() = 0;

  double getRunningTime() {
    chrono::duration<double> diff = stopTime - startTime;
    return diff.count();
  }

  virtual int getFunctionCalls() {
      return functionCalls;
  }

  void printMomentumSum() {
    Vec4 pSum;
    for (int i=0; i<particles.size(); i++) {
      pSum += particles[i];
    }
    cout << pSum << endl;
  }

  mt19937 gen;
  uniform_real_distribution<double> urd;

  int counter; 

  int functionCalls;
  chrono::time_point<chrono::high_resolution_clock> startTime;
  chrono::time_point<chrono::high_resolution_clock> stopTime;
  

  vector<Vec4> particles;
};

struct Shower::Antenna {
  Antenna(Shower* showerPtr, int indexi, int indexj, int indexSelf, int indexLeft, int indexRight) : 
  showerPtr(showerPtr), indexi(indexi), indexj(indexj), 
  indexSelf(indexSelf), indexLeft(indexLeft), indexRight(indexRight) {
    m2 = (showerPtr->particles[indexi] + showerPtr->particles[indexj]).m2Calc();
    m = sqrt(m2);

    if (4*Q2cut/m2 > 1) {
      zPls = log(0.5);
      zMin = log(0.5);
      Iz = 0;
    }
    else if (Q2cut/m2 < 1e-8) {
      zPls = log(1 - Q2cut/m2);
      zMin = log(Q2cut/m2);
      Iz = zPls - zMin;
    }
    else {
      zPls = log(0.5*(1 + sqrt(1 - 4*Q2cut/m2)));
      zMin = log(0.5*(1 - sqrt(1 - 4*Q2cut/m2)));
      Iz = zPls - zMin;
    }
  }

  void kinematics(vector<Vec4>& pThree, double sij, double sjk) const {
    double sik = m2 - sij - sjk;

    // Set up kinematics in rest frame.
    double Ei = 1/m*(sij + sik)/2.;
    double Ej = 1/m*(sij + sjk)/2.;
    double Ek = 1/m*(sik + sjk)/2.;
    double cosij = (Ei*Ej - sij/2)/Ei/Ej;
    double cosik = (Ei*Ek - sik/2)/Ei/Ek;

    // Protection: num. precision loss for small (ultracollinear) invariants.
    if ( 1-abs(cosij) < 1e-15 ) cosij = cosij > 0 ? 1. : -1.;
    if ( 1-abs(cosik) < 1e-15 ) cosik = cosik > 0 ? 1. : -1.;
    double sinij = (abs(cosij) < 1) ? sqrt(abs(1.0 - pow2(cosij))) : 0.0;
    double sinik = (abs(cosik) < 1) ? sqrt(abs(1.0 - pow2(cosik))) : 0.0;

    // Set momenta in CMz frame (with 1 oriented along positive z axis
    // and event in (x,z) plane).
    Vec4 pi(0.0,0.0,Ei,Ei);
    Vec4 pj(-Ej*sinij, 0.0, Ej*cosij, Ej);
    Vec4 pk(Ek*sinik,  0.0, Ek*cosik, Ek);

    double rMap = sjk/(sij+sjk);
    double rho=sqrt(1.0+4.0*rMap*(1.0-rMap)*sij*sjk/sik/m2);
    double s00=-( (1.0-rho)*m2*sik + 2.0*rMap*sij*sjk ) / 2.0 /
      (m2 - sij);
    double psi = acos(1.0+2.0*s00/(m2-sjk));

    double phi = showerPtr->rng()*2*M_PI;

    // Perform global rotations.
    pi.rot(psi,phi);
    pj.rot(psi,phi);
    pk.rot(psi,phi);

    // Rotate and boost to lab frame.
    RotBstMatrix M;
    M.fromCMframe(showerPtr->particles[indexi], showerPtr->particles[indexj]);
    Vec4 total = showerPtr->particles[indexi] + showerPtr->particles[indexj];

    pi.rotbst(M);
    pj.rotbst(M);
    pk.rotbst(M);

    // Save momenta.
    pThree.resize(0);
    pThree.push_back(pi);
    pThree.push_back(pj);
    pThree.push_back(pk);
  }

  bool veto(double sij, double sjk) {
    // Add to functionCalls 
    showerPtr->functionCalls++;

    // Check phase space
    double sik = m2 - sij - sjk;
    if (sik < 0) {return true;}

    // Veto probability
    double pAccept = sik/m2;

    if (showerPtr->rng() < pAccept) {
      return false;
    }
    else {
      return true;
    }
  }

  void print() {
    cout << "Particles: (" << indexi << ", " << indexj << ") m2 = " << m2 << endl;
    cout << "Indices  : (" << indexLeft << ")-(" << indexSelf << ")-(" << indexRight << ")" << endl;
  }

  Shower* showerPtr;
  // Indices of particles
  int indexi, indexj;
  // Index of this antenna
  int indexSelf;
  // Indices of neighboring antennae
  int indexLeft, indexRight;
  double m2, m;
  double Iz;
  double zPls, zMin;
};

#endif