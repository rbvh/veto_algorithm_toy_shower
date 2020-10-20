// Basics.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the Rndm, Vec4,
// RotBstMatrix and Hist classes, and some related global functions.

#include "Basics.h"

// Access time information.
#include <ctime>
#include <limits>

//==========================================================================

// Vec4 class.
// This class implements four-vectors, in energy-momentum space.
// (But could also be used to hold space-time four-vectors.)

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Small number to avoid division by zero.
const double Vec4::TINY = 1e-20;

//--------------------------------------------------------------------------

// Rotation (simple).

void Vec4::rot(double thetaIn, double phiIn) {

  double cthe = cos(thetaIn);
  double sthe = sin(thetaIn);
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double tmpx =  cthe * cphi * xx -    sphi * yy + sthe * cphi * zz;
  double tmpy =  cthe * sphi * xx +    cphi * yy + sthe * sphi * zz;
  double tmpz = -sthe *        xx +                cthe *        zz;
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary axis (nz, ny, nz).

void Vec4::rotaxis(double phiIn, double nx, double ny, double nz) {

  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny         *= norm;
  nz         *= norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Azimuthal rotation phi around an arbitrary (3-vector component of) axis.

void Vec4::rotaxis(double phiIn, const Vec4& n) {

  double nx   = n.xx;
  double ny   = n.yy;
  double nz   = n.zz;
  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx         *= norm;
  ny          *=norm;
  nz          *=norm;
  double cphi = cos(phiIn);
  double sphi = sin(phiIn);
  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
  xx          = tmpx;
  yy          = tmpy;
  zz          = tmpz;

}

//--------------------------------------------------------------------------

// Boost (simple).

void Vec4::bst(double betaX, double betaY, double betaZ) {

  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost (simple, given gamma).

void Vec4::bst(double betaX, double betaY, double betaZ, double gamma) {

  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx += prod2 * betaX;
  yy += prod2 * betaY;
  zz += prod2 * betaZ;
  tt = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p.

void Vec4::bst(const Vec4& pIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m.

void Vec4::bst(const Vec4& pIn, double mIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = pIn.xx / pIn.tt;
  double betaY = pIn.yy / pIn.tt;
  double betaZ = pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
  if (beta2 >= 1.) return;
  double gamma = 1. / sqrt(1. - beta2);
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Boost given by a Vec4 p and double m; boost in opposite direction.

void Vec4::bstback(const Vec4& pIn, double mIn) {

  if (abs(pIn.tt) < Vec4::TINY) return;
  double betaX = -pIn.xx / pIn.tt;
  double betaY = -pIn.yy / pIn.tt;
  double betaZ = -pIn.zz / pIn.tt;
  double gamma = pIn.tt / mIn;
  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
  xx          += prod2 * betaX;
  yy          += prod2 * betaY;
  zz          += prod2 * betaZ;
  tt           = gamma * (tt + prod1);

}

//--------------------------------------------------------------------------

// Arbitrary combination of rotations and boosts defined by 4 * 4 matrix.

void Vec4::rotbst(const RotBstMatrix& M) {

  double x = xx; double y = yy; double z = zz; double t = tt;
  tt = M.M[0][0] * t + M.M[0][1] * x + M.M[0][2] * y +  M.M[0][3] * z;
  xx = M.M[1][0] * t + M.M[1][1] * x + M.M[1][2] * y +  M.M[1][3] * z;
  yy = M.M[2][0] * t + M.M[2][1] * x + M.M[2][2] * y +  M.M[2][3] * z;
  zz = M.M[3][0] * t + M.M[3][1] * x + M.M[3][2] * y +  M.M[3][3] * z;

}

//--------------------------------------------------------------------------

// Print a four-vector: also operator overloading with friend.

ostream& operator<<(ostream& os, const Vec4& v) {
  os << fixed << setprecision(3) << " " << setw(9) << v.xx << " "
     << setw(9) << v.yy << " " << setw(9) << v.zz << " " << setw(9)
     << v.tt << " (" << setw(9) << v.mCalc() << ")\n";
  return os;
}

//--------------------------------------------------------------------------

// The invariant mass of two four-vectors.

double m(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
     - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

//--------------------------------------------------------------------------

// The squared invariant mass of two four-vectors.

double m2(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
     - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
  return m2;
}

//--------------------------------------------------------------------------

// The scalar product of two three-vectors.

double dot3(const Vec4& v1, const Vec4& v2) {
  return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;
}

//--------------------------------------------------------------------------

// The cross product of two three-vectors.

Vec4 cross3(const Vec4& v1, const Vec4& v2) {
  Vec4 v;
  v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
  v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
  v.zz = v1.xx * v2.yy - v1.yy * v2.xx; return v;
}


//--------------------------------------------------------------------------

// Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c)

Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c) {
  Vec4 v(0.,0.,0.,0.);
  v.tt =   a.xx*b.yy*c.zz + a.yy*b.zz*c.xx + a.zz*b.xx*c.yy
         - a.xx*b.zz*c.yy - a.zz*b.yy*c.xx - a.yy*b.xx*c.zz;
  v.xx = -(- a.tt*b.yy*c.zz - a.yy*b.zz*c.tt - a.zz*b.tt*c.yy
           + a.tt*b.zz*c.yy + a.zz*b.yy*c.tt + a.yy*b.tt*c.zz);
  v.yy = -(- a.xx*b.tt*c.zz - a.tt*b.zz*c.xx - a.zz*b.xx*c.tt
           + a.xx*b.zz*c.tt + a.zz*b.tt*c.xx + a.tt*b.xx*c.zz);
  v.zz = -(- a.xx*b.yy*c.tt - a.yy*b.tt*c.xx - a.tt*b.xx*c.yy
           + a.xx*b.tt*c.yy + a.tt*b.yy*c.xx + a.yy*b.xx*c.tt);
  return v;
}

//--------------------------------------------------------------------------

// Opening angle between two three-vectors.

double theta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return acos(cthe);
}

//--------------------------------------------------------------------------

// Cosine of the opening angle between two three-vectors.

double costheta(const Vec4& v1, const Vec4& v2) {
  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
    / sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
    * (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
  cthe = max(-1., min(1., cthe));
  return cthe;
}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors.

double phi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors.

double cosphi(const Vec4& v1, const Vec4& v2) {
  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
    (v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Azimuthal angle between two three-vectors around a third.

double phi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return acos(cphi);
}

//--------------------------------------------------------------------------

// Cosine of the azimuthal angle between two three-vectors around a third.

double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
  double nx = n.xx; double ny = n.yy; double nz = n.zz;
  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
  nx *= norm; ny *=norm; nz *=norm;
  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
    (v1s - v1n*v1n) * (v2s - v2n*v2n) ));
  cphi = max(-1., min(1., cphi));
  return cphi;
}

//--------------------------------------------------------------------------

// Distance in cylindrical (y, phi) coordinates.

double RRapPhi(const Vec4& v1, const Vec4& v2) {
  double dRap = abs(v1.rap() - v2.rap());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dRap*dRap + dPhi*dPhi);
}

//--------------------------------------------------------------------------

// Distance in cylindrical (eta, phi) coordinates.

double REtaPhi(const Vec4& v1, const Vec4& v2) {
  double dEta = abs(v1.eta() - v2.eta());
  double dPhi = abs(v1.phi() - v2.phi());
  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
  return sqrt(dEta*dEta + dPhi*dPhi);
}

//--------------------------------------------------------------------------

// Shift four-momenta within pair from old to new masses.
// Note that p1Move and p2Move change values during operation.

bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New) {

  // Standard kinematics variables.
  double sH  = (p1Move + p2Move).m2Calc();
  double r1  = p1Move.m2Calc() / sH;
  double r2  = p2Move.m2Calc() / sH;
  double r3  = m1New * m1New / sH;
  double r4  = m2New * m2New / sH;
  double l12 = sqrtpos(pow2(1. - r1 - r2) - 4. * r1 * r2);
  double l34 = sqrtpos(pow2(1. - r3 - r4) - 4. * r3 * r4);

  // Check that shift operation possible.
  if (sH <= pow2(m1New + m2New) || l12 < Vec4::TINY || l34 < Vec4::TINY)
    return false;

  // Calculate needed shift and apply it.
  double c1  = 0.5 * ( (1. - r1 + r2) * l34 / l12 - (1. - r3 + r4) );
  double c2  = 0.5 * ( (1. + r1 - r2) * l34 / l12 - (1. + r3 - r4) );
  Vec4   pSh = c1 * p1Move - c2 * p2Move;
  p1Move    += pSh;
  p2Move    -= pSh;
  return true;
}

//--------------------------------------------------------------------------

// Create two vectors that are perpendicular to both input vectors.

pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2) {

  // One perpendicular vector from three-dimensional cross-product.
  Vec4 nPerp( cross3(v1,v2) );
  double TINY = std::numeric_limits<double>::epsilon();
  if ( abs(nPerp.pAbs()) < TINY) {
    Vec4 aux;
    if (v1.px() != 0.)      aux.p(v1.yy,v1.px(),v1.pz(),v1.e());
    else if (v1.py() != 0.) aux.p(v1.px(),v1.pz(),v1.py(),v1.e());
    else if (v1.pz() != 0.) aux.p(v1.pz(),v1.py(),v1.px(),v1.e());
    nPerp.p( cross3(v1,aux) );
  }
  nPerp /= abs(nPerp.pAbs());

  // Second perpendicular vector from four-dimensional cross-product.
  Vec4 lPerp( cross4(v1,v2,nPerp) );
  lPerp /= sqrt(abs(lPerp.m2Calc()));
  return make_pair(nPerp,lPerp);
}

//==========================================================================

// RotBstMatrix class.
// This class implements 4 * 4 matrices that encode an arbitrary combination
// of rotations and boosts, that can be applied to Vec4 four-vectors.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Small number to avoid division by zero.
const double RotBstMatrix::TINY = 1e-20;

//--------------------------------------------------------------------------

// Rotate by polar angle theta and azimuthal angle phi.

void RotBstMatrix::rot(double theta, double phi) {

  // Set up rotation matrix.
  double cthe = cos(theta);
  double sthe = sin(theta);
  double cphi = cos(phi);
  double sphi = sin(phi);
  double Mrot[4][4] = {
    {1.,           0.,         0.,          0.},
    {0.,  cthe * cphi,     - sphi, sthe * cphi},
    {0.,  cthe * sphi,       cphi, sthe * sphi},
    {0., -sthe,                0., cthe       } };

  // Rotate current matrix accordingly.
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mrot[i][0] * Mtmp[0][j] + Mrot[i][1] * Mtmp[1][j]
            + Mrot[i][2] * Mtmp[2][j] + Mrot[i][3] * Mtmp[3][j];

}

//--------------------------------------------------------------------------

// Rotate so that vector originally along z axis becomes parallel with p.

void RotBstMatrix::rot(const Vec4& p) {

  double theta = p.theta();
  double phi = p.phi();
  rot(0., -phi);
  rot(theta, phi);

}

//--------------------------------------------------------------------------

// Boost with velocity vector (betaX, betaY, betaZ).

void RotBstMatrix::bst(double betaX, double betaY, double betaZ) {

  // Set up boost matrix.
  double gm = 1. / sqrt( max( TINY, 1. - betaX*betaX - betaY*betaY
    - betaZ*betaZ ) );
  double gf = gm*gm / (1. + gm);
  double Mbst[4][4] = {
    { gm,           gm*betaX,           gm*betaY,          gm*betaZ },
    { gm*betaX, 1. + gf*betaX*betaX, gf*betaX*betaY, gf*betaX*betaZ },
    { gm*betaY, gf*betaY*betaX, 1. + gf*betaY*betaY, gf*betaY*betaZ },
    { gm*betaZ, gf*betaZ*betaX, gf*betaZ*betaY, 1. + gf*betaZ*betaZ } };

  // Boost current matrix correspondingly.
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mbst[i][0] * Mtmp[0][j] + Mbst[i][1] * Mtmp[1][j]
            + Mbst[i][2] * Mtmp[2][j] + Mbst[i][3] * Mtmp[3][j];

}

//--------------------------------------------------------------------------

// Boost so that vector originally at rest obtains same velocity as p.

void RotBstMatrix::bst(const Vec4& p) {
  double betaX = p.px() / p.e();
  double betaY = p.py() / p.e();
  double betaZ = p.pz() / p.e();
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost so vector originally with same velocity as p is brought to rest.

void RotBstMatrix::bstback(const Vec4& p) {
  double betaX = -p.px() / p.e();
  double betaY = -p.py() / p.e();
  double betaZ = -p.pz() / p.e();
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost that transforms p1 to p2, where p1^2 = p2^2 is assumed.

void RotBstMatrix::bst(const Vec4& p1, const Vec4& p2) {
  double eSum = p1.e() + p2.e();
  double betaX = (p2.px() - p1.px()) / eSum;
  double betaY = (p2.py() - p1.py()) / eSum;
  double betaZ = (p2.pz() - p1.pz()) / eSum;
  double fac = 2. / (1. + betaX*betaX + betaY*betaY + betaZ*betaZ);
  betaX *= fac; betaY *= fac; betaZ *= fac;
  bst(betaX, betaY, betaZ);
}

//--------------------------------------------------------------------------

// Boost and rotation that transforms from p1 and p2
// to their rest frame with p1 along +z axis.

void RotBstMatrix::toCMframe(const Vec4& p1, const Vec4& p2) {
  Vec4 pSum = p1 + p2;
  Vec4 dir  = p1;
  dir.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();
  bstback(pSum);
  rot(0., -phi);
  rot(-theta, phi);
}

//--------------------------------------------------------------------------

// Rotation and boost that transforms from rest frame of p1 and p2
// with p1 along +z axis to actual frame of p1 and p2. (Inverse of above.)

void RotBstMatrix::fromCMframe(const Vec4& p1, const Vec4& p2) {
  Vec4 pSum = p1 + p2;
  Vec4 dir  = p1;
  dir.bstback(pSum);
  double theta = dir.theta();
  double phi   = dir.phi();
  rot(0., -phi);
  rot(theta, phi);
  bst(pSum);
}

//--------------------------------------------------------------------------

// Combine existing rotation/boost matrix with another one.

void RotBstMatrix::rotbst(const RotBstMatrix& Mrb) {
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = Mrb.M[i][0] * Mtmp[0][j] + Mrb.M[i][1] * Mtmp[1][j]
            + Mrb.M[i][2] * Mtmp[2][j] + Mrb.M[i][3] * Mtmp[3][j];
}

//--------------------------------------------------------------------------

// Invert the rotation and boost.

void RotBstMatrix::invert() {
  double Mtmp[4][4];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    Mtmp[i][j] = M[i][j];
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = ( (i == 0 && j > 0) || (i > 0 && j == 0) )
      ? - Mtmp[j][i] : Mtmp[j][i];
}

//--------------------------------------------------------------------------

// Reset to diagonal matrix.

void RotBstMatrix::reset() {
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    M[i][j] = (i==j) ? 1. : 0.;
}

//--------------------------------------------------------------------------

// Crude estimate deviation from unit matrix.

double RotBstMatrix::deviation() const {
  double devSum = 0.;
  for (int i = 0; i < 4; ++i)
  for (int j = 0; j < 4; ++j)
    devSum += (i==j) ? abs(M[i][j] - 1.) : abs(M[i][j]);
  return devSum;
}

//--------------------------------------------------------------------------

// Print a rotation and boost matrix: operator overloading with friend.

ostream& operator<<(ostream& os, const RotBstMatrix& M) {
  os << fixed << setprecision(5) << "    Rotation/boost matrix: \n";
  for (int i = 0; i <4; ++i)
    os << setw(10) << M.M[i][0] << setw(10) << M.M[i][1]
       << setw(10) << M.M[i][2] << setw(10) << M.M[i][3] << "\n";
  return os;
}