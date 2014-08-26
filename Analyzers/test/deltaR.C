#include "TLorentzVector.h"

double deltaR(const double px1, const double py1, const double pz1, const double px2, const double py2, const double pz2){
  TLorentzVector p1(px1,py1,pz1,0.0);
  TLorentzVector p2(px2,py2,pz2,0.0);
  double deta = p1.Eta() - p2.Eta();
  double dphi = fabs(p1.Phi() - p2.Phi());
  if(dphi > 3.1415927) dphi = 2*3.1415927 - dphi;
  return sqrt(deta*deta + dphi*dphi);
}
