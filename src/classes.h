//#include <AnalysisDataFormats/CMGTools/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/UserData.h>
#include <vector>
#include <DataFormats/PatCandidates/interface/PFParticle.h>

edm::Ptr<pat::PFParticle> dummy1;
pat::UserHolder<std::vector<edm::Ptr<pat::PFParticle> > > dummy2;
ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > dummy3;
edm::Wrapper<ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2> > > dummy4;
