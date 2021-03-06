#ifndef PhysicsTools_NanoTuples_PhotonFwd_h
#define PhysicsTools_NanoTuples_PhotonFwd_h

#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include <vector>

typedef pat::PFParticle        Photon;
typedef edm::Ptr<Photon>       PhotonPtr;
typedef std::vector<PhotonPtr> PhotonPtrVector;

#endif
