/* \class ClassicSVfitInterface
**
** This class is a copy of the class SVfitInterface, adapted to be used
** with the SVfit/ClassicSVfit package.
**
** This class provides an interface to the SVfit standalone algorithm
** for the computation of the SVfit mass of the lepton pair candidates.
** 
** The decay mode (e, mu, tauh) of each lepton is the pair is asserted
** from the pdgId associated and is used in the algorithm.
**
** input type is reco::CompositeCandidate for each lepton
** that is coverted to TLorentzVector to be passed to the algorithm
**
** output type is pat::CompositeCandidate, i.e. the original pairs
** plus some userfloats containing the SVfit mass and MET px, px, pt, phi. 
**  
** \date:    10 November 2017
** \author:  F.Brivio (MIB)
*/

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/stream/EDProducer.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include "FWCore/Utilities/interface/StreamID.h"

#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
#include <DataFormats/Common/interface/Ptr.h>
#include <DataFormats/Common/interface/PtrVector.h>
#include <DataFormats/NanoAOD/interface/FlatTable.h>

#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h>
#include <TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrand.h>
#include <TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h>
#include <TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h>

#include <PhysicsTools/NanoTuples/interface/CutSet.h>
#include <PhysicsTools/NanoTuples/interface/LeptonIsoHelper.h>
#include <PhysicsTools/NanoTuples/interface/DaughterDataHelpers.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>
#include <cmath>

using namespace edm;
using namespace std;
using namespace reco;
using namespace classic_svFit;
using METUncertainty = pat::MET::METUncertainty;

// ------------------------------------------------------------------

class ClassicSVfitInterface : public edm::stream::EDProducer<> {
 public:
  /// Constructor
  explicit ClassicSVfitInterface(const edm::ParameterSet&);

  /// Destructor
  ~ClassicSVfitInterface(){};

 private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  classic_svFit::MeasuredTauLepton::kDecayType GetDecayTypeFlag (int pdgId);
  bool Switch (classic_svFit::MeasuredTauLepton::kDecayType type1, double pt1, float l1_iso, classic_svFit::MeasuredTauLepton::kDecayType type2, double pt2, float l2_iso);
  double GetMass (classic_svFit::MeasuredTauLepton::kDecayType type, double candMass);
  bool IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2); // if true, compute SVFit

  edm::EDGetTokenT<View<reco::CompositeCandidate> > theCandidateTag;
  // std::vector <edm::EDGetTokenT<View<pat::MET> > > vtheMETTag; // polymorphism of view --> good for reco::PFMET and pat::MET! 
  edm::EDGetTokenT<View<pat::MET> > theMETTag;
  edm::EDGetTokenT<double> theSigTag;
  edm::EDGetTokenT<math::Error<2>::type> theCovTag;
  bool _usePairMET;
  bool _debug;
  TFile* inputFile_visPtResolution_;
  const std::string SVFitName_; 
 
  // 6,7,8 are expected to be unused
  enum pairType {
    kMuHad  = 0,
    kEHad   = 1,
    kHadHad = 2,
    kMuMu   = 3,
    kEE     = 4,
    kEMu    = 5,
    kEEPrompt = 6, // prompt Z->ee/mumu decays
    kMuMuPrompt = 7,
    kOther  = 8 // for e.g. h->bb
  };

};

// ------------------------------------------------------------------


ClassicSVfitInterface::ClassicSVfitInterface(const edm::ParameterSet& iConfig):
theCandidateTag(consumes<View<reco::CompositeCandidate> >(iConfig.getParameter<InputTag>("srcPairs"))),
theMETTag(consumes<View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
theSigTag(consumes<double>(iConfig.getParameter<edm::InputTag>("srcSig"))),
theCovTag(consumes<math::Error<2>::type>(iConfig.getParameter<edm::InputTag>("srcCov"))),
SVFitName_(iConfig.getParameter<std::string>("SVFitName"))
{
  //theCandidateTag = iConfig.getParameter<InputTag>("srcPairs");
  _usePairMET = iConfig.getParameter<bool>("usePairMET");
  _debug = iConfig.getParameter<bool>("debug");
  
  // const std::vector<edm::InputTag>& inMET = iConfig.getParameter<std::vector<edm::InputTag> >("srcMET");
  // for (std::vector<edm::InputTag>::const_iterator it = inMET.begin(); it != inMET.end(); ++it)
  // {      
  //   // vtheMETTag.emplace_back(consumes<edm::View<reco::MET> >(*it) );
  //   vtheMETTag.emplace_back(consumes<edm::View<pat::MET> >(*it) );
  // }

  /*edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  TH1::AddDirectory(false);  
  inputFile_visPtResolution_ = new TFile(inputFileName_visPtResolution.fullPath().data());*/

  produces<nanoaod::FlatTable>(SVFitName_);
  produces<nanoaod::FlatTable>(SVFitName_+"MET");
  produces<edm::PtrVector<reco::Candidate> >();
  
}

void ClassicSVfitInterface::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  auto selCandPf = std::make_unique<PtrVector<reco::Candidate>>();

  // Get lepton pairs
  Handle<View<reco::CompositeCandidate> > pairHandle;
  iEvent.getByToken(theCandidateTag, pairHandle);
  
  unsigned int pairNumber = pairHandle->size();
  unsigned int metNumber = 0;

  // MET class type changes if using MVA MEt or 'ordinary' MEt
  
  // create handle -- view makes possible to use base class type reco::MET
  // but now everything is produced as pat::MET
  Handle<View<pat::MET> > METHandle;

  // Handle<View<reco::PFMET> > METHandle_PfMET;
  // Handle<View<pat::MET> >    METHandle_PatMET;
  
  // intialize MET
  double METx = 0.;
  double METy = 0.;
  float significance = -999.;
  std::vector<float> METxVec, METyVec, uncorrMETx, uncorrMETy, significanceVec, covMET00, covMET10, covMET01, covMET11, isValid, isInteresting; 
  TMatrixD covMET(2, 2);

  iEvent.getByToken(theMETTag, METHandle);

  // initialize MET once if not using PairMET
  if (!_usePairMET)
  {   
     metNumber = METHandle->size();
     if (metNumber != 1)     
        edm::LogWarning("pfMetHasNotSizeOne") << "(ClassicSVfitInterface) Warning! Using single pf MEt, but input MEt collection size is different from 1"
                                                           << "   --> using MET entry num. 0";
     const pat::MET& patMET = (*METHandle)[0];
     METx = patMET.px();
     METy = patMET.py();
     METxVec.push_back(METx);
     METyVec.push_back(METy);

     Handle<double> significanceHandle;
     Handle<math::Error<2>::type> covHandle;

     iEvent.getByToken (theSigTag, significanceHandle);
     iEvent.getByToken (theCovTag, covHandle);
     
     covMET[0][0] = (*covHandle)(0,0);
     covMET[1][0] = (*covHandle)(1,0);
     covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
     covMET[1][1] = (*covHandle)(1,1);
     covMET00.push_back(covMET[0][0]);
     covMET10.push_back(covMET[1][0]);
     covMET01.push_back(covMET[0][1]);
     covMET11.push_back(covMET[1][1]);

     significance = (float) (*significanceHandle);
     significanceVec.push_back(significance);
     // protection against singular matrices
     if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
        edm::LogWarning("SingularCovarianceMatrix") << "(ClassicSVfitInterface) Warning! Input covariance matrix is singular"
                                                    << "   --> SVfit algorithm will probably crash...";

     if(_debug){      
       cout << " -------- CLASSIC SVIFT MET ---------" << endl;
       cout << "MET       : " << patMET.px() << " / " << patMET.py() << endl;
       cout << " -------- ----------------- ---------" << endl;
     }
  }
  

  // Output collection
  std::vector<float> SVfitMass, SVfitTransverseMass, SVpt, SVeta, SVphi, SVMETRho, SVMETPhi, channel, tau1pt, tau1eta, tau1phi, tau1mass, tau1decaymode, tau1pdgid, tau2pt, tau2eta, tau2phi, tau2mass, tau2decaymode, tau2pdgid;

  // loop on all the pairs
  for (unsigned int i = 0; i < pairNumber; ++i)
  {    
    //Get Candidate names
    edm::Ptr<reco::Candidate> c = pairHandle->ptrAt(i);
    selCandPf->push_back(c);

    // Get the pair and the two leptons composing it
    const CompositeCandidate& pairBuf = (*pairHandle)[i];
    pat::CompositeCandidate pair(pairBuf);
    
    const Candidate *l1 = pair.daughter(0);
    const Candidate *l2 = pair.daughter(1);

    classic_svFit::MeasuredTauLepton::kDecayType l1Type = GetDecayTypeFlag (l1->pdgId());
    classic_svFit::MeasuredTauLepton::kDecayType l2Type = GetDecayTypeFlag (l2->pdgId());
    double mass1 = GetMass (l1Type, l1->mass());
    double mass2 = GetMass (l2Type, l2->mass());
   
    int decay1 = -1;
    int decay2 = -1;
    if (l1Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) decay1 = (int)(userdatahelpers::getUserFloat(l1,"decayMode"));
    if (l2Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) decay2 = (int)(userdatahelpers::getUserFloat(l2,"decayMode"));
    
    float l1_iso = -1.;
    float l2_iso = -1.;
    if (l1Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) l1_iso = userdatahelpers::getUserFloat(l1,"byDeepTau2017v2p1VSjetraw");
    if (l2Type == classic_svFit::MeasuredTauLepton::kTauToHadDecay) l2_iso = userdatahelpers::getUserFloat(l2,"byDeepTau2017v2p1VSjetraw");

    bool swi = Switch (l1Type, l1->pt(), l1_iso, l2Type, l2->pt(), l2_iso);
  
    if (_usePairMET)
    {
      // iEvent.getByToken(theMETTag, METHandle);

      // const PFMET* pfMET = (PFMET*) &((*METHandle)[0]) ; // all this to transform the type of the pointer!
      const pat::MET* patMET = &((*METHandle)[i]);
      const reco::METCovMatrix& covMETbuf = patMET->getSignificanceMatrix();
      significanceVec.push_back((float) patMET->significance());

      METx = patMET->px();
      METy = patMET->py();
      METxVec.push_back(METx);
      METyVec.push_back(METy);

      uncorrMETx.push_back(( patMET->hasUserFloat("uncorrPx") ) ? patMET->userFloat("uncorrPx") : -999);
      uncorrMETy.push_back(( patMET->hasUserFloat("uncorrPy") ) ? patMET->userFloat("uncorrPy") : -999);

      covMET[0][0] = covMETbuf(0,0);
      covMET[1][0] = covMETbuf(1,0);
      covMET[0][1] = covMETbuf(0,1);
      covMET[1][1] = covMETbuf(1,1);

      covMET00.push_back(covMETbuf(0,0));
      covMET10.push_back(covMETbuf(1,0));
      covMET01.push_back(covMETbuf(0,1));
      covMET11.push_back(covMETbuf(1,1));

      // protection against singular matrices
      if (covMET[0][0] == 0 && covMET[1][0] == 0 && covMET[0][1] == 0 && covMET[1][1] == 0)
      {
          edm::LogWarning("SingularCovarianceMatrix") << "(ClassicSVfitInterface) Warning! Input covariance matrix is singular"
                                                    << "   --> SVfit algorithm will probably crash...";
      }
    }

    // prepare tau nominal, up, down candidates            
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;


    // set lepton vector, ordered for SVfit
    if (swi)  // 2 first, 1 second (switch)
    {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));
    }

    else
    {
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l1Type, l1->pt(), l1->eta(), l1->phi(), mass1, decay1 ));
      measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(l2Type, l2->pt(), l2->eta(), l2->phi(), mass2, decay2 ));
    }

    // define algorithm (set the debug level to 3 for testing)
    unsigned int verbosity = 0;
    
    // assessing pair type
    int apdg1 = abs(l1->pdgId());
    int apdg2 = abs(l2->pdgId());

    int nmu = 0;
    int nele = 0;
    int ntau = 0;

    if (apdg1 == 13) nmu++;
    if (apdg1 == 11) nele++;
    if (apdg1 == 15) ntau++;

    if (apdg2 == 13) nmu++;
    if (apdg2 == 11) nele++;
    if (apdg2 == 15) ntau++;

    pairType pType = kOther;
    if (nmu == 1 && nele == 0 && ntau == 1) pType = kMuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) pType = kEHad;
    if (nmu == 0 && nele == 0 && ntau == 2) pType = kHadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) pType = kMuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) pType = kEE;
    if (nmu == 1 && nele == 1 && ntau == 0) pType = kEMu;

    // Define the k factor
    double kappa; // use 3 for emu, 4 for etau and mutau, 5 for tautau channel
    if      (pType == kMuHad ) kappa = 4.;  // mutau
    else if (pType == kEHad  ) kappa = 4.;  // etau
    else if (pType == kHadHad) kappa = 5.;  // tautau
    else                       kappa = 3.;  // ee, emu, mumu
    
    //selCandPf->push_back(c);
    // only run SVfit if taus are passing discriminator, skip mumu and ee pairs, apply very loose quality cuts on objects
    // if (isGoodDR && GoodPairFlag)
    if (IsInteresting(l1, l2))
    {
      ClassicSVfit algo(verbosity);
      algo.addLogM_fixed(false, kappa);
      algo.addLogM_dynamic(false);
      //algo.setLikelihoodFileName("testClassicSVfit.root"); //ROOT file to store histograms of di-tau pT, eta, phi, mass and transverse mass, comment if you don't want it
      //algo.shiftVisPt(true, inputFile_visPtResolution_); //not in Classic_svFit

      if(_debug){      
        cout << "--- SVFit Input Debug ---" << endl;
        cout << "pType     = " << pType << endl;
        cout << "lep1 pt   = " << measuredTauLeptons.at(0).pt() << endl;
        cout << "lep1 eta  = " << measuredTauLeptons.at(0).eta() << endl;
        cout << "lep1 phi  = " << measuredTauLeptons.at(0).phi() << endl;
        cout << "lep1 mass = " << measuredTauLeptons.at(0).mass() << endl;
        cout << "lep1 dm   = " << measuredTauLeptons.at(0).decayMode() << endl;
        cout << "lep1 type = " << measuredTauLeptons.at(0).type() << endl;
        cout << "lep2 pt   = " << measuredTauLeptons.at(1).pt() << endl;
        cout << "lep2 eta  = " << measuredTauLeptons.at(1).eta() << endl;
        cout << "lep2 phi  = " << measuredTauLeptons.at(1).phi() << endl;
        cout << "lep2 mass = " << measuredTauLeptons.at(1).mass() << endl;
        cout << "lep2 dm   = " << measuredTauLeptons.at(1).decayMode() << endl;
        cout << "lep2 type = " << measuredTauLeptons.at(1).type() << endl;
        cout << "METx      = " << METx << endl;
        cout << "METy      = " << METy << endl;
        cout << "covMET00  = " << covMET[0][0]<<endl;
        cout << "covMET01  = " << covMET[0][1]<<endl;
        cout << "covMET10  = " << covMET[1][0]<<endl;
        cout << "covMET11  = " << covMET[1][1]<<endl;
        if (measuredTauLeptons.at(0).type() == classic_svFit::MeasuredTauLepton::kTauToHadDecay && measuredTauLeptons.at(1).type() == classic_svFit::MeasuredTauLepton::kTauToHadDecay)
        {
          if (swi)
          {
            cout << "iso1 = " << l2_iso<< endl;
            cout << "iso2 = " << l1_iso << endl;
          }
          else
          {
            cout << "iso1 = " << l1_iso<< endl;
            cout << "iso2 = " << l2_iso << endl;
          }
        }
      }

      algo.integrate(measuredTauLeptons, METx, METy, covMET);
      
      if ( algo.isValidSolution() )
      {
        SVfitMass.push_back((float)static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getMass()); // full mass of tau lepton pair in units of GeV
        SVfitTransverseMass.push_back((float)static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getTransverseMass());
        SVpt.push_back((float)static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPt());
        SVeta.push_back((float)static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getEta());
        SVphi.push_back((float)static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(algo.getHistogramAdapter())->getPhi());
        
        if(_debug){
          cout << "--- SVFit Output Debug ---" << endl;
          cout << "SVfitMass           = " << SVfitMass.back() << endl;
          cout << "SVfitTransverseMass = " << SVfitTransverseMass.back() << endl;
          cout << "SVpt 	             = " << SVpt.back() << endl;
          cout << "SVeta	             = " << SVeta.back() << endl;
          cout << "SVphi	             = " << SVphi.back() << endl;
        }

        channel.push_back((float)pType);
        tau1pt.push_back(measuredTauLeptons.at(0).pt());
        tau1eta.push_back(measuredTauLeptons.at(0).eta());
        tau1phi.push_back(measuredTauLeptons.at(0).phi());
        tau1mass.push_back(measuredTauLeptons.at(0).mass());
        tau1decaymode.push_back(measuredTauLeptons.at(0).decayMode());
        tau1pdgid.push_back(measuredTauLeptons.at(0).type());
        tau2pt.push_back(measuredTauLeptons.at(1).pt());
        tau2eta.push_back(measuredTauLeptons.at(1).eta());
        tau2phi.push_back(measuredTauLeptons.at(1).phi());
        tau2mass.push_back(measuredTauLeptons.at(1).mass());
        tau2decaymode.push_back(measuredTauLeptons.at(1).decayMode());
        tau2pdgid.push_back(measuredTauLeptons.at(1).type());
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau1(l1->pt(), l1->eta(), l1->phi(), mass1);
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredTau2(l2->pt(), l2->eta(), l2->phi(), mass2);
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> measuredDiTauSystem = measuredTau1 + measuredTau2;
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> fittedDiTauSystem(SVpt.back(), SVeta.back(), SVphi.back(), SVfitMass.back());
        Vector fittedMET = (fittedDiTauSystem.Vect() - measuredDiTauSystem.Vect());
        SVMETRho.push_back((float)fittedMET.Rho());
        SVMETPhi.push_back((float)fittedMET.Phi());
      }
      else{
        SVfitMass.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVfitTransverseMass.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVpt.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVeta.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVphi.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVMETRho.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        SVMETPhi.push_back(-111); // -111: SVfit failed (cfr: -999: SVfit not computed)
        channel.push_back(-111);
        tau1pdgid.push_back(-999);
        tau1pt.push_back(-111);
        tau1eta.push_back(-111);
        tau1phi.push_back(-111);
        tau1mass.push_back(-111);
        tau1decaymode.push_back(-111);
        tau2pdgid.push_back(-999);
        tau2pt.push_back(-111);
        tau2eta.push_back(-111);
        tau2phi.push_back(-111);
        tau2mass.push_back(-111);
        tau2decaymode.push_back(-111);
      }
    } else {
      SVfitMass.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVfitTransverseMass.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVpt.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVeta.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVphi.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVMETRho.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      SVMETPhi.push_back(-999); // -111: SVfit failed (cfr: -999: SVfit not computed)
      channel.push_back(-999);
      tau1pdgid.push_back(-999);
      tau1pt.push_back(-999);
      tau1eta.push_back(-999);
      tau1phi.push_back(-999);
      tau1mass.push_back(-999);
      tau1decaymode.push_back(-999);
      tau2pdgid.push_back(-999);
      tau2pt.push_back(-999);
      tau2eta.push_back(-999);
      tau2phi.push_back(-999);
      tau2mass.push_back(-999);
      tau2decaymode.push_back(-999);
    } // end of quality checks IF
  }// end for loop over pairs
  int _isValid = 0;
  int _isInteresting = 0;
  for(long unsigned int i = 0; i != SVfitMass.size(); i++){
    if(SVfitMass[i] >= 0){
      _isValid = 1;
      _isInteresting = 1;
      break;
    } else if(SVfitMass[i] >= -200){
      _isInteresting = 1;
    }
  }
  isValid.push_back(_isValid);
  isInteresting.push_back(_isInteresting);

  auto candTable = std::make_unique<nanoaod::FlatTable>(selCandPf->size(), SVFitName_, false); 
  candTable->setDoc("Save SV Fit candidate reconstruction");
  // add user floats: SVfit mass, met properties, etc..  
  candTable->addColumn<float>("Mass", 		SVfitMass 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("TransverseMass", SVfitTransverseMass, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("Pt", 		SVpt 		, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("Eta", 		SVeta 		, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("Phi", 		SVphi 		, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("METRho", 	SVMETRho 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("METPhi", 	SVMETPhi 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("channel", 	channel 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1Mass", 	tau1mass 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1pdgId", 	tau1pdgid 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1Pt", 	tau1pt 		, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1Eta", 	tau1eta 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1Phi", 	tau1phi 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau1DM", 	tau1decaymode 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2Mass", 	tau2mass 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2pdgId", 	tau2pdgid 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2Pt", 	tau2pt 		, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2Eta", 	tau2eta 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2Phi", 	tau2phi 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  candTable->addColumn<float>("tau2DM", 	tau2decaymode 	, 	"", nanoaod::FlatTable::FloatColumn, 10);

  auto singleTable = std::make_unique<nanoaod::FlatTable>(1, SVFitName_+"MET", false); 
  singleTable->addColumn<float>("px", 		METxVec 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("py", 		METyVec 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("cov00",	covMET00 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("cov01", 	covMET01 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("cov10", 	covMET10 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("cov11", 	covMET11 	, 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<float>("significance", significanceVec , 	"", nanoaod::FlatTable::FloatColumn, 10);
  singleTable->addColumn<int>("isInteresting",  isInteresting 	, 	"", nanoaod::FlatTable::IntColumn);
  singleTable->addColumn<int>("isValid", 	isValid 	, 	"", nanoaod::FlatTable::IntColumn);
 
  iEvent.put(std::move(singleTable),SVFitName_+"MET");
  iEvent.put(std::move(candTable),SVFitName_);
  iEvent.put(std::move(selCandPf)); 
}


classic_svFit::MeasuredTauLepton::kDecayType ClassicSVfitInterface::GetDecayTypeFlag (int pdgId)
{
    if (abs(pdgId) == 11) return classic_svFit::MeasuredTauLepton::kTauToElecDecay;
    if (abs(pdgId) == 13) return classic_svFit::MeasuredTauLepton::kTauToMuDecay;
    if (abs(pdgId) == 15) return classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    
    edm::LogWarning("WrongDecayModePdgID")
       << "(ClassicSVfitInterface): unable to identify decay type from pdgId"
       << "     ---> Decay will be treated as an hadronic decay";
    return classic_svFit::MeasuredTauLepton::kTauToHadDecay;
}

// decide if leptons 1 and 2 must be switched to respect SVfit conventions
bool ClassicSVfitInterface::Switch (classic_svFit::MeasuredTauLepton::kDecayType type1, double pt1, float l1_iso, classic_svFit::MeasuredTauLepton::kDecayType type2, double pt2, float l2_iso)
{
    // e e, mu mu, tau tau
    if (type1 == type2)
    {
      if (type1 == classic_svFit::MeasuredTauLepton::kTauToHadDecay && type2 == classic_svFit::MeasuredTauLepton::kTauToHadDecay)
        return (l1_iso < l2_iso);
      else
        return (pt1 < pt2);
    }
    
    // e tau, mu tau
    if ( (type1 == classic_svFit::MeasuredTauLepton::kTauToElecDecay || type1 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) &&
         type2 == classic_svFit::MeasuredTauLepton::kTauToHadDecay ) {return false;}
    if ( (type2 == classic_svFit::MeasuredTauLepton::kTauToElecDecay || type2 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) &&
         type1 == classic_svFit::MeasuredTauLepton::kTauToHadDecay ) {return true;}

    // e mu
    if (type1 == classic_svFit::MeasuredTauLepton::kTauToElecDecay && type2 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) {return false;}
    if (type2 == classic_svFit::MeasuredTauLepton::kTauToElecDecay && type1 == classic_svFit::MeasuredTauLepton::kTauToMuDecay) {return true;}
    
    cout << "SVfit Standalone: ordering not done (should never happen)" << endl;
    return false;
}

// set mass (pdg ele/mu or leave cand mass for tauh
double ClassicSVfitInterface::GetMass (classic_svFit::MeasuredTauLepton::kDecayType type, double candMass)
{
    if (type == classic_svFit::MeasuredTauLepton::kTauToElecDecay) return 0.51100e-3;
    if (type == classic_svFit::MeasuredTauLepton::kTauToMuDecay)   return 105.658e-3;

    return candMass; // for tauh and all exceptions return cand mass
}

bool ClassicSVfitInterface::IsInteresting (const reco::Candidate *l1, const reco::Candidate *l2)
{
  int apdg1 = abs(l1->pdgId());
  int apdg2 = abs(l2->pdgId());

  int nmu = 0;
  int nele = 0;
  int ntau = 0;

  if (apdg1 == 13) nmu++;
  if (apdg1 == 11) nele++;
  if (apdg1 == 15) ntau++;

  if (apdg2 == 13) nmu++;
  if (apdg2 == 11) nele++;
  if (apdg2 == 15) ntau++;

  pairType pType = kOther;
  if (nmu == 1 && nele == 0 && ntau == 1) pType = kMuHad;
  if (nmu == 0 && nele == 1 && ntau == 1) pType = kEHad;
  if (nmu == 0 && nele == 0 && ntau == 2) pType = kHadHad;
  if (nmu == 2 && nele == 0 && ntau == 0) pType = kMuMu;
  if (nmu == 0 && nele == 2 && ntau == 0) pType = kEE;
  if (nmu == 1 && nele == 1 && ntau == 0) pType = kEMu;

  ///////

  // switch to apply different requirements to the objects
  //if (deltaR(l1->p4(), l2->p4()) < 0.1)
  if (deltaR(l1->p4(), l2->p4()) < 0.4)
    return false; // for overlap removal

  ///////

  // create pointers with usual pair ordering -- easier to apply cuts later
  const reco::Candidate* dau1;
  const reco::Candidate* dau2;

  if (pType == kMuHad)
  {
    dau1 = (apdg1 == 13 ? l1 : l2);
    dau2 = (apdg1 == 13 ? l2 : l1);

    //if (dau1->pt() < 17.)
    if (dau1->pt() < 20.)
      return false;

    if (dau2->pt() < 20.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1) // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.2); //Commented during March 2020 sync: inconsistency with KLUB
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    //if (!iso1 || !iso2)
    if (!iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kEHad)
  {
    dau1 = (apdg1 == 11 ? l1 : l2);
    dau2 = (apdg1 == 11 ? l2 : l1);

    //if (dau1->pt() < 19.)
    if (dau1->pt() < 20.)
      return false;

    if (dau2->pt() < 20.)
      return false;

    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.3);
    //bool iso1 = (userdatahelpers::getUserFloat(l1,"combRelIsoPF") < 0.2); //Commented during March 2020 sync: inconsistency with KLUB
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    //if (!iso1 || !iso2)
    if (!iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kHadHad)
  {
    dau1 = ((l1->pt() > l2->pt()) ? l1 : l2);
    dau2 = ((l1->pt() > l2->pt()) ? l2 : l1);

    //if (dau1->pt() < 30.)
    if (dau1->pt() < 20.)
      return false;
    
    //if (dau2->pt() < 30.)
    if (dau2->pt() < 20.)
      return false;
    
    if (userdatahelpers::getUserInt(l1,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;
    
    if (userdatahelpers::getUserInt(l2,"decayModeFindingNewDMs") != 1)  // decayModeFinding == decayModeFindingOldDMs
      return false;

    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVLooseIsolationMVArun2v1DBoldDMwLT") == 1);
    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") == 1); //FRA 2017
    //bool iso1 = (userdatahelpers::getUserInt(l1,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    //bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso1 = (userdatahelpers::getUserInt(l1,"byVVVLooseDeepTau2017v2p1VSjet") == 1);
    bool iso2 = (userdatahelpers::getUserInt(l2,"byVVVLooseDeepTau2017v2p1VSjet") == 1);

    if (!iso1 || !iso2)
      return false;

    return true; // passed all requirements
  }

  else if (pType == kMuMu)
    return false;
  
  else if (pType == kEE)
    return false;
  
  else if (pType == kEMu)
    return false;
  
  else
  {
    // should never happen
    edm::LogWarning("Unrecognised pair") << "(ClassicSVfitInterface) Warning! could not assess the pair type, won't compute SVFit";
    return false;
  }
}

//------------ method called once each stream before processing any runs, lumis or events  ------------
void ClassicSVfitInterface::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void ClassicSVfitInterface::endStream() 
{
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ClassicSVfitInterface);
