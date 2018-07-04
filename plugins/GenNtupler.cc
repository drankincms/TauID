// -*- C++ -*-
//
// Package:    MonteCarloInfo/GenNtupler
// Class:      GenNtupler
// 
/**\class GenNtupler GenNtupler.cc UserCode/GenNtupler/plugins/GenNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  dylan rankin
//         Created:  Fri, 25 May 2018 16:29:00 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
//
// class declaration
//

class GenNtupler : public edm::EDAnalyzer {
   public:
      explicit GenNtupler(const edm::ParameterSet&);
      ~GenNtupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void addMother(const reco::GenParticle& iPart,int &iPartId,int &ig,  std::vector<edm::Ptr<reco::GenParticle>> &iMothers,  std::vector<std::pair<edm::Ptr<reco::GenParticle>,int> >  &iMothersAdded);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree* gentree;

      std::vector<float> genpt_;
      std::vector<float> geneta_;
      std::vector<float> genphi_;
      std::vector<float> genet_;
      std::vector<int>   genid_;
      std::vector<int>   genpartid_;
      std::vector<int>   genstatus_;
      std::vector<int>   genparent_;
      std::vector<int>   genindex_;

      std::vector<float> genjetpt_;
      std::vector<float> genjeteta_;
      std::vector<float> genjetphi_;
      std::vector<float> genjetet_;
      std::vector<int>   genjetid_;

      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genJetToken_;
      double minGenPt_;
      double minJetPt_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenNtupler::GenNtupler(const edm::ParameterSet& iConfig):
	genToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("generator"))),
	genJetToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"))),
        minGenPt_(iConfig.getParameter<double>("minGenParticlePt")),
        minJetPt_(iConfig.getParameter<double>("minGenJetPt"))


{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  //weightHist = fs->make<TH1F>("weightHist","MC Weight",10000,-10,10);
  gentree = fs->make<TTree>("gentree","gentree");

  gentree->Branch("genpt", &genpt_);
  gentree->Branch("geneta", &geneta_);
  gentree->Branch("genphi", &genphi_);
  gentree->Branch("genet", &genet_);
  gentree->Branch("genid", &genid_);
  gentree->Branch("genpartid", &genpartid_);
  gentree->Branch("genindex", &genindex_);
  gentree->Branch("genparent", &genparent_);
  gentree->Branch("genstatus", &genstatus_);

  gentree->Branch("genjetpt", &genjetpt_);
  gentree->Branch("genjeteta", &genjeteta_);
  gentree->Branch("genjetphi", &genjetphi_);
  gentree->Branch("genjetet", &genjetet_);
  gentree->Branch("genjetid", &genjetid_);
}


GenNtupler::~GenNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //load genparticles collection
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genToken_, genParticles);

   //load genjets collection
   edm::Handle<std::vector<reco::GenJet>> genJets;
   iEvent.getByToken(genJetToken_, genJets);

   genpt_.clear();
   geneta_.clear();
   genphi_.clear();
   genet_.clear();
   genid_.clear();
   genpartid_.clear();
   genindex_.clear();
   genstatus_.clear();
   genparent_.clear();
   genjetpt_.clear();
   genjeteta_.clear();
   genjetphi_.clear();
   genjetet_.clear();
   genjetid_.clear();

   std::vector<edm::Ptr<reco::GenParticle>> lMothers;
   std::vector<std::pair<edm::Ptr<reco::GenParticle>,int> > lMothersAdded;
   std::vector<std::pair<TLorentzVector,int>> promptgenv;
   std::pair<TLorentzVector,int> tmpv;
   for (reco::GenParticleCollection::const_iterator itGenP = genParticles->begin(); itGenP!=genParticles->end(); ++itGenP) {
     const reco::GenParticle & p = *itGenP;
     tmpv.first.SetPtEtaPhiE(p.pt(),p.eta(),p.phi(),p.energy());
     tmpv.second = p.pdgId();
     promptgenv.push_back(tmpv);
     edm::Ptr<reco::GenParticle> thePtr(genParticles, itGenP - genParticles->begin());
     lMothers.emplace_back(thePtr);
   }

   int ig = 0;
   int lPartId = 0;
   for (const reco::GenJet &j : *genJets){
       if (j.pt()<minJetPt_) continue;
       genjetpt_.push_back(j.pt());
       genjeteta_.push_back(j.eta());
       genjetphi_.push_back(j.phi());
       genjetet_.push_back(j.et());
       int jetid = 0;
       tmpv.first.SetPtEtaPhiE(j.pt(),j.eta(),j.phi(),j.energy());
       for (unsigned int i = 0; i<promptgenv.size(); i++) {
           if (tmpv.first.DeltaR(promptgenv[i].first)<0.3 && fabs((tmpv.first.Pt()/promptgenv[i].first.Pt())-1.25)<0.75) {
               jetid = promptgenv[i].second;
               break;
           }
       }
       genjetid_.push_back(jetid);
       for (unsigned int id = 0, nd = j.numberOfDaughters(); id < nd; ++id) {
	 if ( !(j.daughterPtr(id).isNonnull()) ) continue;
	 if ( !(j.daughterPtr(id).isAvailable()) ) continue;
	 const reco::Candidate &_ijet_const = dynamic_cast<const reco::Candidate &>(*j.daughter(id));
	 if (_ijet_const.pt()<minGenPt_) continue;
	 genpt_.push_back(_ijet_const.pt());
	 geneta_.push_back(_ijet_const.eta());
	 genphi_.push_back(_ijet_const.phi());
	 genet_.push_back(_ijet_const.et());
	 genid_.push_back(_ijet_const.pdgId());
	 genstatus_.push_back(_ijet_const.status());
	 genpartid_.push_back(lPartId);
	 genindex_.push_back(ig);
	 lPartId++;
	 for (reco::GenParticleCollection::const_iterator itGenP = genParticles->begin(); itGenP!=genParticles->end(); ++itGenP) {
	   //if(!itGenP->statusFlags().isPrompt()) continue;
	   //if(itGenP->pdgId()  != _ijet_const.pdgId()) continue;
	   if(fabs(itGenP->pt()-_ijet_const.pt())   > 0.1)  continue;
	   if(fabs(itGenP->eta()-_ijet_const.eta()) > 0.01) continue;
	   addMother(*itGenP,lPartId,ig,lMothers,lMothersAdded);
	   //std::cout << "==> " << _ijet_const.pdgId() << " -- " << itGenP->pdgId() << " -- " << itGenP->eta() << " -- " << _ijet_const.eta() << " -- " << _ijet_const.phi() << " -- " << itGenP->phi() << " -- " << itGenP->status()  << " -- " << _ijet_const.status()  << " -- " << _ijet_const.pt() << " -- " << itGenP->pt() << std::endl;
	   break;
	 }
       }
       ig++;
   }

   gentree->Fill();

}
void GenNtupler::addMother(const reco::GenParticle& iPart,int &iPartId,int &ig,  std::vector<edm::Ptr<reco::GenParticle>> &iMothers,  std::vector<std::pair<edm::Ptr<reco::GenParticle>,int> >  &iMothersAdded) { 
  int pId = -2;
  genparent_.push_back(pId);
  int genIndex =  genparent_.size()-1;
  if(iPart.numberOfMothers() >  0 ) {
    edm::Ptr<reco::GenParticle> lMomPtr = edm::refToPtr(iPart.motherRef()); 
    for(unsigned int im=0; im < iMothers.size(); ++im) { 
      if(iMothers[im] != lMomPtr) continue;
      for(unsigned int im1=0; im1 < iMothersAdded.size(); ++im1) { 
	if(iMothers[im] != iMothersAdded[im1].first) continue;
	pId = iMothersAdded[im1].second;
      }
      if(pId == -2) { //not added
	genpt_.push_back(iMothers[im]->pt());
	geneta_.push_back(iMothers[im]->eta());
	genphi_.push_back(iMothers[im]->phi());
	genet_.push_back(iMothers[im]->et());
	genid_.push_back(iMothers[im]->pdgId());
	genstatus_.push_back(iMothers[im]->status());
	genpartid_.push_back(iPartId);
	genindex_.push_back(ig);
	iMothersAdded.emplace_back(std::pair<edm::Ptr<reco::GenParticle>,int>(lMomPtr,iPartId));
	iPartId++;
	addMother(*lMomPtr,iPartId,ig,iMothers,iMothersAdded);
      }
    }
    genparent_[genIndex] = pId;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GenNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GenNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GenNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GenNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenNtupler);
