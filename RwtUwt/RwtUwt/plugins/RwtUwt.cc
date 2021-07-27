// -*- C++ -*-
//
// Package:    RwtUwt/RwtUwt
// Class:      RwtUwt
//
/**\class RwtUwt RwtUwt.cc RwtUwt/RwtUwt/plugins/RwtUwt.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  stqian
//         Created:  Wed, 21 Jul 2021 07:23:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include <string>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"

using std::unique_ptr;
using std::auto_ptr;

//
// class declaration
//

class RwtUwt : public edm::stream::EDProducer<> {
public:
  explicit RwtUwt(const edm::ParameterSet&);
  ~RwtUwt();

  // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  TH1D* hist4Rwt;
  TH1D* hist4RwtRaw;

  TH1D* hist4RwtRatio;

  TFile* rootFile;
  TDirectory* rootDir;

  std::string fileName;
  std::string dirName;
  std::string name4Rwt;
  std::string name4Raw;

  edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;

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

struct ifKeepAndProb{
  int ifKeep;
  double prob;
};

RwtUwt::RwtUwt(const edm::ParameterSet& iConfig)
{
  //register your products
/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
  //now do what ever other initialization is needed
  hist4Rwt = nullptr;
  hist4RwtRaw = nullptr;
  hist4RwtRatio = nullptr;
  rootFile = nullptr;

  fileName = "/home/pku/stqian/NegWgt/CMSSW_11_3_0_pre2/src/TestRwt/TestRwt/python/TestRwt.root";
  dirName = "TestRwt";
  name4Rwt = "h_leppt_cpp";
  name4Raw = "h_leppt_cpp_raw";

  genSrc_ = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
  GenToken_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;


  // produces<std::vector<int>>( "ifkeepAndProb" ).setBranchAlias( "ifkeepAndProb");
  produces<std::vector<int>>( "ifkeep" ).setBranchAlias( "ifkeep");
  // produces<std::vector<double>>("Probs").setBranchAlias("Probs");
}

RwtUwt::~RwtUwt() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}


bool check_lep(int id){
    if (id !=  13 && id != 11 && id != 15){
        return false;
    }
    else {
        return true;
    }
}

bool check_neu(int id){
    if (id !=  14 && id != 12 && id != 16){
        return false;
    }
    else {
        return true;
    }
}

double inv_mass(double px1, double py1, double pz1, double e1, double px2, double py2, double pz2, double e2){
    double e = e1 + e2;
    double px = px1 + px2;
    double py = py1 + py2;
    double pz = pz1 + pz2;
    
    double ret = e*e - px*px - py*py - pz*pz;
    if (ret < 0){
        return 0.;
    }
    else{
        return sqrt(ret);
    }
}



//
// member functions
//

// ------------ method called to produce the data  ------------
void RwtUwt::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;

  edm::Handle<GenEventInfoProduct> genInfo;
        //iEvent.getByLabel( "generator", genEvtInfo );
        iEvent.getByToken(GenToken_,genInfo);
  edm::Handle<edm::View<reco::GenParticle> > genParticles;//define genParticle
    //iEvent.getByLabel(InputTag("prunedGenParticles"), genParticles);
    iEvent.getByToken(genSrc_, genParticles);



    int idx_lep = -1;
    double pt_lep = -1;
    bool found_lep = false;
    int idx_neu = -1; 
    double pt_neu = -1;
    bool found_neu = false;



  //  iEvent.getByLabel("genParticles", genParticles);
   for(size_t i = 0; i < genParticles->size(); ++ i) {
    //  int idx_lep = -1;
    //     double pt_lep = -1;
    //     bool found_lep = false;
    //     // int idx_neu = -1; 
    //     // double pt_neu = -1;
    //     // bool found_neu = false;

        for(size_t i = 0; i < genParticles->size(); ++ i) {
     auto p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();
     if ((st == 1) && check_lep(abs(id))){
        if (pt_lep < p.pt()){
            pt_lep = p.pt();
            idx_lep = i;
            found_lep = true;
        }
     }
    //  if ((st == 1) && check_neu(abs(id))){
    //     if (pt_neu < p.pt()){
    //         pt_neu = p.pt();
    //         idx_neu = i;
    //         found_neu = true;
    //     }
    //  }

     
         
     }
   }
     auto ptLep = pt_lep;
    //  auto ptNeu = pt_neu;

    //  double wMass = 0;
     
    //  if (found_lep && found_neu){
    //       auto p_lep = (*genParticles)[idx_lep];
    //       auto p_neu = (*genParticles)[idx_neu];
         
        
    //      wMass = inv_mass(p_lep.px(),p_lep.py(),p_lep.pz(),p_lep.energy(),p_neu.px(),p_neu.py(),p_neu.pz(),p_neu.energy());
    //  }
    //  else{
    //      wMass = 0;
    //  }
     const std::vector<double>& evtWeights = genInfo->weights();
     double theWeight = genInfo->weight();
     double genWeight = theWeight / abs(theWeight);

     auto idx = hist4RwtRatio->FindBin(ptLep);
     double prob = hist4RwtRatio->GetBinContent(idx);
     
     auto rndm = new TRandom();
     int ifKeep = 0;

     auto prob_rndm = rndm->Rndm();
     
     if(rndm->Rndm() < prob){
       ifKeep = 1;
     }

     std::cout<<ifKeep<<"\t"<<prob<<std::endl;

    //  auto test = new ifKeepAndProb;
    //  test->prob = prob;
    //  test->ifKeep = ifKeep;

     unique_ptr<std::vector<int>> ifkeep(new std::vector<int>);
     ifkeep->reserve(1);
     ifkeep->push_back(ifKeep);
     iEvent.put(std::move(ifkeep),"ifkeep");

    //  for(size_t j = 0; j < n; ++ j) {
    //    const Candidate * d = p.daughter( j );
    //    int dauId = d->pdgId();
    //    // . . . 
    //  }
     // . . . 
/* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
*/

/* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/

  


}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void RwtUwt::beginStream(edm::StreamID) {
   rootFile = new TFile(fileName.c_str(),"read");
   rootDir = (TDirectory*)rootFile->Get(dirName.c_str());
   rootDir->cd();
  hist4Rwt = (TH1D*)rootDir->Get(name4Rwt.c_str());
  hist4RwtRaw = (TH1D*)rootDir->Get(name4Raw.c_str());

  hist4RwtRatio = (TH1D*)hist4Rwt->Clone();
  hist4RwtRatio->Divide(hist4RwtRaw);
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void RwtUwt::endStream() {
  // please remove this method if not needed
  rootFile->Close();
}

// ------------ method called when starting to processes a run  ------------

void
RwtUwt::beginRun(edm::Run const&, edm::EventSetup const&)
{
 
}


// ------------ method called when ending the processing of a run  ------------

void
RwtUwt::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void
RwtUwt::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
RwtUwt::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void RwtUwt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

//define this as a plug-in
DEFINE_FWK_MODULE(RwtUwt);
