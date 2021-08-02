// -*- C++ -*-
//
// Package:    TestRwt/TestRwt
// Class:      TestRwt
//
/**\class TestRwt TestRwt.cc TestRwt/TestRwt/plugins/TestRwt.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  stqian
//         Created:  Tue, 20 Jul 2021 09:10:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>
#include "TH1D.h"
#include "TFile.h"

using std::cout;
using std::endl;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TestRwt : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TestRwt(const edm::ParameterSet&);
  ~TestRwt();

  // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GenEventInfoProduct> GenToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genSrc_;


  TH1D* h_leppt_cpp;
  TH1D* h_leppt_cpp_raw;
  TH1D* h_neupt_cpp;
  TH1D* h_neupt_cpp_raw;
  TH1D* h_wmass_cpp;
  TH1D* h_wmass_cpp_raw;
};


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
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TestRwt::TestRwt(const edm::ParameterSet& iConfig){
  genSrc_ = consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>( "genSrc") ) ;
  GenToken_ = consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "generator") ) ;

  h_leppt_cpp = nullptr;
  h_leppt_cpp_raw = nullptr;
  h_neupt_cpp = nullptr;
  h_neupt_cpp_raw = nullptr;
  h_wmass_cpp = nullptr;
  h_wmass_cpp_raw = nullptr;

}
//     : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))) {
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   setupDataToken_ = esConsumes<SetupData, SetupRecord>();
// #endif
  //now do what ever initialization is needed
// }

TestRwt::~TestRwt() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void TestRwt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;

//   for (const auto& track : iEvent.get(tracksToken_)) {
//     // do something with track parameters, e.g, plot the charge.
//     // int charge = track.charge();
//   }

// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   // if the SetupData is always needed
//   auto setup = iSetup.getData(setupToken_);
//   // if need the ESHandle to check if the SetupData was there or not
//   auto pSetup = iSetup.getHandle(setupToken_);
// #endif
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

    // cout<<genParticles<<endl;



  //  iEvent.getByLabel("genParticles", genParticles);
   for(size_t i = 0; i < genParticles->size(); ++ i) {
    
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
     if ((st == 1) && check_neu(abs(id))){
        if (pt_neu < p.pt()){
            pt_neu = p.pt();
            idx_neu = i;
            found_neu = true;
        }
     }

     
         
     }
   }
     auto ptLep = pt_lep;
     auto ptNeu = pt_neu;

     double wMass = 0.;
     
     if (found_lep && found_neu){
          auto p_lep = (*genParticles)[idx_lep];
          auto p_neu = (*genParticles)[idx_neu];
         
        
         wMass = inv_mass(p_lep.px(),p_lep.py(),p_lep.pz(),p_lep.energy(),p_neu.px(),p_neu.py(),p_neu.pz(),p_neu.energy());
     }
     else{
         wMass = 0;
     }
     const std::vector<double>& evtWeights = genInfo->weights();
     double theWeight = genInfo->weight();
     double genWeight = theWeight / abs(theWeight);
    //  for(size_t j = 0; j < n; ++ j) {
    //    const Candidate * d = p.daughter( j );
    //    int dauId = d->pdgId();
    //    // . . . 
    //  }
     // . . . 

    //  cout<<"This is info from event\n"<<ptLep<<"\t"<<ptNeu<<"\t"<<wMass<<"\t"<<genWeight<<endl;

     h_leppt_cpp->Fill(ptLep,genWeight);
     h_leppt_cpp_raw->Fill(ptLep);
     h_neupt_cpp->Fill(ptNeu,genWeight);
     h_neupt_cpp_raw->Fill(ptNeu);
     h_wmass_cpp->Fill(wMass,genWeight);
     h_wmass_cpp_raw->Fill(wMass);

    //  h_leppt_cpp->Print();
     
   
}

// ------------ method called once each job just before starting event loop  ------------
void TestRwt::beginJob() {
  // please remove this method if not needed


  cout<<"begin job\n";
  edm::Service<TFileService> fs;

  h_leppt_cpp = fs->make<TH1D>("h_leppt_cpp","h_leppt_cpp",50,0,50);
  h_leppt_cpp_raw = fs->make<TH1D>("h_leppt_cpp_raw","h_leppt_cpp_raw",50,0,50);
  h_neupt_cpp = fs->make<TH1D>("h_neupt_cpp","h_neupt_cpp",50,0,50);
  h_neupt_cpp_raw = fs->make<TH1D>("h_neupt_cpp_raw","h_neupt_cpp_raw",50,0,50);
  h_wmass_cpp = fs->make<TH1D>("h_wmass_cpp","h_wmass_cpp",50,50,100);
  h_wmass_cpp_raw = fs->make<TH1D>("h_wmass_cpp_raw","h_wmass_cpp_raw",50,50,100);
}

// ------------ method called once each job just after ending the event loop  ------------
void TestRwt::endJob() {
  // please remove this method if not needed
  h_leppt_cpp->Print();
  // h_leppt_cpp->Write();
  // h_leppt_cpp_raw->Write();
  // h_neupt_cpp->Write();
  // h_neupt_cpp_raw->Write();
  // h_wmass_cpp->Write();
  // h_wmass_cpp_raw->Write();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void TestRwt::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);

//   //Specify that only 'tracks' is allowed
//   //To use, remove the default given above and uncomment below
//   //ParameterSetDescription desc;
//   //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
//   //descriptions.addWithDefaultLabel(desc);
// }

//define this as a plug-in
DEFINE_FWK_MODULE(TestRwt);
