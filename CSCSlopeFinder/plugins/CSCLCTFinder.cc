#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//user include files below
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;


struct LCTData
{
  void init();
  TTree* book(TTree *t);
  //============ Muon Info ================//
  int muon_charge; float muon_pt; float muon_eta;
  unsigned long long evtNum; unsigned long long lumiBlock; int runNum;

};

void LCTData::init()
{
  //=========== Muon Info ===============//
  float value = 99999;
  muon_charge = value; muon_pt = value; muon_eta = value;
  evtNum = value; lumiBlock = value; runNum = value;

}

TTree* LCTData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("GEM_CSC_Trigger", "GEM_CSC_Trigger");

  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt); t->Branch("muon_eta", &muon_eta);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("runNum", &runNum);

  return t;
}

class CSCLCTFinder : public edm::one::EDAnalyzer<> {
public:
  explicit CSCLCTFinder(const edm::ParameterSet&);
  ~CSCLCTFinder(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;



  //edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> cscCorrLCTs_;
  //edm::Handle<CSCCorrelatedLCTDigiCollection> cscCorrLCTs;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> co_token;

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;

  edm::Service<TFileService> fs;
  MuonServiceProxy* theService_;
  
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  
  edm::ESHandle<GEMGeometry> GEMGeometry_;
  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool debug;

  LCTData data_;
  TTree* tree;

  bool isMC;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};


CSCLCTFinder::CSCLCTFinder(const edm::ParameterSet& iConfig)
  : gemGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin CSCLCTFinder" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
  co_token = consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter<edm::InputTag>("corrlctDigiTag"));
  //cscCorrLCTs_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag("cscCorrLCTs"));


  debug = iConfig.getParameter<bool>("debug");
  std::cout << "debug " << debug << std::endl;

  tree = data_.book(tree);

}


void
CSCLCTFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GEMGeometry_ = &iSetup.getData(gemGeomToken_);
  CSCGeometry_ = &iSetup.getData(cscGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  iEvent.getByToken(gemRecHits_, gemRecHits);
  edm::Handle<View<reco::Muon> > muons;
  //iEvent.getByToken(cscCorrLCTs_, cscCorrLCTs);
  edm::Handle<CSCCorrelatedLCTDigiCollection> correlatedlcts;
  iEvent.getByToken(co_token, correlatedlcts);
  //cout << "Size of muonDigis? " << muonDigis->size() << endl;

  if (debug) cout << "New! EventNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (debug) cout << "There are " << muons->size() << " muons in event" << endl;
  if (muons->size() == 0) return;

  edm::Handle<CSCSegmentCollection> cscSegments;
  if (! iEvent.getByToken(cscSegments_, cscSegments)){std::cout << "Bad segments" << std::endl;}

  /*
  //Lets loop over all ME1/1 chambers? :D
  for (const auto& ch : CSCGeometry_->layers()) {
    if (!(ch->id().station() == 1 and (ch->id().ring() == 1 or ch->id().ring() == 4))) continue;
    cout << "At new ME1/1! Lets print out the info" << endl;
    auto chID = ch->id();
    cout << chID.endcap() << chID.station() << chID.ring() << chID.chamber() << chID.layer() << endl;
  }
  */

  /*
  for(CSCCorrelatedLCTDigiCollection::const_iterator cscCorrLCT = cscCorrLCTs->begin(); cscCorrLCT != cscCorrLCTs->end(); cscCorrLCT++){
    auto cscDetID = (cscCorrLCT)->getCSCID();
    cout << "At new LCT " << endl;
    cout << cscDetID.endcap() << cscDetID.station() << cscDetID.ring() << cscDetID.chamber() << cscDetID.layer() << endl;
  }
  */
  //cout << "Starting the corr LCT search" << endl;
  cout << "New Event" << endl;
  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j = correlatedlcts->begin(); j != correlatedlcts->end(); j++){
    //cout << "Looping corr lcts" << endl;
    cout << "New Chamber " << (*j).first << endl;
    std::vector<CSCCorrelatedLCTDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCCorrelatedLCTDigi>::const_iterator last = (*j).second.second;
    for (; digiItr != last; ++digiItr) {
      cout << "Found a LCT" << endl;
      cout << "Strip:        " << digiItr->getStrip() << endl;
      cout << "Slope:        " << digiItr->getSlope() << endl;
      cout << "Bend          " << digiItr->getBend() << endl;
      cout << "Run2 Pattern: " << digiItr->getPattern() << endl;
      cout << "Run3 Pattern: " << digiItr->getRun3Pattern() << endl;
    }
  }

  /*
  for (auto cscCorrLCT = cscCorrLCTs->begin(); cscCorrLCT != cscCorrLCTs->end(); cscCorrLCT++) {
    cout << "New 'corrLCT'" << endl;
    auto range = cscCorrLCTs->get((*cscCorrLCT).first);
      
    for (auto cluster = range.first; cluster != range.second; cluster++) {
      cout << "New 'cluster'" << endl;
      cout << (cluster->isValid()) << endl;
      cout << cluster->getStrip() << endl;
      cout << cluster->getSlope() << endl;
      cout << cluster->getRun3Pattern() << endl;
      cout << cluster->getCSCID() << endl;
    }
  }
  */

  /*
  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    if (debug) cout << "new muon, i = " << i << " is it global? " << mu->isGlobalMuon() << endl;

    if (not mu->isGlobalMuon()) continue;
    
    data_.init();
    data_.muon_pt = mu->pt();
    tree->Fill();
  }
  */
}




void CSCLCTFinder::beginJob(){}
void CSCLCTFinder::endJob(){}

DEFINE_FWK_MODULE(CSCLCTFinder);
