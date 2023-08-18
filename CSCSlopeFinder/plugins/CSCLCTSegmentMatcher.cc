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


struct SegmentLCTData
{
  void init();
  TTree* book(TTree *t);
  //============ Muon Info ================//
  int muon_charge; float muon_pt; float muon_eta;
  unsigned long long evtNum; unsigned long long lumiBlock; int runNum;

};

void SegmentLCTData::init()
{
  //=========== Muon Info ===============//
  float value = 99999;
  muon_charge = value; muon_pt = value; muon_eta = value;
  evtNum = value; lumiBlock = value; runNum = value;

}

TTree* SegmentLCTData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("GEM_CSC_Trigger", "GEM_CSC_Trigger");

  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt); t->Branch("muon_eta", &muon_eta);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("runNum", &runNum);

  return t;
}

class CSCLCTSegmentMatcher : public edm::one::EDAnalyzer<> {
public:
  explicit CSCLCTSegmentMatcher(const edm::ParameterSet&);
  ~CSCLCTSegmentMatcher(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;



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

  SegmentLCTData data_;
  TTree* tree;

  bool isMC;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};


CSCLCTSegmentMatcher::CSCLCTSegmentMatcher(const edm::ParameterSet& iConfig)
  : gemGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin CSC_LCT_Segment_Matcher" << endl;
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
CSCLCTSegmentMatcher::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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

  cout << "New Event" << endl;


  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    if (not mu->isGlobalMuon()) continue;

    //Lets find the ME1/1 segments from this muon
    if(!(mu->isStandAloneMuon())){continue;}
    cout << "New Muon" << endl;
    const reco::Track* Track = mu->outerTrack().get();
    if (Track->validFraction() > 0.0) continue;
    for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
      const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId();
      uint16_t RecHitDetId = RecHitId.det();
      if (RecHitDetId != DetId::Muon) continue;
      uint16_t RecHitSubDet = RecHitId.subdetId();
      if (RecHitSubDet != (uint16_t)MuonSubdetId::CSC) continue;
      CSCDetId SegmentDetId = CSCDetId(RecHitId);
      if (not(SegmentDetId.station() == 1 and SegmentDetId.ring() == 1 and RecHit->dimension() == 4)) continue;
      RecSegment* Rec_segment = (RecSegment*)RecHit;
      const CSCSegment* ME11_segment = (CSCSegment*)Rec_segment;
      //We have a CSCDet with a segment, lets now try to find a matching LCT
      for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j = correlatedlcts->begin(); j != correlatedlcts->end(); j++){
        CSCDetId LCTDetId = (*j).first;
        if (not(LCTDetId.station() == 1 and LCTDetId.ring() == 1)) continue;
        if (not(SegmentDetId == LCTDetId)) continue;
        cout << "Found a Match!!!" << endl;
        cout << "Segment DetID: " << SegmentDetId << endl;
        cout << "LCT DetID    : " << LCTDetId << endl;
        cout << "Lets look at some higher level info" << endl;

        std::vector<CSCCorrelatedLCTDigi>::const_iterator digiItr = (*j).second.first;
        std::vector<CSCCorrelatedLCTDigi>::const_iterator last = (*j).second.second;
        for (; digiItr != last; ++digiItr){
          auto SegmentDetIdL4 = CSCDetId(SegmentDetId.endcap(), SegmentDetId.station(), SegmentDetId.ring(), SegmentDetId.chamber(), 4);
          const CSCLayer* ME11_layer = CSCGeometry_->layer(SegmentDetIdL4);
          const CSCLayerGeometry* ME11_layer_geo = ME11_layer->geometry();
          int SegmentStrip = ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
          cout << "SegStrip:LCTStrip -- " << SegmentStrip << ":" << digiItr->getStrip() << endl;
        }
      }
    }
    //Now find the LCTs from this muon and check if the chambers match
  }

  /*
  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j = correlatedlcts->begin(); j != correlatedlcts->end(); j++){
    CSCDetId SegmentDetId = (*j).first;
    if (not(SegmentDetId.station() == 1 and SegmentDetId.ring() == 1)) continue;
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
  */
}




void CSCLCTSegmentMatcher::beginJob(){}
void CSCLCTSegmentMatcher::endJob(){}

DEFINE_FWK_MODULE(CSCLCTSegmentMatcher);
