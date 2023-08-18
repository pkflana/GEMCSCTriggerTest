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

  //============ Segment Info =============//
  float Segment_GP_x; float Segment_GP_y; float Segment_GP_z;
  float Segment_LP_x; float Segment_LP_y; float Segment_LP_z;
  float Segment_slope_angle;
  int Segment_ME11_endcap; int Segment_ME11_station; 
  int Segment_ME11_ring; int Segment_ME11_chamber;

  //============ Track Info ===============//
  float Track_GP_x; float Track_GP_y; float Track_GP_z;
  float Track_LP_x; float Track_LP_y; float Track_LP_z;
  float Track_slope_angle;
  int Track_ME11_endcap; int Track_ME11_station;
  int Track_ME11_ring; int Track_ME11_chamber;

  //============ LCT Info =================//
  float LCT_slope; int LCT_quality; float LCT_strip;
  int LCT_ME11_endcap; int LCT_ME11_station;
  int LCT_ME11_ring; int LCT_ME11_chamber;

  //============ Segment Prop Info ========//
  float SegmentProp_GP_x; float SegmentProp_GP_y; float SegmentProp_GP_z;
  float SegmentProp_LP_x; float SegmentProp_LP_y; float SegmentProp_LP_z;
  int SegmentProp_GE11_endcap; int SegmentProp_GE11_station;
  int SegmentProp_GE11_ring; int SegmentProp_GE11_chamber;
  int SegmentProp_GE11_layer; int SegmentProp_GE11_eta;

  //============ Track Prop Info ==========//
  float TrackProp_GP_x; float TrackProp_GP_y; float TrackProp_GP_z;
  float TrackProp_LP_x; float TrackProp_LP_y; float TrackProp_LP_z;
  int TrackProp_GE11_endcap; int TrackProp_GE11_station;
  int TrackProp_GE11_ring; int TrackProp_GE11_chamber;
  int TrackProp_GE11_layer; int TrackProp_GE11_eta;

  //============ Segment Hit Info =========//
  float SegmentHit_GP_x; float SegmentHit_GP_y; float SegmentHit_GP_z;
  float SegmentHit_LP_x; float SegmentHit_LP_y; float SegmentHit_LP_z;
  float SegmentHit_residual;
  int SegmentHit_GE11_endcap; int SegmentHit_GE11_station;
  int SegmentHit_GE11_ring; int SegmentHit_GE11_chamber;
  int SegmentHit_GE11_layer; int SegmentHit_GE11_eta;

  //============ Track Hit Info ===========//
  float TrackHit_GP_x; float TrackHit_GP_y; float TrackHit_GP_z;
  float TrackHit_LP_x; float TrackHit_LP_y; float TrackHit_LP_z;
  float TrackHit_residual;
  int TrackHit_GE11_endcap; int TrackHit_GE11_station;
  int TrackHit_GE11_ring; int TrackHit_GE11_chamber;
  int TrackHit_GE11_layer; int TrackHit_GE11_eta;

};

void SegmentLCTData::init()
{
  //=========== Muon Info ===============//
  float value = 99999;
  muon_charge = value; muon_pt = value; muon_eta = value;
  evtNum = value; lumiBlock = value; runNum = value;

  //============ Segment Info =============//
  Segment_GP_x = value; Segment_GP_y = value; Segment_GP_z = value;
  Segment_LP_x = value; Segment_LP_y = value; Segment_LP_z = value;
  Segment_slope_angle = value;
  Segment_ME11_endcap = value; Segment_ME11_station = value;
  Segment_ME11_ring = value; Segment_ME11_chamber = value;

  //============ Track Info ===============//
  Track_GP_x = value; Track_GP_y = value; Track_GP_z = value;
  Track_LP_x = value; Track_LP_y = value; Track_LP_z = value;
  Track_slope_angle = value;
  Track_ME11_endcap = value; Segment_ME11_station = value;
  Track_ME11_ring = value; Segment_ME11_chamber = value;

  //============ LCT Info =================//
  LCT_slope = value; LCT_quality = value; LCT_strip = value;
  LCT_ME11_endcap = value; LCT_ME11_station = value;
  LCT_ME11_ring = value; LCT_ME11_chamber = value;

  //============ Segment Prop Info ========//
  SegmentProp_GP_x = value; SegmentProp_GP_y = value; SegmentProp_GP_z = value;
  SegmentProp_LP_x = value; SegmentProp_LP_y = value; SegmentProp_LP_z = value;
  SegmentProp_GE11_endcap = value; SegmentProp_GE11_station = value;
  SegmentProp_GE11_ring = value; SegmentProp_GE11_chamber = value;
  SegmentProp_GE11_layer = value; SegmentProp_GE11_eta = value;

  //============ Track Prop Info ==========//
  TrackProp_GP_x = value; TrackProp_GP_y = value; TrackProp_GP_z = value;
  TrackProp_LP_x = value; TrackProp_LP_y = value; TrackProp_LP_z = value;
  TrackProp_GE11_endcap = value; TrackProp_GE11_station = value;
  TrackProp_GE11_ring = value; TrackProp_GE11_chamber = value;
  TrackProp_GE11_layer = value; TrackProp_GE11_eta = value;

  //============ Segment Hit Info =========//
  SegmentHit_GP_x = value; SegmentHit_GP_y = value; SegmentHit_GP_z = value;
  SegmentHit_LP_x = value; SegmentHit_LP_y = value; SegmentHit_LP_z = value;
  SegmentHit_residual = value;
  SegmentHit_GE11_endcap = value; SegmentHit_GE11_station = value;
  SegmentHit_GE11_ring = value; SegmentHit_GE11_chamber = value;
  SegmentHit_GE11_layer = value; SegmentHit_GE11_eta = value;

  //============ Track Hit Info ===========//
  TrackHit_GP_x = value; TrackHit_GP_y = value; TrackHit_GP_z = value;
  TrackHit_LP_x = value; TrackHit_LP_y = value; TrackHit_LP_z = value;
  TrackHit_residual = value;
  TrackHit_GE11_endcap = value; TrackHit_GE11_station = value;
  TrackHit_GE11_ring = value; TrackHit_GE11_chamber = value;
  TrackHit_GE11_layer = value; TrackHit_GE11_eta = value;
}

TTree* SegmentLCTData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("GEM_CSC_Trigger", "GEM_CSC_Trigger");

  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt); t->Branch("muon_eta", &muon_eta);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("runNum", &runNum);

  //============ Segment Info =============//
  t->Branch("Segment_GP_x", &Segment_GP_x); t->Branch("Segment_GP_y", &Segment_GP_y); t->Branch("Segment_GP_z", &Segment_GP_z);
  t->Branch("Segment_LP_x", &Segment_LP_x); t->Branch("Segment_LP_y", &Segment_LP_y); t->Branch("Segment_LP_z", &Segment_LP_z);
  t->Branch("Segment_slope_angle", &Segment_slope_angle);
  t->Branch("Segment_ME11_endcap", &Segment_ME11_endcap); t->Branch("Segment_ME11_station", &Segment_ME11_station);
  t->Branch("Segment_ME11_ring", &Segment_ME11_ring); t->Branch("Segment_ME11_chamber", &Segment_ME11_chamber);

  //============ Track Info ===============//
  t->Branch("Track_GP_x", &Track_GP_x); t->Branch("Track_GP_y", &Track_GP_y); t->Branch("Track_GP_z", &Track_GP_z);
  t->Branch("Track_LP_x", &Track_LP_x); t->Branch("Track_LP_y", &Track_LP_y); t->Branch("Track_LP_z", &Track_LP_z);
  t->Branch("Track_slope_angle", &Track_slope_angle);
  t->Branch("Track_ME11_endcap", &Track_ME11_endcap); t->Branch("Track_ME11_station", &Track_ME11_station);
  t->Branch("Track_ME11_ring", &Track_ME11_ring); t->Branch("Track_ME11_chamber", &Track_ME11_chamber);

  //============ LCT Info =================//
  t->Branch("LCT_slope", &LCT_slope); t->Branch("LCT_quality", &LCT_quality); t->Branch("LCT_strip", &LCT_strip);
  t->Branch("LCT_ME11_endcap", &LCT_ME11_endcap); t->Branch("LCT_ME11_station", &LCT_ME11_station);
  t->Branch("LCT_ME11_ring", &LCT_ME11_ring); t->Branch("LCT_ME11_chamber", &LCT_ME11_chamber);

  //============ Segment Prop Info ========//
  t->Branch("SegmentProp_GP_x", &SegmentProp_GP_x); t->Branch("SegmentProp_GP_y", &SegmentProp_GP_y); t->Branch("SegmentProp_GP_z", &SegmentProp_GP_z);
  t->Branch("SegmentProp_LP_x", &SegmentProp_LP_x); t->Branch("SegmentProp_LP_y", &SegmentProp_LP_y); t->Branch("SegmentProp_LP_z", &SegmentProp_LP_z);
  t->Branch("SegmentProp_GE11_endcap", &SegmentProp_GE11_endcap); t->Branch("SegmentProp_GE11_station", &SegmentProp_GE11_station);
  t->Branch("SegmentProp_GE11_ring", &SegmentProp_GE11_ring); t->Branch("SegmentProp_GE11_chamber", &SegmentProp_GE11_chamber);
  t->Branch("SegmentProp_GE11_layer", &SegmentProp_GE11_layer); t->Branch("SegmentProp_GE11_eta", &SegmentProp_GE11_eta);

  //============ Track Prop Info ==========//
  t->Branch("TrackProp_GP_x", &TrackProp_GP_x); t->Branch("TrackProp_GP_y", &TrackProp_GP_y); t->Branch("TrackProp_GP_z", &TrackProp_GP_z);
  t->Branch("TrackProp_LP_x", &TrackProp_LP_x); t->Branch("TrackProp_LP_y", &TrackProp_LP_y); t->Branch("TrackProp_LP_z", &TrackProp_LP_z);
  t->Branch("TrackProp_GE11_endcap", &TrackProp_GE11_endcap); t->Branch("TrackProp_GE11_station", &TrackProp_GE11_station);
  t->Branch("TrackProp_GE11_ring", &TrackProp_GE11_ring); t->Branch("TrackProp_GE11_chamber", &TrackProp_GE11_chamber);
  t->Branch("TrackProp_GE11_layer", &TrackProp_GE11_layer); t->Branch("TrackProp_GE11_eta", &TrackProp_GE11_eta);

  //============ Segment Hit Info =========//
  t->Branch("SegmentHit_GP_x", &SegmentHit_GP_x); t->Branch("SegmentHit_GP_y", &SegmentHit_GP_y); t->Branch("SegmentHit_GP_z", &SegmentHit_GP_z);
  t->Branch("SegmentHit_LP_x", &SegmentHit_LP_x); t->Branch("SegmentHit_LP_y", &SegmentHit_LP_y); t->Branch("SegmentHit_LP_z", &SegmentHit_LP_z);
  t->Branch("SegmentHit_residual", &SegmentHit_residual);
  t->Branch("SegmentHit_GE11_endcap", &SegmentHit_GE11_endcap); t->Branch("SegmentHit_GE11_station", &SegmentHit_GE11_station);
  t->Branch("SegmentHit_GE11_ring", &SegmentHit_GE11_ring); t->Branch("SegmentHit_GE11_chamber", &SegmentHit_GE11_chamber);
  t->Branch("SegmentHit_GE11_layer", &SegmentHit_GE11_layer); t->Branch("SegmentHit_GE11_eta", &SegmentHit_GE11_eta);

  //============ Track Hit Info ===========//
  t->Branch("TrackHit_GP_x", &TrackHit_GP_x); t->Branch("TrackHit_GP_y", &TrackHit_GP_y); t->Branch("TrackHit_GP_z", &TrackHit_GP_z);
  t->Branch("TrackHit_LP_x", &TrackHit_LP_x); t->Branch("TrackHit_LP_y", &TrackHit_LP_y); t->Branch("TrackHit_LP_z", &TrackHit_LP_z);
  t->Branch("TrackHit_residual", &TrackHit_residual);
  t->Branch("TrackHit_GE11_endcap", &TrackHit_GE11_endcap); t->Branch("TrackHit_GE11_station", &TrackHit_GE11_station);
  t->Branch("TrackHit_GE11_ring", &TrackHit_GE11_ring); t->Branch("TrackHit_GE11_chamber", &TrackHit_GE11_chamber);
  t->Branch("TrackHit_GE11_layer", &TrackHit_GE11_layer); t->Branch("TrackHit_GE11_eta", &TrackHit_GE11_eta);

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
  edm::Handle<CSCCorrelatedLCTDigiCollection> correlatedlcts;
  iEvent.getByToken(co_token, correlatedlcts);

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
    if (debug) cout << "New Muon" << endl;
    const reco::Track* Track = mu->outerTrack().get();
    if (Track->validFraction() > 0.0) continue;
    for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){

      //The tree fill will happen at this level (segment level), lets initialize our variables
      data_.init();
      const CSCSegment* ME11_segment;
      TrajectoryStateOnSurface TSOS_Segment;
      TrajectoryStateOnSurface TSOS_Track;
      CSCCorrelatedLCTDigi LCTDigiMatch;
      CSCDetId SegmentCSCDetId;
      CSCDetId TrackCSCDetId;
      CSCDetId LCTDetId;
      TrajectoryStateOnSurface TSOS_Segment_On_GEM;
      TrajectoryStateOnSurface TSOS_Track_On_GEM;
      GEMRecHit GEM_Hit_Match_Segment;
      GEMRecHit GEM_Hit_Match_Track;
      GEMDetId SegmentGEMDetId;
      GEMDetId TrackGEMDetId;
      

      const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId();
      uint16_t RecHitDetId = RecHitId.det();
      if (RecHitDetId != DetId::Muon) continue;
      uint16_t RecHitSubDet = RecHitId.subdetId();
      if (RecHitSubDet != (uint16_t)MuonSubdetId::CSC) continue;
      SegmentCSCDetId = CSCDetId(RecHitId);
      TrackCSCDetId = SegmentCSCDetId;
      if (not(SegmentCSCDetId.station() == 1 and SegmentCSCDetId.ring() == 1 and RecHit->dimension() == 4)) continue;
      RecSegment* Rec_segment = (RecSegment*)RecHit;
      ME11_segment = (CSCSegment*)Rec_segment;

      auto SegmentCSCDetIdL4 = CSCDetId(SegmentCSCDetId.endcap(), SegmentCSCDetId.station(), SegmentCSCDetId.ring(), SegmentCSCDetId.chamber(), 4);
      const CSCLayer* ME11_layer = CSCGeometry_->layer(SegmentCSCDetIdL4);
      const CSCLayerGeometry* ME11_layer_geo = ME11_layer->geometry();

      //We have a segment, lets now try to find the track at this same position
      reco::TransientTrack TTrack = ttrackBuilder_->build(Track);
      GlobalPoint ME11_segment_GP = ME11_layer->toGlobal(ME11_segment->localPosition());
      TrajectoryStateOnSurface track_at_segment = TTrack.stateOnSurface(ME11_segment_GP);
      if (!(track_at_segment.isValid())) continue;
      LocalPoint track_at_segment_LP = ME11_layer->toLocal(track_at_segment.globalPosition());
      if (debug) cout << "Debug track strip issue? TrackLP:SegmentLP " << track_at_segment_LP << ":" << ME11_segment->localPosition() << endl;
      if (debug) cout << "Debug track strip issue? TrackGP:SegmentGP " << track_at_segment.globalPosition() << ":" << ME11_segment_GP << endl;

      //We have a CSCDet with a segment, lets now try to find a matching LCT
      for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j = correlatedlcts->begin(); j != correlatedlcts->end(); j++){
        LCTDetId = (*j).first;
        if (not(LCTDetId.station() == 1 and LCTDetId.ring() == 1)) continue;
        if (not(SegmentCSCDetId == LCTDetId)) continue;
        if (debug){
          cout << "Found a Match!!!" << endl;
          cout << "Segment DetID: " << SegmentCSCDetId << endl;
          cout << "LCT DetID    : " << LCTDetId << endl;
          cout << "Lets look at some higher level info" << endl;
        }

        std::vector<CSCCorrelatedLCTDigi>::const_iterator digiItr = (*j).second.first;
        std::vector<CSCCorrelatedLCTDigi>::const_iterator last = (*j).second.second;
        //Pick the best LCT Digi to match the segment
        cout << "nope ):" << endl;
        for (; digiItr != last; ++digiItr){
          float SegmentStrip = ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
          float SegmentFracStrip = ME11_layer_geo->strip(ME11_segment->localPosition());
          const float TrackStrip = ME11_layer_geo->nearestStrip(track_at_segment_LP);
          const float TrackFracStrip = ME11_layer_geo->strip(track_at_segment_LP);
          if (debug)  cout << "SegStrip:TrackStrip:LCTStrip Frac -- " << SegmentFracStrip << ":" << TrackFracStrip << ":" << digiItr->getFractionalStrip() << endl;

          //If new Itr is closer to the segment strip than the old one, replace
          if (abs(SegmentFracStrip - digiItr->getFractionalStrip()) < abs(SegmentFracStrip - LCTDigiMatch.getFractionalStrip())) LCTDigiMatch = *digiItr;
        }
        //Found the best DigiItr, now lets work on propagating to GEM
        //Set Up Segment Propagation
        DetId segDetId = ME11_segment->geographicalId();
        const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId());
        LocalVector momentum_at_surface = (Track->outerP()) * (ME11_segment->localDirection());
        LocalTrajectoryParameters param(ME11_segment->localPosition(), momentum_at_surface, mu->charge());
        AlgebraicSymMatrix mat(5,0);
        mat = ME11_segment->parametersError().similarityT(ME11_segment->projectionMatrix());
        LocalTrajectoryError error(asSMatrix<5>(mat));
        TrajectoryStateOnSurface TSOS_Segment(param, error, segDet->surface(), &*theService_->magneticField());

        //Set Up Track Propagation
        TSOS_Track = track_at_segment;

        auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
        for (const auto& ch : GEMGeometry_->etaPartitions()){
          if (ch->id().station() != 1) continue; //Only GE1/1
          const BoundPlane& bps(ch->surface());

          TSOS_Segment_On_GEM = propagator->propagate(track_at_segment, ch->surface());
          if (TSOS_Segment_On_GEM.isValid()){
            SegmentGEMDetId = ch->id();
            const GlobalPoint Prop_GEM_GP = TSOS_Segment_On_GEM.globalPosition();
            const LocalPoint Prop_GEM_LP = ch->toLocal(Prop_GEM_GP);
            const LocalPoint Prop_GEM_LP_2D(Prop_GEM_LP.x(), Prop_GEM_LP.y(), 0.0);
            if (((TSOS_Segment.globalPosition().z() * Prop_GEM_GP.z()) > 0) and (bps.bounds().inside(Prop_GEM_LP_2D)) and (ch->id().station() == 1 and ch->id().ring() == 1)){
              cout << "Found a good track prop" << endl;

              float tmp_delta_x = 999.9;
              for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
                GEMDetId gemid((hit)->geographicalId());
                if (!(gemid.det() == DetId::Detector::Muon && gemid.subdetId() == MuonSubdetId::GEM)) continue;
                if (!(gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region())) continue;
                float new_delta_x = abs(hit->localPosition().x() - ch->toLocal(TSOS_Segment_On_GEM.globalPosition()).x());
                if (not (new_delta_x < tmp_delta_x)) continue;
                tmp_delta_x = new_delta_x;
                GEM_Hit_Match_Segment = *hit;
              }
            }
          }
          TSOS_Track_On_GEM = propagator->propagate(TSOS_Track, ch->surface());
          if (TSOS_Track_On_GEM.isValid()){
            TrackGEMDetId = ch->id();
            const GlobalPoint Prop_GEM_GP = TSOS_Track_On_GEM.globalPosition();
            const LocalPoint Prop_GEM_LP = ch->toLocal(Prop_GEM_GP);
            const LocalPoint Prop_GEM_LP_2D(Prop_GEM_LP.x(), Prop_GEM_LP.y(), 0.0);
            if (((TSOS_Segment.globalPosition().z() * Prop_GEM_GP.z()) > 0) and (bps.bounds().inside(Prop_GEM_LP_2D)) and (ch->id().station() == 1 and ch->id().ring() == 1)){
              cout << "Found a good segment prop" << endl;

              float tmp_delta_x = 999.9;
              for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
                GEMDetId gemid((hit)->geographicalId());
                if (!(gemid.det() == DetId::Detector::Muon && gemid.subdetId() == MuonSubdetId::GEM)) continue;
                if (!(gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region())) continue;
                float new_delta_x = abs(hit->localPosition().x() - ch->toLocal(TSOS_Track_On_GEM.globalPosition()).x());
                if (not (new_delta_x < tmp_delta_x)) continue;
                tmp_delta_x = new_delta_x;
                GEM_Hit_Match_Track = *hit;
              }
            }
          }
        }
      }
      //This is where the fill will happen
      /*
      const CSCSegment* ME11_segment;
      TrajectoryStateOnSurface TSOS_Segment;
      TrajectoryStateOnSurface TSOS_Track;
      CSCCorrelatedLCTDigi LCTDigiMatch;
      CSCDetId SegmentCSCDetId;
      CSCDetId TrackCSCDetId;
      CSCDetId LCTDetId;
      TrajectoryStateOnSurface TSOS_Segment_On_GEM;
      TrajectoryStateOnSurface TSOS_Track_On_GEM;
      GEMRecHit GEM_Hit_Match_Segment;
      GEMRecHit GEM_Hit_Match_Track;
      GEMDetId SegmentGEMDetId;
      GEMDetId TrackGEMDetId;
      */
      tree->Fill();
    }
  }
}




void CSCLCTSegmentMatcher::beginJob(){}
void CSCLCTSegmentMatcher::endJob(){}

DEFINE_FWK_MODULE(CSCLCTSegmentMatcher);
