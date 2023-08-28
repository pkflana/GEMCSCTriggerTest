#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
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
#include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCluster.h"
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


struct GEMCSCTriggerData
{
  void init();
  TTree* book(TTree *t);
  //============ Muon Info ================//
  int muon_charge; float muon_pt; float muon_eta;
  unsigned long long evtNum; unsigned long long lumiBlock; int runNum;

  //============ Segment Info =============//
  bool has_segment;
  int segment_CSC_endcap; int segment_CSC_station; 
  int segment_CSC_ring;   int segment_CSC_chamber;
  float segment_CSC_LP_x; float segment_CSC_LP_y; float segment_CSC_LP_z;
  float segment_CSC_GP_x; float segment_CSC_GP_y; float segment_CSC_GP_z;

  bool has_segment_prop1;  bool has_segment_prop2;
  int segment_GE1_region;  int segment_GE2_region;
  int segment_GE1_station; int segment_GE2_station;
  int segment_GE1_ring;    int segment_GE2_ring;
  int segment_GE1_chamber; int segment_GE2_chamber;
  int segment_GE1_layer;   int segment_GE2_layer;
  int segment_GE1_roll;    int segment_GE2_roll;

  //============ Track Info ===============//
  bool has_track;
  int track_CSC_endcap; int track_CSC_station;
  int track_CSC_ring;   int track_CSC_chamber;
  float track_CSC_LP_x; float track_CSC_LP_y; float track_CSC_LP_z;
  float track_CSC_GP_x; float track_CSC_GP_y; float track_CSC_GP_z;

  bool has_track_prop1;  bool has_track_prop2;
  int track_GE1_region;  int track_GE2_region;
  int track_GE1_station; int track_GE2_station;
  int track_GE1_ring;    int track_GE2_ring;
  int track_GE1_chamber; int track_GE2_chamber;
  int track_GE1_layer;   int track_GE2_layer;
  int track_GE1_roll;    int track_GE2_roll;

  //============ LCT Info =================//
  bool has_LCT;
  int LCT_CSC_endcap;    int LCT_CSC_station;
  int LCT_CSC_ring;      int LCT_CSC_chamber;
  bool LCT_CSC_ME1a;     bool LCT_CSC_ME1b;
  int LCT_eighthstrip;   int LCT_slope;
  int LCT_wiregroup;     int LCT_quality;
  int LCT_bend;          int LCT_BX;

  bool has_LCT_prop1;    bool has_LCT_prop2;
  int LCT_GE1_region;    int LCT_GE2_region;
  int LCT_GE1_station;   int LCT_GE2_station;
  int LCT_GE1_ring;      int LCT_GE2_ring;
  int LCT_GE1_chamber;   int LCT_GE2_chamber;
  int LCT_GE1_layer;     int LCT_GE2_layer;
  int LCT_GE1_strip;     int LCT_GE2_strip;
  int LCT_GE1_wiregroup; int LCT_GE2_wiregroup;


};

void GEMCSCTriggerData::init()
{
  //=========== Muon Info ===============//
  float value = 99999;
  muon_charge = value; muon_pt = value; muon_eta = value;
  evtNum = value; lumiBlock = value; runNum = value;

  //============ Segment Info =============//
  has_segment = false;
  segment_CSC_endcap = value; segment_CSC_station = value;
  segment_CSC_ring = value; segment_CSC_chamber = value;
  segment_CSC_LP_x = value; segment_CSC_LP_y = value; segment_CSC_LP_z = value;
  segment_CSC_GP_x = value; segment_CSC_GP_y = value; segment_CSC_GP_z = value;

  //============ Segment Prop =============//
  has_segment_prop1 = false;   has_segment_prop2 = false;
  segment_GE1_region = value;  segment_GE2_region = value;
  segment_GE1_station = value; segment_GE2_station = value;
  segment_GE1_ring = value;    segment_GE2_ring = value;
  segment_GE1_chamber = value; segment_GE2_chamber = value;
  segment_GE1_layer = value;   segment_GE2_layer = value;
  segment_GE1_roll = value;    segment_GE2_roll = value;

  //============ Track Info ===============//
  has_track = false;
  track_CSC_endcap = value; track_CSC_station = value;
  track_CSC_ring = value; track_CSC_chamber = value;
  track_CSC_LP_x = value; track_CSC_LP_y = value; track_CSC_LP_z = value;
  track_CSC_GP_x = value; track_CSC_GP_y = value; track_CSC_GP_z = value;

  //============ Track Prop ===============//
  has_track_prop1 = false; has_track_prop2 = false;
  track_GE1_region = value;  track_GE2_region = value;
  track_GE1_station = value; track_GE2_station = value;
  track_GE1_ring = value;    track_GE2_ring = value;
  track_GE1_chamber = value; track_GE2_chamber = value;
  track_GE1_layer = value;   track_GE2_layer = value;
  track_GE1_roll = value;    track_GE2_roll = value;

  //============ LCT Info =================//
  has_LCT = false;
  LCT_CSC_endcap = value; LCT_CSC_station = value;
  LCT_CSC_ring = value; LCT_CSC_chamber = value;
  LCT_CSC_ME1a = false; LCT_CSC_ME1b = false;
  LCT_eighthstrip = value; LCT_slope = value;
  LCT_wiregroup = value; LCT_quality = value;
  LCT_bend = value; LCT_BX = value;

  //============ LCT Prop =================//
  has_LCT_prop1 = false; has_LCT_prop2 = false;
  LCT_GE1_region = value;  LCT_GE2_region = value;
  LCT_GE1_station = value; LCT_GE2_station = value;
  LCT_GE1_ring = value;    LCT_GE2_ring = value;
  LCT_GE1_chamber = value; LCT_GE2_chamber = value;
  LCT_GE1_layer = value;   LCT_GE2_layer = value;
  LCT_GE1_strip = value;   LCT_GE2_strip = value;


}

TTree* GEMCSCTriggerData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("GEM_CSC_Trigger", "GEM_CSC_Trigger");

  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt); t->Branch("muon_eta", &muon_eta);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("runNum", &runNum);

  //============ Segment Info =============//

  return t;
}

class GEMCSCTriggerTester : public edm::one::EDAnalyzer<> {
public:
  explicit GEMCSCTriggerTester(const edm::ParameterSet&);
  ~GEMCSCTriggerTester(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void startMuonAnalysis(const reco::Muon* mu);
  void findSegmentTrackLCT(const TrackingRecHit* RecHit, const reco::Track* Track, CSCSegment*& ME11_Segment, TrajectoryStateOnSurface& Track_At_Segment);
  void propagateSegmentTrackLCT(const reco::Track* Track, CSCSegment* ME11_Segment, TrajectoryStateOnSurface Track_At_Segment, TrajectoryStateOnSurface& TSOS_Segment_On_GE1, TrajectoryStateOnSurface& TSOS_Segment_On_GE2, TrajectoryStateOnSurface& TSOS_Track_On_GE1, TrajectoryStateOnSurface& TSOS_Track_On_GE2, GEMDetId& TSOS_Segment_On_GE1_ch, GEMDetId& TSOS_Segment_On_GE2_ch, GEMDetId& TSOS_Track_On_GE1_ch, GEMDetId& TSOS_Track_On_GE2_ch);

  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> co_token;
  edm::Handle<CSCCorrelatedLCTDigiCollection> correlatedlcts;
  edm::EDGetTokenT<GEMPadDigiClusterCollection> gemdigi_token;
  edm::Handle<GEMPadDigiClusterCollection> gemPadDigis;

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

  GEMCSCTriggerData data_;
  TTree* tree;

  bool isMC;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;

  //Maps for the CSC LCT Slope Extrapolation (CSC LCT -> GEM Layers)
  //Split by Even/Odd, by L1/L2, and by ME11a/ME11b (8 cases)
  string SlopeExtrapolationME11aEvenL1Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11a_even_GEMlayer1.txt";
  map<int, int> SlopeExtrapolationME11aEvenL1_Map;
  string SlopeExtrapolationME11bEvenL1Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11b_even_GEMlayer1.txt";
  map<int, int> SlopeExtrapolationME11bEvenL1_Map;
  string SlopeExtrapolationME11aOddL1Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11a_odd_GEMlayer1.txt";
  map<int, int> SlopeExtrapolationME11aOddL1_Map;
  string SlopeExtrapolationME11bOddL1Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11b_odd_GEMlayer1.txt";
  map<int, int> SlopeExtrapolationME11bOddL1_Map;

  string SlopeExtrapolationME11aEvenL2Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11a_even_GEMlayer2.txt";
  map<int, int> SlopeExtrapolationME11aEvenL2_Map;
  string SlopeExtrapolationME11bEvenL2Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11b_even_GEMlayer2.txt";
  map<int, int> SlopeExtrapolationME11bEvenL2_Map;
  string SlopeExtrapolationME11aOddL2Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11a_odd_GEMlayer2.txt";
  map<int, int> SlopeExtrapolationME11aOddL2_Map;
  string SlopeExtrapolationME11bOddL2Name = "../luts/GEMCSC/SlopeCorrection/FacingChambers/ExtrapolationBySlope_ME11b_odd_GEMlayer2.txt";
  map<int, int> SlopeExtrapolationME11bOddL2_Map;

  //Maps for the GEM Pad Digi Clusters to be converted into CSC eighth strip units
  //Split by Even/Odd, and by ME11a/ME11b (4 cases)
  string GEMPadDigiToCSCEightStripME11aEvenName = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_pad_es_ME1a_even.txt";
  map<int, int> GEMPadDigiToCSCEigthStripME11aEven_Map;
  string GEMPadDigiToCSCEightStripME11bEvenName = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_pad_es_ME1b_even.txt";
  map<int, int> GEMPadDigiToCSCEigthStripME11bEven_Map;
  string GEMPadDigiToCSCEightStripME11aOddName = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_pad_es_ME1a_odd.txt";
  map<int, int> GEMPadDigiToCSCEigthStripME11aOdd_Map;
  string GEMPadDigiToCSCEightStripME11bOddName = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_pad_es_ME1b_odd.txt";
  map<int, int> GEMPadDigiToCSCEigthStripME11bOdd_Map;

  //Maps for the GEM Pad Digi Clusters to be converted into CSC Min and Max WireGroups
  //Split by Even/Odd, Layer1/Layer2, and by Min/Max
  string GEMPadDigiToCSCWGMinEvenL1Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l1_min_wg_ME11_even.txt";
  map<int, int> GEMPadDigiToCSCWGMinEvenL1_Map;
  string GEMPadDigiToCSCWGMaxEvenL1Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l1_max_wg_ME11_even.txt";
  map<int, int> GEMPadDigiToCSCWGMaxEvenL1_Map;

  string GEMPadDigiToCSCWGMinOddL1Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l1_min_wg_ME11_odd.txt";
  map<int, int> GEMPadDigiToCSCWGMinOddL1_Map;
  string GEMPadDigiToCSCWGMaxOddL1Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l1_max_wg_ME11_odd.txt";
  map<int, int> GEMPadDigiToCSCWGMaxOddL1_Map;

  string GEMPadDigiToCSCWGMinEvenL2Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l2_min_wg_ME11_even.txt";
  map<int, int> GEMPadDigiToCSCWGMinEvenL2_Map;
  string GEMPadDigiToCSCWGMaxEvenL2Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l2_max_wg_ME11_even.txt";
  map<int, int> GEMPadDigiToCSCWGMaxEvenL2_Map;

  string GEMPadDigiToCSCWGMinOddL2Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l2_min_wg_ME11_odd.txt";
  map<int, int> GEMPadDigiToCSCWGMinOddL2_Map;
  string GEMPadDigiToCSCWGMaxOddL2Name = "../luts/GEMCSC/CoordinateConversion/GEMCSCLUT_roll_l2_max_wg_ME11_odd.txt";
  map<int, int> GEMPadDigiToCSCWGMaxOddL2_Map;

};


GEMCSCTriggerTester::GEMCSCTriggerTester(const edm::ParameterSet& iConfig)
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
  gemdigi_token = consumes<GEMPadDigiClusterCollection>(iConfig.getParameter<edm::InputTag>("gemPadDigiCluster"));

  debug = iConfig.getParameter<bool>("debug");
  std::cout << "debug " << debug << std::endl;

  tree = data_.book(tree);

}


void
GEMCSCTriggerTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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
  iEvent.getByToken(co_token, correlatedlcts);
  iEvent.getByToken(gemdigi_token, gemPadDigis);

  if (debug) cout << "New! EventNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (debug) cout << "There are " << muons->size() << " muons in event" << endl;
  if (muons->size() == 0) return;

  edm::Handle<CSCSegmentCollection> cscSegments;
  if (! iEvent.getByToken(cscSegments_, cscSegments)){std::cout << "Bad segments" << std::endl;}

  if (debug) cout << "New Event" << endl;


  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    if (not mu->isGlobalMuon()) continue;

    //Lets find the ME1/1 segments from this muon
    if (!(mu->isStandAloneMuon())) continue;
    if (debug) cout << "New Muon" << endl;
    startMuonAnalysis(mu);
  }
}


void GEMCSCTriggerTester::startMuonAnalysis(const reco::Muon* mu){
  const reco::Track* Track = mu->outerTrack().get();
  if (Track->validFraction() > 0.0) return;

  //Find ME11 Segment from Track
  for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
    data_.init();
    data_.muon_charge = mu->charge();
    const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();

    CSCSegment* ME11_Segment;
    TrajectoryStateOnSurface Track_At_Segment = TrajectoryStateOnSurface();

    findSegmentTrackLCT(RecHit, Track, ME11_Segment, Track_At_Segment);
    if ((!data_.has_segment) or (!data_.has_track) or (!data_.has_LCT)) continue;
    cout << "Results of findSegmentTrackLCT S:T:L " << data_.has_segment << ":" << data_.has_track << ":" << data_.has_LCT << endl;

    TrajectoryStateOnSurface TSOS_Segment_On_GE1 = TrajectoryStateOnSurface();
    TrajectoryStateOnSurface TSOS_Segment_On_GE2 = TrajectoryStateOnSurface();
    TrajectoryStateOnSurface TSOS_Track_On_GE1 = TrajectoryStateOnSurface();
    TrajectoryStateOnSurface TSOS_Track_On_GE2 = TrajectoryStateOnSurface();

    GEMDetId TSOS_Segment_On_GE1_ch = GEMDetId();
    GEMDetId TSOS_Segment_On_GE2_ch = GEMDetId();
    GEMDetId TSOS_Track_On_GE1_ch = GEMDetId();
    GEMDetId TSOS_Track_On_GE2_ch = GEMDetId();

    propagateSegmentTrackLCT(Track, ME11_Segment, Track_At_Segment, TSOS_Segment_On_GE1, TSOS_Segment_On_GE2, TSOS_Track_On_GE1, TSOS_Track_On_GE2, TSOS_Segment_On_GE1_ch, TSOS_Segment_On_GE2_ch, TSOS_Track_On_GE1_ch, TSOS_Track_On_GE2_ch);

    
  } 
}


void GEMCSCTriggerTester::findSegmentTrackLCT(const TrackingRecHit* RecHit, const reco::Track* Track, CSCSegment*& ME11_Segment, TrajectoryStateOnSurface& Track_At_Segment){
  //if (debug) cout << "Starting findSegmentTrackLCT" << endl; //Even for debug this print out is obnoxious
  DetId RecHitId = RecHit->geographicalId();
  uint16_t RecHitDetId = RecHitId.det();
  if (RecHitDetId != DetId::Muon) return;
  uint16_t RecHitSubDet = RecHitId.subdetId();
  if (RecHitSubDet != (uint16_t)MuonSubdetId::CSC) return;
  CSCDetId SegmentCSCDetId = CSCDetId(RecHitId);
  if (!(SegmentCSCDetId.station() == 1 and SegmentCSCDetId.ring() == 1 and RecHit->dimension() == 4)) return;
  //Found a rechit on Muon Detector, CSC, Station/Ring 1/1, and Segment (dimensions 4)
  RecSegment* Rec_Segment = (RecSegment*)RecHit;
  ME11_Segment = (CSCSegment*)Rec_Segment;
  CSCDetId SegmentCSCDetIdL4 = CSCDetId(SegmentCSCDetId.endcap(), SegmentCSCDetId.station(), SegmentCSCDetId.ring(), SegmentCSCDetId.chamber(), 4);
  //Build a new CSCDetId at layer 4 (key layer)
  const CSCLayer* ME11_Layer = CSCGeometry_->layer(SegmentCSCDetIdL4);
  const CSCLayerGeometry* ME11_Layer_Geo = ME11_Layer->geometry();
  reco::TransientTrack TTrack = ttrackBuilder_->build(Track);
  LocalPoint ME11_Segment_LP = ME11_Segment->localPosition();
  GlobalPoint ME11_Segment_GP = ME11_Layer->toGlobal(ME11_Segment_LP);

  data_.has_segment = true;
  data_.segment_CSC_endcap = SegmentCSCDetId.endcap();
  data_.segment_CSC_station = SegmentCSCDetId.station();
  data_.segment_CSC_ring = SegmentCSCDetId.ring();
  data_.segment_CSC_chamber = SegmentCSCDetId.chamber();
  data_.segment_CSC_LP_x = ME11_Segment_LP.x();
  data_.segment_CSC_LP_y = ME11_Segment_LP.y();
  data_.segment_CSC_LP_z = ME11_Segment_LP.z();
  data_.segment_CSC_GP_x = ME11_Segment_GP.x();
  data_.segment_CSC_GP_y = ME11_Segment_GP.y();
  data_.segment_CSC_GP_z = ME11_Segment_GP.z();
  if (debug) cout << "Found an ME11 Segment" << endl;

  //Find the Track at this Segment Position
  Track_At_Segment = TTrack.stateOnSurface(ME11_Segment_GP);
  if (!(Track_At_Segment.isValid())) return;

  GlobalPoint Track_At_Segment_GP = Track_At_Segment.globalPosition();
  LocalPoint Track_At_Segment_LP = ME11_Layer->toLocal(Track_At_Segment_GP);
  data_.has_track = true;
  data_.track_CSC_endcap = SegmentCSCDetId.endcap();
  data_.track_CSC_station = SegmentCSCDetId.station();
  data_.track_CSC_ring = SegmentCSCDetId.ring();
  data_.track_CSC_chamber = SegmentCSCDetId.chamber();
  data_.track_CSC_LP_x = Track_At_Segment_GP.x();
  data_.track_CSC_LP_y = Track_At_Segment_GP.y();
  data_.track_CSC_LP_z = Track_At_Segment_GP.z();
  data_.track_CSC_GP_x = Track_At_Segment_GP.x();
  data_.track_CSC_GP_y = Track_At_Segment_GP.y();
  data_.track_CSC_GP_z = Track_At_Segment_GP.z();
  if (debug) cout << "Found a matching Track" << endl;
  if (debug) cout << "Segment LP : GP " << ME11_Segment_LP << " : " << ME11_Segment_GP << endl;
  if (debug) cout << "Track   LP : GP " << Track_At_Segment_LP << " : " << Track_At_Segment_LP << endl;

  //Prepare FracStrips for LCT Matching
  float SegmentFracStrip = ME11_Layer_Geo->strip(ME11_Segment->localPosition());
  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j = correlatedlcts->begin(); j != correlatedlcts->end(); j++){
    CSCDetId LCTDetId = (*j).first;
    if (!(LCTDetId.station() == 1 and LCTDetId.ring() == 1)) continue;
    if (!(LCTDetId == SegmentCSCDetId)) continue;
    data_.LCT_CSC_endcap = LCTDetId.endcap();
    data_.LCT_CSC_station = LCTDetId.station();
    data_.LCT_CSC_ring = LCTDetId.ring();
    data_.LCT_CSC_chamber = LCTDetId.chamber();
    data_.LCT_CSC_ME1a = LCTDetId.isME1a();
    data_.LCT_CSC_ME1b = LCTDetId.isME1b();
    if (debug) cout << "Found a LCTDet to SegmentDet Match on Det " << LCTDetId << endl;
    //Loop over LCT Clusters on the Chamber and find best match
    std::vector<CSCCorrelatedLCTDigi>::const_iterator CSCCorrLCT = (*j).second.first;
    std::vector<CSCCorrelatedLCTDigi>::const_iterator last = (*j).second.second;
    CSCCorrelatedLCTDigi LCTDigiMatch = CSCCorrelatedLCTDigi();
    for (; CSCCorrLCT != last; ++CSCCorrLCT){
      float LCTFracStrip = CSCCorrLCT->getFractionalStrip();
      if (abs(SegmentFracStrip - LCTFracStrip) < abs(SegmentFracStrip - LCTDigiMatch.getFractionalStrip())){
        //This is the best LCT Match so far
        LCTDigiMatch = *CSCCorrLCT;
        data_.has_LCT = true;
        if (debug) cout << "Found a matching LCT (checking for better one)" << endl;
        data_.LCT_eighthstrip = (CSCCorrLCT->getStrip())*4 + (CSCCorrLCT->getQuartStripBit())*2 + (CSCCorrLCT->getEighthStripBit());
        data_.LCT_slope = CSCCorrLCT->getSlope();
        data_.LCT_wiregroup = CSCCorrLCT->getKeyWG();
        data_.LCT_quality = CSCCorrLCT->getQuality();
        data_.LCT_bend = CSCCorrLCT->getBend();
        data_.LCT_BX = CSCCorrLCT->getBX();
      }
    }
  }
}


void GEMCSCTriggerTester::propagateSegmentTrackLCT(const reco::Track* Track, CSCSegment* ME11_Segment, TrajectoryStateOnSurface Track_At_Segment, TrajectoryStateOnSurface& TSOS_Segment_On_GE1, TrajectoryStateOnSurface& TSOS_Segment_On_GE2, TrajectoryStateOnSurface& TSOS_Track_On_GE1, TrajectoryStateOnSurface& TSOS_Track_On_GE2, GEMDetId& TSOS_Segment_On_GE1_ch, GEMDetId& TSOS_Segment_On_GE2_ch, GEMDetId& TSOS_Track_On_GE1_ch, GEMDetId& TSOS_Track_On_GE2_ch){
  if (debug) cout << "Starting propagateSegmentTrackLCT" << endl;
  //Set up Segment propagation
  DetId segDetId = ME11_Segment->geographicalId();
  const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId());
  LocalVector momentum_at_surface = (Track->outerP()) * (ME11_Segment->localDirection());
  LocalTrajectoryParameters param(ME11_Segment->localPosition(), momentum_at_surface, data_.muon_charge);
  AlgebraicSymMatrix mat(5,0);
  mat = ME11_Segment->parametersError().similarityT(ME11_Segment->projectionMatrix());
  LocalTrajectoryError error(asSMatrix<5>(mat));
  TrajectoryStateOnSurface TSOS_Segment(param, error, segDet->surface(), &*theService_->magneticField());

  //Set up Track propagation
  TrajectoryStateOnSurface TSOS_Track = Track_At_Segment;


  if (debug) cout << "Found Starting Trajectories ***NOTE Track Local Coords are not the same!!!***" << endl;
  if (debug) cout << "Seg TSOS   LP : GP : LocalDir " << TSOS_Segment.localPosition() << " : " << TSOS_Segment.globalPosition() << " : " << TSOS_Segment.localDirection() << endl;
  if (debug) cout << "Track TSOS LP : GP : LocalDir " << TSOS_Track.localPosition() << " : " << TSOS_Track.globalPosition() << " : " << TSOS_Track.localDirection() << endl;

  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  for (const auto& ch : GEMGeometry_->etaPartitions()){
    if (!(ch->id().station() == 1 and ch->id().ring() == 1)) continue; //Only look at GE1/1
    CSCDetId SegmentDetId = CSCDetId(segDetId);
    if (((SegmentDetId.endcap() == 2 and ch->id().region() == -1) or (SegmentDetId.endcap() == 1 and ch->id().region() == 1)) and (SegmentDetId.station() == ch->id().station()) and (SegmentDetId.ring() == ch->id().ring()) and (SegmentDetId.chamber() == ch->id().chamber())) continue; //Only look at same detector as segment
    bool good_segment_prop = false;
    const BoundPlane& bps(ch->surface());
    TrajectoryStateOnSurface TSOS_Segment_On_GEM = propagator->propagate(TSOS_Segment, ch->surface());
    if (TSOS_Segment_On_GEM.isValid()){
      const GlobalPoint Prop_GEM_GP = TSOS_Segment_On_GEM.globalPosition();
      const LocalPoint Prop_GEM_LP = TSOS_Segment_On_GEM.localPosition();
      const LocalPoint Prop_GEM_LP_2D(Prop_GEM_LP.x(), Prop_GEM_LP.y(), 0.0);
      if (((TSOS_Segment.globalPosition().z() * Prop_GEM_GP.z()) > 0.0) and (bps.bounds().inside(Prop_GEM_LP_2D))){
        good_segment_prop = true;
      }
    }
    bool good_track_prop = false;
    TrajectoryStateOnSurface TSOS_Track_On_GEM = propagator->propagate(TSOS_Track, ch->surface());
    if (TSOS_Track_On_GEM.isValid()){
      const GlobalPoint Prop_GEM_GP = TSOS_Track_On_GEM.globalPosition();
      const LocalPoint Prop_GEM_LP = TSOS_Track_On_GEM.localPosition();
      const LocalPoint Prop_GEM_LP_2D(Prop_GEM_LP.x(), Prop_GEM_LP.y(), 0.0);
      if (((TSOS_Track.globalPosition().z() * Prop_GEM_GP.z()) > 0.0) and (bps.bounds().inside(Prop_GEM_LP_2D))){
        good_track_prop = true;
      }
    }
    //End of chamber loop, lets check if props were successful and fill correct references
    int ch_region = ch->id().region();
    int ch_station = ch->id().station();
    int ch_ring = ch->id().ring();
    int ch_chamber = ch->id().chamber();
    int ch_layer = ch->id().layer();
    int ch_roll = ch->id().roll();
    if (ch->id().layer() == 1){
      if (good_segment_prop){
        data_.has_segment_prop1 = true;
        TSOS_Segment_On_GE1 = TSOS_Segment_On_GEM;
        if (debug) cout << "Found segprop at L1 LP/GP " << TSOS_Segment_On_GEM.localPosition() << "/" << TSOS_Segment_On_GEM.globalPosition() << endl;
        data_.segment_GE1_region = ch_region;
        data_.segment_GE1_station = ch_station;
        data_.segment_GE1_ring = ch_ring;
        data_.segment_GE1_chamber = ch_chamber;
        data_.segment_GE1_layer = ch_layer;
        data_.segment_GE1_roll = ch_roll;
      }
      if (good_track_prop){
        data_.has_track_prop1 = true;
        TSOS_Track_On_GE1 = TSOS_Track_On_GEM;
        if (debug) cout << "Found traprop at L1 LP/GP " << TSOS_Track_On_GEM.localPosition() << "/" << TSOS_Track_On_GEM.globalPosition() << endl;
        data_.track_GE1_region = ch_region;
        data_.track_GE1_station = ch_station;
        data_.track_GE1_ring = ch_ring;
        data_.track_GE1_chamber = ch_chamber;
        data_.track_GE1_layer = ch_layer;
        data_.track_GE1_roll = ch_roll;
      }
    }
    if (ch->id().layer() == 2){
      if (good_segment_prop){
        data_.has_segment_prop2 = true;
        TSOS_Segment_On_GE2 = TSOS_Segment_On_GEM;
        if (debug) cout << "Found segprop at L2 LP/GP " << TSOS_Segment_On_GEM.localPosition() << "/" << TSOS_Segment_On_GEM.globalPosition() << endl;
        data_.segment_GE2_region = ch_region;
        data_.segment_GE2_station = ch_station;
        data_.segment_GE2_ring = ch_ring;
        data_.segment_GE2_chamber = ch_chamber;
        data_.segment_GE2_layer = ch_layer;
        data_.segment_GE2_roll = ch_roll;
      }
      if (good_track_prop){
        data_.has_track_prop2 = true;
        TSOS_Track_On_GE2 = TSOS_Track_On_GEM;
        if (debug) cout << "Found traprop at L2 LP/GP " << TSOS_Track_On_GEM.localPosition() << "/" << TSOS_Track_On_GEM.globalPosition() << endl;
        data_.track_GE2_region = ch_region;
        data_.track_GE2_station = ch_station;
        data_.track_GE2_ring = ch_ring;
        data_.track_GE2_chamber = ch_chamber;
        data_.track_GE2_layer = ch_layer;
        data_.track_GE2_roll = ch_roll;
      }
    }
  }
  //Find LCT Propagation
  int LCT_eighth_strip = data_.LCT_eighthstrip;
  int LCT_slope = data_.LCT_slope;
  int LCT_bend = data_.LCT_bend;
  if (debug) cout << "Found the LCT Eighth Strip at " << LCT_eighth_strip << endl;
  if (debug) cout << "And slope at " << LCT_slope << endl;
  map<int, int> CSCLCTPropL1Map;
  map<int, int> CSCLCTPropL2Map;
  bool LCTChamber_Even = (data_.LCT_CSC_chamber % 2 == 0);
  if (data_.LCT_CSC_ME1a){
    CSCLCTPropL1Map = (LCTChamber_Even) ? SlopeExtrapolationME11aEvenL1_Map : SlopeExtrapolationME11aOddL1_Map;
    CSCLCTPropL2Map = (LCTChamber_Even) ? SlopeExtrapolationME11aEvenL2_Map : SlopeExtrapolationME11aOddL2_Map;
  }
  if (data_.LCT_CSC_ME1b){
    CSCLCTPropL1Map = (LCTChamber_Even) ? SlopeExtrapolationME11bEvenL1_Map : SlopeExtrapolationME11bOddL1_Map;
    CSCLCTPropL2Map = (LCTChamber_Even) ? SlopeExtrapolationME11bEvenL2_Map : SlopeExtrapolationME11bOddL2_Map;
  }
  int slope_propagationL1 = CSCLCTPropL1Map[LCT_slope];
  int slope_propagationL2 = CSCLCTPropL2Map[LCT_slope];
  if (debug) cout << "LCT Eighth Strip Prop Adjust L1/L2 = " << slope_propagationL1 << "/" << slope_propagationL2 << endl;
  int LCTToGEML1EighthStrip = LCT_eighth_strip + slope_propagationL1*((LCT_bend*2)-1);
  int LCTToGEML2EighthStrip = LCT_eighth_strip + slope_propagationL2*((LCT_bend*2)-1);
  if (debug) cout << "LCT Prop from " << LCT_eighth_strip << " to L1 " << LCTToGEML1EighthStrip << " and L2 " << LCTToGEML2EighthStrip << endl;
  data_.LCT_GE1_region = (data_.LCT_CSC_endcap*2)-1;  data_.LCT_GE2_region = (data_.LCT_CSC_endcap*2)-1;
  data_.LCT_GE1_station = data_.LCT_CSC_station; data_.LCT_GE2_station = data_.LCT_CSC_station;
  data_.LCT_GE1_ring = data_.LCT_CSC_ring;    data_.LCT_GE2_ring = data_.LCT_CSC_ring;
  data_.LCT_GE1_chamber = data_.LCT_CSC_chamber; data_.LCT_GE2_chamber = data_.LCT_CSC_chamber;
  data_.LCT_GE1_layer = 1;   data_.LCT_GE2_layer = 2;
  data_.LCT_GE1_strip = LCTToGEML1EighthStrip;   data_.LCT_GE2_strip = LCTToGEML2EighthStrip;
}


void GEMCSCTriggerTester::beginJob(){
  //Lets make the SlopeExtrapolationLUTMaps
  cout << "Begin job!" << endl;

  string delimiter = " ";
  string line;

  ifstream SlopeExtrapolationME11aEvenL1File;
  SlopeExtrapolationME11aEvenL1File.open(SlopeExtrapolationME11aEvenL1Name);
  if (SlopeExtrapolationME11aEvenL1File.is_open()){
    while(getline(SlopeExtrapolationME11aEvenL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11aEvenL1_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11bEvenL1File;
  SlopeExtrapolationME11bEvenL1File.open(SlopeExtrapolationME11bEvenL1Name);
  if (SlopeExtrapolationME11bEvenL1File.is_open()){
    while(getline(SlopeExtrapolationME11bEvenL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11bEvenL1_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11aOddL1File;
  SlopeExtrapolationME11aOddL1File.open(SlopeExtrapolationME11aOddL1Name);
  if (SlopeExtrapolationME11aOddL1File.is_open()){
    while(getline(SlopeExtrapolationME11aOddL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11aOddL1_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11bOddL1File;
  SlopeExtrapolationME11bOddL1File.open(SlopeExtrapolationME11bOddL1Name);
  if (SlopeExtrapolationME11bOddL1File.is_open()){
    while(getline(SlopeExtrapolationME11bOddL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11bOddL1_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11aEvenL2File;
  SlopeExtrapolationME11aEvenL2File.open(SlopeExtrapolationME11aEvenL2Name);
  if (SlopeExtrapolationME11aEvenL2File.is_open()){
    while(getline(SlopeExtrapolationME11aEvenL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11aEvenL2_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11bEvenL2File;
  SlopeExtrapolationME11bEvenL2File.open(SlopeExtrapolationME11bEvenL2Name);
  if (SlopeExtrapolationME11bEvenL2File.is_open()){
    while(getline(SlopeExtrapolationME11bEvenL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11bEvenL2_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11aOddL2File;
  SlopeExtrapolationME11aOddL2File.open(SlopeExtrapolationME11aOddL2Name);
  if (SlopeExtrapolationME11aOddL2File.is_open()){
    while(getline(SlopeExtrapolationME11aOddL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11aOddL2_Map[key] = value;
      }
    }
  }

  ifstream SlopeExtrapolationME11bOddL2File;
  SlopeExtrapolationME11bOddL2File.open(SlopeExtrapolationME11bOddL2Name);
  if (SlopeExtrapolationME11bOddL2File.is_open()){
    while(getline(SlopeExtrapolationME11bOddL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        SlopeExtrapolationME11bOddL2_Map[key] = value;
      }
    }
  }

  //GEMPadDigi to CSC Eighth Strip LUTs
  ifstream GEMPadDigiToCSCEightStripME11aEvenFile;
  GEMPadDigiToCSCEightStripME11aEvenFile.open(GEMPadDigiToCSCEightStripME11aEvenName);
  if (GEMPadDigiToCSCEightStripME11aEvenFile.is_open()){
    while(getline(GEMPadDigiToCSCEightStripME11aEvenFile, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCEigthStripME11aEven_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCEightStripME11bEvenFile;
  GEMPadDigiToCSCEightStripME11bEvenFile.open(GEMPadDigiToCSCEightStripME11bEvenName);
  if (GEMPadDigiToCSCEightStripME11bEvenFile.is_open()){
    while(getline(GEMPadDigiToCSCEightStripME11bEvenFile, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCEigthStripME11bEven_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCEightStripME11aOddFile;
  GEMPadDigiToCSCEightStripME11aOddFile.open(GEMPadDigiToCSCEightStripME11aOddName);
  if (GEMPadDigiToCSCEightStripME11aOddFile.is_open()){
    while(getline(GEMPadDigiToCSCEightStripME11aOddFile, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCEigthStripME11aOdd_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCEightStripME11bOddFile;
  GEMPadDigiToCSCEightStripME11bOddFile.open(GEMPadDigiToCSCEightStripME11bOddName);
  if (GEMPadDigiToCSCEightStripME11bOddFile.is_open()){
    while(getline(GEMPadDigiToCSCEightStripME11bOddFile, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCEigthStripME11bOdd_Map[key] = value;
      }
    }
  }

  //GEMPadDigi to CSC WG Min/Max LUTs
  ifstream GEMPadDigiToCSCWGMinEvenL1File;
  GEMPadDigiToCSCWGMinEvenL1File.open(GEMPadDigiToCSCWGMinEvenL1Name);
  if (GEMPadDigiToCSCWGMinEvenL1File.is_open()){
    while(getline(GEMPadDigiToCSCWGMinEvenL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMinEvenL1_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMaxEvenL1File;
  GEMPadDigiToCSCWGMaxEvenL1File.open(GEMPadDigiToCSCWGMaxEvenL1Name);
  if (GEMPadDigiToCSCWGMaxEvenL1File.is_open()){
    while(getline(GEMPadDigiToCSCWGMaxEvenL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMaxEvenL1_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMinOddL1File;
  GEMPadDigiToCSCWGMinOddL1File.open(GEMPadDigiToCSCWGMinOddL1Name);
  if (GEMPadDigiToCSCWGMinOddL1File.is_open()){
    while(getline(GEMPadDigiToCSCWGMinOddL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMinOddL1_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMaxOddL1File;
  GEMPadDigiToCSCWGMaxOddL1File.open(GEMPadDigiToCSCWGMaxOddL1Name);
  if (GEMPadDigiToCSCWGMaxOddL1File.is_open()){
    while(getline(GEMPadDigiToCSCWGMaxOddL1File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMaxOddL1_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMinEvenL2File;
  GEMPadDigiToCSCWGMinEvenL2File.open(GEMPadDigiToCSCWGMinEvenL2Name);
  if (GEMPadDigiToCSCWGMinEvenL2File.is_open()){
    while(getline(GEMPadDigiToCSCWGMinEvenL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMinEvenL2_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMaxEvenL2File;
  GEMPadDigiToCSCWGMaxEvenL2File.open(GEMPadDigiToCSCWGMaxEvenL2Name);
  if (GEMPadDigiToCSCWGMaxEvenL2File.is_open()){
    while(getline(GEMPadDigiToCSCWGMaxEvenL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMaxEvenL2_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMinOddL2File;
  GEMPadDigiToCSCWGMinOddL2File.open(GEMPadDigiToCSCWGMinOddL2Name);
  if (GEMPadDigiToCSCWGMinOddL2File.is_open()){
    while(getline(GEMPadDigiToCSCWGMinOddL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMinOddL2_Map[key] = value;
      }
    }
  }

  ifstream GEMPadDigiToCSCWGMaxOddL2File;
  GEMPadDigiToCSCWGMaxOddL2File.open(GEMPadDigiToCSCWGMaxOddL2Name);
  if (GEMPadDigiToCSCWGMaxOddL2File.is_open()){
    while(getline(GEMPadDigiToCSCWGMaxOddL2File, line)){
      int key = atoi(line.substr(0, line.find(delimiter)).c_str());
      int value = atoi(line.substr(line.find(delimiter), -1).c_str());
      if (!((key == 0) and (value == 0))){
        GEMPadDigiToCSCWGMaxOddL2_Map[key] = value;
      }
    }
  }



  cout << "Created all LUTs" << endl;
  cout << "Ended Begin Job, starting Event Loop" << endl;
}
void GEMCSCTriggerTester::endJob(){}

DEFINE_FWK_MODULE(GEMCSCTriggerTester);
