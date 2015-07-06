////////////////////////////////////////////////
// rpc_tuplizer.h
//
// Creates the class rpc_tuplizer
// Works in L1Tmuon Framework
//
// Imported by rpc_tuplizer.cc
//
// Written by David Curry
//
///////////////////////////////////////////////

#ifndef rpc_tuplizer_h
#define rpc_tuplizer_h

// Includes
#include <stdlib.h>
#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include <typeinfo>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "L1TriggerDPGUpgrade/L1TMuon/interface/SubsystemCollectorFactory.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitive.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonTriggerPrimitiveFwd.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrack.h"
#include "L1TriggerDPGUpgrade/DataFormats/interface/L1TMuonInternalTrackFwd.h"

//CSCTF
#include "DataFormats/L1CSCTrackFinder/interface/L1CSCTrackCollection.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerPtScale.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerPtScaleRcd.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "PtAddress.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"



 
// Useful Definitions
#define MAX_CSCTF_TRK 36    // max # of CSCTF tracks per BX
#define MAX_LCTS_PER_TRK 4  // max # of LCTS which form a CSCTF track
#define MAX_RPC_LCTS_EVT 36 // 
#define MAX_MUONS 100       // max# of reco muons in event
#define MAX_SEGS_STD 16

using namespace std;
using namespace edm;
using namespace L1TMuon;
using namespace reco;

typedef ParameterSet PSet;
typedef EventSetup iSetup;


// Create the main class that wraps up the ntuplizer
class rpc_tuplizer: public edm::EDAnalyzer {

 public:
  explicit rpc_tuplizer(const edm::ParameterSet&);
  ~rpc_tuplizer(void); 
  
  void analyze(const edm::Event&, const edm::EventSetup&);

  int printLevel;

  
 private:
  
  //------- Specify all inputs the class will need to make ntuple-------
  std::string outputFile;
  edm::InputTag csctfTag;
  edm::InputTag rpcTPTag;
  edm::ParameterSet LUTparam;
  edm::InputTag genSrc;
  edm::InputTag muonsTag;
  edm::InputTag cscSegTag;


  int isMC;

  //Define the file to be filled
  TFile *file;
  TTree *CSCtree;
  //TTree *RPCtree;

  // Needed for CSCTF LUTs
  CSCSectorReceiverLUT* srLUTs_[5][2];
  const L1MuTriggerScales  *scale;
  const L1MuTriggerPtScale *ptScale;  

  // ptLUT. Input: packed pt/eta bit from track converter
  const float ptscale[33] = {
    -1.,   0.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
    4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,
    16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,
    50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 1.E6 };

  const float etabins[16] =
    {0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
     1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};


  //------------Create methods to fill tree------------------
  void fillCSCTF(const edm::Handle< vector<L1TMuon::InternalTrack> > tracks,
		 const edm::Handle< vector<L1TMuon::TriggerPrimitive> > hits,
		 //const edm::Handle< vector<reco::GenParticle> > genps,
		 const L1MuTriggerScales  *scales,
		 const L1MuTriggerPtScale *ptscales,
		 CSCSectorReceiverLUT* srLUTs_[5][2]);
  
  void fillRPC(const edm::Handle< vector<L1TMuon::TriggerPrimitive> > rpc_hits);
  
  int check_track_pt(vector<vector<int>> cschits,
		      float track_pt);

  
  void three_track_pt(vector<vector<int>> cschits,
		      int nTrk,
		      int trkPt,
		      int rpc_cluster_phiBit,
		      double rpc_cluster_phi,
		      double rpc_cluster_eta,
		      const edm::Handle< vector<L1TMuon::TriggerPrimitive> > rpc_hits);
  
  void fillGenParticles(const edm::Handle<reco::GenParticleCollection> genps,
			edm::Handle<CSCSegmentCollection> cscSegments,
			edm::Handle< vector<L1TMuon::InternalTrack> > CSCTFtracks);

  void fillMuons(const edm::Handle<reco::MuonCollection> muons,
		 edm::Handle<CSCSegmentCollection> cscSegments,
		 edm::Handle< vector<L1TMuon::InternalTrack> > CSCTFtracks);

  int closest( vector<int> const& vec, int value);
  
  int isMuonMatched(__gnu_cxx::__normal_iterator<const reco::Muon*, std::vector<reco::Muon> >& muon,
		     edm::Handle<CSCSegmentCollection> cscSegments,
		     edm::Handle< vector<L1TMuon::InternalTrack> > CSCTFtracks,
		     int printLevel);
  
  
  
 //-------Define objects to be saved in TTree--------------------------------
 //--------------------------------------------------------------------------

  int Run, Event, Lumi;

  // Put generated muons in CSC tree
  float gen_pt;
  float gen_eta;
  float gen_phi;
  unsigned int gen_chg;
  
  // Reco Muon variables

  
  double  muon_pt[MAX_MUONS];
  double muon_eta[MAX_MUONS];
  double muon_phi[MAX_MUONS];
  float  muon_chg[MAX_MUONS];
  int  muon_match[MAX_MUONS]; // is muon matched to csctf track?  Gives -999 if no match, and track number if matched.

  int muonSize; // num muons in event

  // CSC objects ------------------------------------------------------------
  int SizeTrk; // # of tracks in event
  
  float PtTrk   [MAX_CSCTF_TRK];
  double EtaTrk [MAX_CSCTF_TRK];
  double PhiTrk [MAX_CSCTF_TRK];
  int EtaBitTrk [MAX_CSCTF_TRK];
  int PhiBitTrk [MAX_CSCTF_TRK];
  int PtBitTrk  [MAX_CSCTF_TRK];  
  int ModeTrk   [MAX_CSCTF_TRK];
  int ChgTrk    [MAX_CSCTF_TRK];

  // These store track pT info from PtAddtess.h when calculating two hit tracks and adding one rpc hit

  int PtTrk_reco_front     [MAX_CSCTF_TRK][MAX_RPC_LCTS_EVT]; // pT for front layer
  int PtTrk_reco_rear      [MAX_CSCTF_TRK][MAX_RPC_LCTS_EVT]; // pT for rear layer
  int PtTrk_reco_best_two  [MAX_CSCTF_TRK][MAX_RPC_LCTS_EVT]; // reco pt for rpc into the track
  int rpcIsmatched         [MAX_CSCTF_TRK][MAX_RPC_LCTS_EVT];
  int trkIsMatched         [MAX_CSCTF_TRK]; // Is this track matched to an rpc
  int RpcMatchedID         [MAX_CSCTF_TRK]; // Which RPC hit is macthed to track
  int isGoodTrk            [MAX_CSCTF_TRK]; // Does track have good pT assignment
  
  // NEW decluster RPC phi pt assignment
  int PtTrk_cluster_front [MAX_CSCTF_TRK]; 
  int PtTrk_cluster_rear [MAX_CSCTF_TRK];
  
  // These store track pT when using three hit tracks and replacing one with an rpc hit
  // these have been phased out
  int PtTrk_reco_front_three [MAX_CSCTF_TRK]; // pT for front layer
  int PtTrk_reco_rear_three  [MAX_CSCTF_TRK]; // pT for rear layer
  int PtTrk_reco_best        [MAX_CSCTF_TRK]; // Used in place of front/rear.  ALready has best pt value stored
  
  // THis gives best pt value from PtAddress.h for each time an rpc hit replaced a csc hit in 3 hit tracks.
  // [track][csc lct][rpc lct]
  int PtTrk_reco_rpc_csc [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int csc_rpc_match      [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int csc_rpc_id         [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];

  // New three track branches for cluster tests.  Replace a csc hit with an rpc cluster hit.
  int PtTrk_threeHit_rpc_cluster [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];


  // Info about LCTs forming the track
  // Stored in 2D array [track x] [LCT y] - each event has x tracks containing y Lcts
  int NumLctsTrk[MAX_CSCTF_TRK];  // How many Lcts made up track

  int rpc_NumLctsTrk;   // How many rpc lcts in event

  
  double trLctglobalPhi [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK]; 
  double trLctglobalEta [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctPhiBit       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEtaBit       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEtaBitMatt   [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctSector       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctSubSector    [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctStation      [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctChamber      [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctEndcap       [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctBx           [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctCSCId        [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  int trLctRing         [MAX_CSCTF_TRK][MAX_LCTS_PER_TRK];
  
  // RPC objects ------------------------------------------------------------
  
  double rpc_gblEta      [MAX_RPC_LCTS_EVT];
  double rpc_gblPhi      [MAX_RPC_LCTS_EVT];
  //unsigned int rpc_strip [MAX_RPC_LCTS_EVT];
  unsigned rpc_Layer [MAX_RPC_LCTS_EVT];
  int rpc_bx             [MAX_RPC_LCTS_EVT];
  int rpc_Station        [MAX_RPC_LCTS_EVT];
  int rpc_Sector         [MAX_RPC_LCTS_EVT];
  int rpc_Phibit         [MAX_RPC_LCTS_EVT]; 
  int rpc_Etabit         [MAX_RPC_LCTS_EVT];

  // New cluster objects.  Per event each station has one cluster value
  double rpc_stat1_cluster_phi;
  double rpc_stat2_cluster_phi;
  double rpc_stat3_cluster_phi;
  double rpc_stat4_cluster_phi;

  double rpc_stat1_cluster_eta;
  double rpc_stat2_cluster_eta;
  double rpc_stat3_cluster_eta;
  double rpc_stat4_cluster_eta;

  int rpc_stat1_cluster_phiBit;
  int rpc_stat2_cluster_phiBit;
  int rpc_stat3_cluster_phiBit;
  int rpc_stat4_cluster_phiBit;

  // Method to set the tree branches to our delcared variables
  void set_tree() {

    CSCtree = new TTree("CSCtree", "CSCtree");

    CSCtree -> Branch("Event" ,  &Event , "Event/I");
    CSCtree -> Branch("Run"   ,  &Run   , "Run/I");
    CSCtree -> Branch("Lumi"  ,  &Lumi  , "Lumi/I");

    // Fill gen particles
    CSCtree -> Branch("gen_pt"    ,  &gen_pt   , "gen_pt/F");
    CSCtree -> Branch("gen_eta"   ,  &gen_eta  , "gen_eta/F");
    CSCtree -> Branch("gen_phi"   ,  &gen_phi  , "gen_phi/F");
    CSCtree -> Branch("gen_chg"   ,  &gen_chg  , "gen_chg/I");

    // Fill Reco variables
    CSCtree -> Branch("muonSize"   ,  &muonSize , "muonSize/I");
    CSCtree -> Branch("muon_pt"    ,  muon_pt   , "muon_pt[100]/D");
    CSCtree -> Branch("muon_eta"   ,  muon_eta  , "muon_eta[100]/D");
    CSCtree -> Branch("muon_phi"   ,  muon_phi  , "muon_phi[100]/D");
    CSCtree -> Branch("muon_chg"   ,  muon_chg  , "muon_chg[100]/F");
    CSCtree -> Branch("muon_match" ,  muon_match, "muon_match[100]/I");

    // CSCTF track variables
    CSCtree -> Branch("SizeTrk"   ,  &SizeTrk  , "SizeTrk/I");
    CSCtree -> Branch("EtaTrk"    ,  EtaTrk    , "EtaTrk[SizeTrk]/D");
    CSCtree -> Branch("ChgTrk"    ,  ChgTrk    , "ChgTrk[SizeTrk]/I");
    CSCtree -> Branch("PhiTrk"    ,  PhiTrk    , "PhiTrk[SizeTrk]/D");
    CSCtree -> Branch("PtTrk"     ,  PtTrk     , "PtTrk[SizeTrk]/F");
    CSCtree -> Branch("PtBitTrk"  ,  PtBitTrk  , "PtBitTrk[SizeTrk]/I");
    CSCtree -> Branch("EtaBitTrk" ,  EtaBitTrk , "EtaBitTrk[SizeTrk]/I");
    CSCtree -> Branch("PhiBitTrk" ,  PhiBitTrk , "PhiBitTrk[SizeTrk]/I");
    CSCtree -> Branch("ModeTrk"   ,  ModeTrk   , "ModeTrk[SizeTrk]/I");
    
    CSCtree -> Branch("PtTrk_reco_front"      , PtTrk_reco_front     , "PtTrk_reco_front[SizeTrk][36]/I");
    CSCtree -> Branch("PtTrk_reco_rear"       , PtTrk_reco_rear      , "PtTrk_reco_rear[SizeTrk][36]/I");
    CSCtree -> Branch("PtTrk_reco_best_two"   , PtTrk_reco_best_two  , "PtTrk_reco_best_two[SizeTrk][36]/I");
    CSCtree -> Branch("rpcIsmatched"          , rpcIsmatched         , "rpcIsmatched[SizeTrk][36]/I");
    
    CSCtree -> Branch("PtTrk_reco_best"  , PtTrk_reco_best  , "PtTrk_reco_best[SizeTrk]/I");
    CSCtree -> Branch("trkIsMatched"     , trkIsMatched     , "trkIsMatched[SizeTrk]/I");
    CSCtree -> Branch("RpcMatchedID"     , RpcMatchedID     , "RpcMatchedID[SizeTrk]/I");
    CSCtree -> Branch("isGoodTrk"        , isGoodTrk        , "isGoodTrk[SizeTrk]/I");
    
    CSCtree -> Branch("PtTrk_reco_front_three" , PtTrk_reco_front_three , "PtTrk_reco_front_three[SizeTrk]/I");
    CSCtree -> Branch("PtTrk_reco_rear_three"  , PtTrk_reco_rear_three  , "PtTrk_reco_rear_three[SizeTrk]/I");

    CSCtree -> Branch("PtTrk_threeHit_rpc_cluster"  , PtTrk_threeHit_rpc_cluster, "PtTrk_threeHit_rpc_cluster[SizeTrk][4]/I");

    // These are the variables LCTs that belong to the CSCTF track
    CSCtree -> Branch("NumLctsTrk"     , NumLctsTrk   , "NumLctsTrk[SizeTrk]/I");
    CSCtree -> Branch("trLctglobalPhi" , trLctglobalPhi , "trLctglobalPhi[SizeTrk][4]/D");
    CSCtree -> Branch("trLctglobalEta" , trLctglobalEta , "trLctglobalEta[SizeTrk][4]/D");
    CSCtree -> Branch("trLctPhiBit"    , trLctPhiBit    , "trLctPhiBit[SizeTrk][4]/I");
    CSCtree -> Branch("trLctEtaBit"    , trLctEtaBit    , "trLctEtaBit[SizeTrk][4]/I");
    CSCtree -> Branch("trLctSector"    , trLctSector    , "trLctSector[SizeTrk][4]/I");
    CSCtree -> Branch("trLctSubSector" , trLctSubSector , "trLctSubSector[SizeTrk][4]/I");
    CSCtree -> Branch("trLctStation"   , trLctStation   , "trLctStation[SizeTrk][4]/I");
    CSCtree -> Branch("trLctChamber"   , trLctChamber   , "trLctChamber[SizeTrk][4]/I");
    CSCtree -> Branch("trLctEndcap"    , trLctEndcap    , "trLctEndcap[SizeTrk][4]/I");
    CSCtree -> Branch("trLctBx"        , trLctBx        , "trLctBx[SizeTrk][4]/I");
    CSCtree -> Branch("trLctCSCId"     , trLctCSCId     , "trLctCSCId[SizeTrk][4]/I");    
    CSCtree -> Branch("trLctRing"      , trLctRing      , "trLctRing[SizeTrk][4]/I");    
    
    // RPC Trigger Primitive Branches(all hits in event record)
    CSCtree -> Branch("rpc_NumLctsTrk", &rpc_NumLctsTrk , "rpc_NumLctsTrk/I");
    
    CSCtree -> Branch("PtTrk_cluster_front", &PtTrk_cluster_front , "PtTrk_cluster_front[SizeTrk]/I");
    CSCtree -> Branch("PtTrk_cluster_rear", &PtTrk_cluster_rear , "PtTrk_cluster_rear[SizeTrk]/I");
    
    CSCtree -> Branch("rpc_stat1_cluster_phi", &rpc_stat1_cluster_phi , "rpc_stat1_cluster_phi/D");
    CSCtree -> Branch("rpc_stat2_cluster_phi", &rpc_stat2_cluster_phi , "rpc_stat2_cluster_phi/D");
    CSCtree -> Branch("rpc_stat3_cluster_phi", &rpc_stat3_cluster_phi , "rpc_stat3_cluster_phi/D");
    CSCtree -> Branch("rpc_stat4_cluster_phi", &rpc_stat4_cluster_phi , "rpc_stat4_cluster_phi/D");

    CSCtree -> Branch("rpc_stat1_cluster_eta", &rpc_stat1_cluster_eta , "rpc_stat1_cluster_eta/D");
    CSCtree -> Branch("rpc_stat2_cluster_eta", &rpc_stat2_cluster_eta , "rpc_stat2_cluster_eta/D");
    CSCtree -> Branch("rpc_stat3_cluster_eta", &rpc_stat3_cluster_eta , "rpc_stat3_cluster_eta/D");
    CSCtree -> Branch("rpc_stat4_cluster_eta", &rpc_stat4_cluster_eta , "rpc_stat4_cluster_eta/D");
    
    CSCtree -> Branch("rpc_stat1_cluster_phiBit", &rpc_stat1_cluster_phiBit , "rpc_stat1_cluster_phiBit/I");
    CSCtree -> Branch("rpc_stat2_cluster_phiBit", &rpc_stat2_cluster_phiBit , "rpc_stat2_cluster_phiBit/I");
    CSCtree -> Branch("rpc_stat3_cluster_phiBit", &rpc_stat3_cluster_phiBit , "rpc_stat3_cluster_phiBit/I");
    CSCtree -> Branch("rpc_stat4_cluster_phiBit", &rpc_stat4_cluster_phiBit , "rpc_stat4_cluster_phiBit/I");
    

    CSCtree -> Branch("csc_rpc_match"       , csc_rpc_match       , "csc_rpc_match[SizeTrk][4]/I");
    CSCtree -> Branch("PtTrk_reco_rpc_csc"  , PtTrk_reco_rpc_csc  , "PtTrk_reco_rpc_csc[SizeTrk][4]/I");
    CSCtree -> Branch("csc_rpc_id"          , csc_rpc_id          , "csc_rpc_id[SizeTrk][4]/I");

    CSCtree -> Branch("rpc_gblEta"    ,  rpc_gblEta     , "rpc_gblEta[36]/D");
    CSCtree -> Branch("rpc_gblPhi"    ,  rpc_gblPhi     , "rpc_gblPhi[36]/D");
    CSCtree -> Branch("rpc_Layer"     ,  rpc_Layer      , "rpc_Layer[36]/I");
    //CSCtree -> Branch("rpc_strip"     ,  rpc_strip      , "rpc_strip[36]/I");
    CSCtree -> Branch("rpc_bx"        ,  rpc_bx         , "rpc_bx[36]/I");
    CSCtree -> Branch("rpc_Station"   ,  rpc_Station    , "rpc_Station[36]/I");
    CSCtree -> Branch("rpc_Sector"    ,  rpc_Sector     , "rpc_Sector[36]/I");
    CSCtree -> Branch("rpc_Phibit"    ,  rpc_Phibit     , "rpc_Phibit[36]/I");
    
    
    return;
  }
  
  
  void fill_CSCtree() { 
    CSCtree -> Fill();}


  // ------- Function to init the TTree every event; needed in conjuction with del_tree()
  // ------------------------------------------------------------------------------------
  void init_tree() {
    

  }


  void del_tree() {

  }
  

  // -------- Function to set phi/eta CSCTF LUTs
  // --------------------------------------------
  void set_etaLUT() {

    bzero(srLUTs_ , sizeof(srLUTs_));
    int sector=1;    // assume SR LUTs are all same for every sector
    bool TMB07=true; // specific TMB firmware
    edm::ParameterSet srLUTset;
    srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
    srLUTset.addUntrackedParameter<bool>("Binary",   false);
    srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
 
    // positive endcap
    int endcap = 1;
    for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
      if(station==1)
        for(int subSector=0; subSector<2 && fpga<5; subSector++)
          srLUTs_[fpga++][1] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
                                                        station, srLUTset, TMB07);
      else
        srLUTs_[fpga++][1]   = new CSCSectorReceiverLUT(endcap,  sector,   0,
                                                        station, srLUTset, TMB07);
    }
    // negative endcap
    endcap = 2;
    for(int station=1,fpga=0; station<=4 && fpga<5; station++) {
      if(station==1)
        for(int subSector=0; subSector<2 && fpga<5; subSector++)
          srLUTs_[fpga++][0] = new CSCSectorReceiverLUT(endcap,sector,subSector+1,
                                                        station, srLUTset, TMB07);
      else
        srLUTs_[fpga++][0]   = new CSCSectorReceiverLUT(endcap,  sector,   0,
                                                        station, srLUTset, TMB07);
    }

    
    return;

  }
  
  // Initialize the CSCTF pt LUT
  void pt_lut_init() {        
    return;
  }


  // pT scaling used with Matt's PtAddress.h
  int scaling(double a)//original pt scale
  {
    if (a >= 0 && a < 1.5){ a = 0.0;} else if (a >= 1.5 && a < 2.0){ a = 1.5;} else if (a >= 2.0 && a < 2.5){ a = 2.5;} else if (a >= 2.5 && a < 3.0){ a = 3.0;}
    else if (a >= 3.0 && a < 3.5){ a = 3.0;} else if (a >= 3.5 && a < 4.0){ a = 3.5;} else if(a >= 4.0 && a < 4.5){a = 4.0;} else if(a >= 4.5 && a < 5.0){a = 4.5;}
    else if (a >= 5.0 && a < 6.0){ a = 5.0;} else if (a >= 6.0 && a < 7.0){ a = 6.0;} else if(a >= 7.0 && a < 8.0){a = 7.0;} else if(a >= 8.0 && a < 10.0){a = 8.0;}
    else if (a >= 10.0 && a < 12.0){ a = 10.0;} else if (a >= 12.0 && a < 14.0){ a = 12.0;} else if(a >= 14.0 && a < 16.0){a = 14.0;} else if(a >= 16.0 && a < 18.0){a = 16.0\
	;}
    else if (a >= 18.0 && a < 20.0){ a = 18.0;} else if (a >= 20.0 && a < 25.0){ a = 20.0;} else if(a >= 25.0 && a < 30.0){a = 25.0;} else if(a >= 30.0 && a < 35.0){a = 30.0\
	;}
    else if (a >= 35.0 && a < 40.0){ a = 35.0;} else if (a >= 40.0 && a < 45.0){ a = 40.0;} else if(a >= 45.0 && a < 50.0){a = 45.0;} else if(a >= 50.0 && a < 60.0){a = 50.0\
	;}
    else if (a >= 60.0 && a < 70.0){ a = 60.0;} else if (a >= 70.0 && a < 80.0){ a = 70.0;} else if(a >= 80.0 && a < 90.0){a = 80.0;} else if(a >= 90.0 && a < 100.0){a = 90.0\
	;}
    else if (a >= 100.0 && a < 120.0){ a = 100.0;} else if (a >= 120.0 && a < 140.0){ a = 120.0;} else if(a >= 140.0 && a < 1E6){a = 140.0;}
    return a;
  }

  
  // -------------------------------------------------------------------------------------------
  // RPC phi: calculating rpc global phi in CSC bit form
  /*
  int CalcIntegerPhi(float GblPhi,int sector,int station,int ring){
    
    
    float slope[3][4][6] = { { {3783.79, 3784.91, 3784.06, 3784.56, 3784.05, 3784.86}, {3783.79, 3784.91, 3784.06, 3784.56, 3784.05, 3784.86}, \
			     {3784.54,3784.61,3784.54,3785.22,3784.56,3784.51}, {3783.79,3784.91,3784.06,3784.56,3784.05,3784.86} },
			    { {0,0,0,0,0,0}, {3784.07, 3784.65, 3784.17, 3784.59, 3784.17, 3785.08}, {0,0,0,0,0,0}, {0,0,0,0,0,0} },
			    { {0,0,0,0,0,0}, {3784.36,3784.68,3784.09,3784.61,3784.14,3785.02}, {0,0,0,0,0,0}, {0,0,0,0,0,0} } };
    


    float offset[3][4][6] = { { {-924.362, -4888.94, -8849.91, 10964.5, 7000.71, 3038.25}, {-924.362, -4888.94, -8849.91, 10964.5, 7000.71, 3038.25}, \
				{-925.02, -4888.23, -8851.1, 10966.2, 7001.31, 3038.31}, {-925.02, -4888.23, -8851.1, 10966.2, 7001.31, 3038.31} },
			      { {0,0,0,0,0,0}, {-924.479, -4888.28, -8850.11, 10964.6, 7001.03, 3038.3}, {0,0,0,0,0,0}, {0,0,0,0,0,0} },
			      { {0,0,0,0,0,0}, {-924.573, -4888.41, -8849.93, 10964.6, 7000.99, 3038.31}, {0,0,0,0,0,0}, {0,0,0,0,0,0} } };
    int phi_out = -1;
    
    if (sector == 3 and GblPhi < 0) GblPhi += 2*3.14159;
    
    phi_out = (GblPhi*slope[ring-1][station-1][sector-1]) + offset[ring-1][station-1][sector-1];
    
    if (phi_out < 0) phi_out += 4096;
    
    if (phi_out > 4096) phi_out -= 4096;
    
    return phi_out;
  }
  */

  // -------------------------------------------------------------------------------------------
  

  int CalcIntegerPhi(float GblPhi, int sector) {
    
    float phiBit = -999;

    if (sector == 1)
      phiBit = (GblPhi - 0.243) / 1.0835;

    if (sector == 2)
      phiBit = (GblPhi - 1.2914) / 1.0835;

    if (sector == 3) {
      if (GblPhi > 0)
	phiBit = (GblPhi - 2.338) / 1.0835;
      if (GblPhi < 0) {
	float sector_distance = abs(GblPhi + 3.1416) + (3.1416 - 2.338);
        phiBit = sector_distance / 1.0835;
      }
    }

    if (sector == 4)
      phiBit = (GblPhi + 2.898) / 1.0835;

    if (sector == 5)
      phiBit = (GblPhi + 1.8507) / 1.0835;
    
    if (sector == 6) {
      if (GblPhi < 0)
	phiBit = (GblPhi + 0.803) / 1.0835;
      if (GblPhi > 0) {
	float sector_distance = GblPhi + 0.803;
	phiBit = sector_distance / 1.0835;
      }	  
    }
    
    phiBit = phiBit*4096;
    
    //cout << "PhiBit = " << phiBit << endl;

    return static_cast<int>(phiBit);
      
  }



}; // end class declaration

#endif



