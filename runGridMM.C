// Grid run task for Diffractive Cross Sections at \sqrt{s} = 13 TeV
//
//
// Submit with    runmode = "full"
// Test with      runmode = "test"
// Terminate with runmode = "terminate"
//
// mikael.mieskolainen@cern.ch, 2018
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// Skeleton of this class is based on ALICE collaboration common grid code and macros

TString gGridDataDir;
TString gRunPrefix;
TString gGridWorkingDir;
TString gDataPattern;


void runGridMM(TString runMode, TString period="LHC15f", TString pass="pass2") {

  // Addpaths for local root
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$/ROOTSYS/include -I./");

  // Addpaths
// gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ROOTSYS")));
// gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
// gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_PHYSICS")));
// gROOT->ProcessLine(".include ./");

  // Load own files
  gROOT->LoadMacro("AliAnalysisTaskDiffCrossSectionsMM.cxx+");
  gROOT->LoadMacro("AddAnalysisTaskDiffCrossSectionsMM.cxx");

  //---------------------------------------------------------------------------
  // 0. Initialize variables

  TString mcTypes[5] = {
    "",         // 0
    "Pythia6",  // 1
    "Pythia8",  // 2
    "PHOJET",   // 3
    "EPOS"      // 4
  };

  // The MC numbers below match the numbers above
  Int_t indexMC_2015 = (1*period.Contains("LHC15g3c2") +  // pythia6
                        2*period.Contains("LHC15g3a2") +  // pythia8
                        2*period.Contains("LHC15g3a3"));  // pythia8

  Int_t indexMC_2016 = (1*period.Contains("LHC16a2a2") +  // pythia6
                        2*period.Contains("LHC16a2b2") +  // pythia8
                        3*period.Contains("LHC16a2c2") +  // phojet
                        4*period.Contains("LHC16a2d2"));  // epos

  Int_t indexMC_2017 = (1*period.Contains("LHC17h7a") +   // pythia6  
                        3*period.Contains("LHC17h7b"));   // phojet
  
  
  Bool_t isMC     = kFALSE;
  TString trigSel = "";
  TString mcType  = "";

  // Data 2015
  if (period.Contains("LHC15f") || period.Contains("LHC15h")) {
    isMC            = kFALSE;
    gRunPrefix      = "000";
    gGridDataDir    = "/alice/data/2015/"+period;
    gDataPattern    = "/"+pass+"/*/AliESDs.root";
    gGridWorkingDir = period+"_"+pass+"_DiffCS";

  // Data 2017
  } else if (period.Contains("LHC17j")) {

    isMC            = kFALSE;
    gRunPrefix      = "000";
    gGridDataDir    = "/alice/data/2017/"+period;
    gDataPattern    = "/"+pass+"/*/AliESDs.root";
    gGridWorkingDir = period+"_"+pass+"_DiffCS";

  // 2015 produced MC
  } else if (indexMC_2015) {
    isMC            = kTRUE;
    gRunPrefix      = "";
    gGridDataDir    = "/alice/sim/2015/"+period;
    gDataPattern    = "/*/AliESDs.root";
    gGridWorkingDir = period+"_DiffCS";
    mcType          = mcTypes[indexMC_2015];

  // 2016 produced MC
  } else if (indexMC_2016) {
    isMC            = kTRUE;
    gRunPrefix      = "";
    gGridDataDir    = "/alice/sim/2016/"+period;
    gDataPattern    = "/*/AliESDs.root";
    gGridWorkingDir = period+"_DiffCS";
    mcType          = mcTypes[indexMC_2016];
  
  // 2017 produced MC
  } else if (indexMC_2017) {
    isMC            = kTRUE;
    gRunPrefix      = "";
    gGridDataDir    = "/alice/sim/2017/"+period;
    gDataPattern    = "/*/AliESDs.root";
    gGridWorkingDir = period+"_DiffCS";
    mcType          = mcTypes[indexMC_2017];
  }


  //---------------------------------------------------------------------------
  // 1. Create Analysis Manager
  AliAnalysisManager* mgr = new AliAnalysisManager("TestManager");
  mgr->SetCacheSize(100*1000*1000);


  //---------------------------------------------------------------------------
  // Create Analysis object (wrapper function in this file)

  AliAnalysisAlien* alienHandler = CreateAlienHandler(runMode);
  if (NULL == alienHandler) {
    Error("Error:: alienHandler == NULL");
    return;
  }

  // ---------------------------------------------------------------------
  // Setup paths for grid

  // Declare the analysis source files names separated by blancs. To be
  // compiled runtime using ACLiC on the worker nodes.

  alienHandler->SetAnalysisSource("AliAnalysisTaskDiffCrossSectionsMM.cxx");

  // Documentation says that (.cxx/.h) needs to be added here also
  alienHandler->SetAdditionalLibs("AliAnalysisTaskDiffCrossSectionsMM.h AliAnalysisTaskDiffCrossSectionsMM.cxx");


  //---------------------------------------------------------------------------
  // Read run numbers from text files and add the run numbers

  TTree t;
  t.ReadFile(period+".txt", "rn/I");

  Int_t rn = 0;
  t.SetBranchAddress("rn", &rn);
  for (Int_t i = 0, n(t.GetEntries()); i < n; ++i) {
    t.GetEntry(i);
    alienHandler->AddRunNumber(rn);
  }
  mgr->SetGridHandler(alienHandler);


  //---------------------------------------------------------------------------
  // 2. Create ESD input handler

  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetNeedField();
  mgr->SetInputEventHandler(esdH);


  //---------------------------------------------------------------------------
  // Create "Physics Selection" task
  // (we are not going to use it, but use our own selection)

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
  AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
  physSel->SetUseBXNumbers(kFALSE);

/* // AOD
  
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC);
  if(!physSelTask) { printf("No physSelTask"); return; }
  AliOADBPhysicsSelection* ps = new AliOADBPhysicsSelection("custom");
  ps->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+CINT10-B-NOPF-ALLNOTRD","B",0);
  ps->SetHardwareTrigger(0,"V0A || V0C || ADA || ADC");
  ps->SetOfflineTrigger(0,"(V0A || V0C || ADA || ADC) && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");
  ps->AddCollisionTriggerClass(AliVEvent::kUserDefined,"+C0SMB-B-NOPF-ALLNOTRD","B",1);
  ps->SetHardwareTrigger(1,"SPDGFO >= 1");
  ps->SetOfflineTrigger(1,"SPDGFO >= 1 && !V0ABG && !V0CBG && !ADABG && !ADCBG && !TPCLaserWarmUp");
  physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(ps,0);
*/
  
  //---------------------------------------------------------------------------
  // Add analysis tasks. Find out the spesific trigger names from the Alice e-logbook
  //
  //
  // AddAnalysisTaskDiff...() spans new AnalysisManagers.
  // To be verified if this is a 100% valid procedure.

  // Monte Carlo
  if (isMC) {

    AddAnalysisTaskDiffCrossSectionsMM(isMC, mcType, trigSel);

  // Data
  } else {

    if (period.Contains("LHC15f")) {
      TString trigSelData[4] = {
	       "C0SMB-B-NOPF-ALLNOTRD|CINT10-B-NOPF-ALLNOTRD",
	       "C0SMB-A-NOPF-ALLNOTRD|CINT10-A-NOPF-ALLNOTRD",
	       "C0SMB-C-NOPF-ALLNOTRD|CINT10-C-NOPF-ALLNOTRD",
	       "C0SMB-E-NOPF-ALLNOTRD|CINT10-E-NOPF-ALLNOTRD"};

      for (Int_t i = 0; i < 4; ++i) {
	      AddAnalysisTaskDiffCrossSectionsMM(isMC, mcType, trigSelData[i]);
      }

    } else if (period.Contains("LHC15h"))  {
      TString trigSelData[4] = {
	       "CINT11-B-NOPF-ALLNOTRD",
	       "CINT11-A-NOPF-ALLNOTRD",
	       "CINT11-C-NOPF-ALLNOTRD",
	       "CINT11-E-NOPF-ALLNOTRD"};

      for (Int_t i = 0; i < 4; ++i) {
	      AddAnalysisTaskDiffCrossSectionsMM(isMC, mcType, trigSelData[i]);
      }

    } else if (period.Contains("LHC17j"))  {
      TString trigSelData[4] = {
         "CINT11-B-NOPF-CENTNOTRD",
         "CINT11-A-NOPF-CENTNOTRD",
         "CINT11-C-NOPF-CENTNOTRD",
         "CINT11-E-NOPF-CENTNOTRD"};

      for (Int_t i = 0; i < 4; ++i) {
        AddAnalysisTaskDiffCrossSectionsMM(isMC, mcType, trigSelData[i]);
      }
    } else {

      printf("Error:: Unknown period: %s \n", period.Data());
      return;
    }

  }

  //---------------------------------------------------------------------------
  // 3. Start Analysis
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->StartAnalysis("grid");

}


// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
AliAnalysisGrid* CreateAlienHandler(TString runMode) {

  AliAnalysisAlien* plugin = new AliAnalysisAlien();

  plugin->SetOverwriteMode(kTRUE);                // Overwrite all previous generated files
  plugin->SetRunMode(runMode);                    // Run mode

  plugin->SetAPIVersion("V1.1x");
  
  //plugin->SetAliPhysicsVersion("vAN-20160720-1"); // *** AliPhysics version e.g. "vAN-20160720-1" ***
  plugin->SetAliPhysicsVersion("v5-09-13-01-1"); // 2017 productions

  plugin->SetGridDataDir(gGridDataDir);
  plugin->SetDataPattern(gDataPattern);
  plugin->SetRunPrefix(gRunPrefix);
  plugin->SetGridWorkingDir(gGridWorkingDir);
  plugin->SetGridOutputDir("output");             // Set output directory name
  plugin->SetOutputToRunNo(kTRUE);                // Output to run number folder

  plugin->SetAnalysisMacro("AnalysisMM.C");       // Set generated analysis macro name
  plugin->SetExecutable("AnalysisMM.sh");         // Set generated shell script name
  plugin->SetJDLName("AnalysisMM.jdl");           // Name of the .jdl to be generated

  plugin->SetCheckCopy(kFALSE);
  plugin->SetKeepLogs(kTRUE);
  plugin->SetMergeViaJDL();                       // Do merging on the grid

  plugin->SetDefaultOutputs();
  plugin->SetNtestFiles(1);                       // Number of files in "test" mode
  plugin->SetSplitMaxInputFileNumber(90);         // Maximum number of input files/subjob
  plugin->SetTTL(9*3600);                         // Job time to live (TTL)
  plugin->SetInputFormat("xml-single");           // Input format, default "xml-single"

  plugin->SetPrice(1);                            // Job price
  plugin->SetSplitMode("se");                     // Split mode, default "se", "file" mode splits too much!
  plugin->SetDropToShell(kFALSE);                 // ?

  // Addpaths for grid runtime task compilation (with ACLiC)
  plugin->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB");

  return plugin;

}
