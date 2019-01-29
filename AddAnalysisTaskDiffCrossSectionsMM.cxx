// Grid code steering macro
//
// mikael.mieskolainen@cern.ch, 2018
// Based on ALICE collaboration common grid code and macros

AliAnalysisTaskDiffCrossSectionsMM* AddAnalysisTaskDiffCrossSectionsMM(Bool_t isMC, TString mcType, TString triggerSelection) {
  
  //---------------------------------------------------------------------------
  // Create new Analysis Manager
  
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    mgr = new AliAnalysisManager("TestRunMM");
    
  //---------------------------------------------------------------------------
  // Create Analysis Monte Carlo handler 
  if (isMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  
  //---------------------------------------------------------------------------
  // Create the Analysis Task  
  AliAnalysisTaskDiffCrossSectionsMM* task = new AliAnalysisTaskDiffCrossSectionsMM;

  // Set the setup
  // task->SelectCollisionCandidates(AliVEvent::kMB); // This would force processing of only kMB tagged events
  task->SetIsMC(isMC);
  task->SetMCType(mcType);
  task->SetTriggerSelection(triggerSelection);
  
  
  //---------------------------------------------------------------------------
  // Create Data Container
  AliAnalysisDataContainer* output =
  mgr->CreateContainer(task->GetTreeName(), TTree::Class(), AliAnalysisManager::kOutputContainer,
                       TString(AliAnalysisManager::GetCommonFileName())+":"+task->GetResultsFileName());

  // Add the task
  mgr->AddTask(task);

  // Connect input and output
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

	
  return task;
}
