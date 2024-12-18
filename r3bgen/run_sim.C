void run_sim()
{
    TString transport = "TGeant4";

    TString outFile = "sim.root";
    TString parFile = "par.root";
    Bool_t magnet = kFALSE;
    Float_t fieldScale = -0.6;

    TString generator1 = "box";
    TString generator2 = "ascii";
    TString generator3 = "r3b";
    TString generator4 = "qfs";
    TString generator = generator4;
    TString inputFile = "";

    Int_t nEvents = 100000;
    Bool_t storeTrajectories = kTRUE;
    Int_t randomSeed = 335566; // 0 for time-dependent random numbers

    // Target type
    TString target1 = "LeadTarget";
    TString target2 = "Para";
    TString target3 = "Para45";
    TString target4 = "LiH";
    TString targetType = target4;

    // ------------------------------------------------------------------------
    // Stable part ------------------------------------------------------------

    TString dir = getenv("VMCWORKDIR");

    // ----    Debug option   -------------------------------------------------
    gDebug = 0;

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();

    // -----   Create simulation run   ----------------------------------------
    FairRunSim* run = new FairRunSim();
    run->SetName(transport);            // Transport engine
    run->SetOutputFile(outFile.Data()); // Output file
    FairRuntimeDb* rtdb = run->GetRuntimeDb();

    // -----   Create media   -------------------------------------------------
    run->SetMaterials("media_r3b.geo"); // Materials

    // -----   Create R3B geometry --------------------------------------------
    // R3B Cave definition
    FairModule* cave = new R3BCave("CAVE");
    cave->SetGeometryFileName("r3b_cave_vacuum.geo");
    run->AddModule(cave);

    // To skip the detector comment out the line with: run->AddModule(...

    // Target
    run->AddModule(new R3BTarget(targetType, "target_" + targetType + ".geo.root"));

    // GLAD
    //run->AddModule(new R3BGladMagnet("glad_v2023.1.geo.root")); // GLAD should not be moved or rotated

    // PSP

    // CALIFA
    /*R3BCalifa* califa = new R3BCalifa("califa_full.geo.root");
      califa->SelectGeometryVersion(2020);
      run->AddModule(califa);*/


    
    // run->AddModule(new R3BNeuland("neuland_test.geo.root", { 0., 0., 1400. + 12 * 5. }));

    // -----   Create R3B  magnetic field ----------------------------------------
    // NB: <D.B>
    // If the Global Position of the Magnet is changed
    // the Field Map has to be transformed accordingly
    R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap");
    magField->SetScale(fieldScale);

    if (magnet == kTRUE)
    {
	    run->SetField(magField);
    }
    else
    {
	    run->SetField(NULL);
    }

    // -----   Create PrimaryGenerator   --------------------------------------
    // 1 - Create the Main API class for the Generator
    FairPrimaryGenerator* primGen = new FairPrimaryGenerator();

    if (generator.CompareTo("box") == 0)
    {
	    // 2- Define the BOX generator
	    Int_t pdgId = 2212;     // proton beam <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<?????????
	    Double32_t theta1 = 0.; // polar angle distribution
	    Double32_t theta2 = 2.;
	    Double32_t momentum = 1.5;
	    FairBoxGenerator* boxGen = new FairBoxGenerator(pdgId, 3); //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<?????
	    boxGen->SetThetaRange(theta1, theta2);
	    boxGen->SetPRange(momentum, momentum * 1.2);
	    boxGen->SetPhiRange(0, 360);
	    boxGen->SetXYZ(0.0, 0.0, 0);
	    primGen->AddGenerator(boxGen);

	    // 128-Sn fragment
	    R3BIonGenerator* ionGen = new R3BIonGenerator(50, 128, 50, 10, 1.3);
	    ionGen->Beam.SetVertexDistribution(
			    R3BDistribution3D::Prism(R3BDistribution2D::Circle({ 0., 0. }, 0.1), R3BDistribution1D::Delta(-300)));
	    primGen->AddGenerator(ionGen);

	    // neutrons
	    FairBoxGenerator* boxGen_n = new FairBoxGenerator(2112, 3);
	    boxGen_n->SetThetaRange(theta1, theta2);
	    boxGen_n->SetPRange(momentum, momentum * 1.2);
	    boxGen_n->SetPhiRange(0, 360);
	    boxGen_n->SetXYZ(0.0, 0.0, 0);
	    primGen->AddGenerator(boxGen_n);
    }
    
    if (generator.CompareTo("qfs") == 0)
    {
       auto qfsGen = new R3BQFSGenerator();
       //qfsGen->SetProjectile(12, 6, 1250); //A, Z, energy per nucleon
       qfsGen->SetProjectile(1, 1, 10000000);
       qfsGen->SetHeavyFragment(2,1); //A, Z, energy per nucleon -> 12C(p,2p)11B
       qfsGen->SetInverse(true);                                                                                       
       qfsGen->SetIsotropic(false);
	                                                                                               
       qfsGen->SetAddGamma(0);
       qfsGen->SetMomDistrib(0);                                                                     
      // qfsGen->SetExcitation(0);//<<<<<<<<<<<<<<<<<<<<<<<<<Setting for excitation energy, equals to gamma energy                                                              
       //qfsGen->SetLightNucleus(938.272, 938.272,false); //The last boolean is for settting (p,pn)=true and (p,2p)= false  beam,target                                                                                
       qfsGen->SetLightNucleus(938.272, 938.272,false);
       qfsGen->Print();    
       primGen->AddGenerator(qfsGen);
    }
    run->SetGenerator(primGen);

    run->SetStoreTraj(storeTrajectories);

    FairLogger::GetLogger()->SetLogVerbosityLevel("LOW");
    FairLogger::GetLogger()->SetLogScreenLevel("info");
    //FairLogger::GetLogger()->SetLogScreenLevel("debug");

    // -----   Initialize simulation run   ------------------------------------
    TRandom3 random(randomSeed);
    gRandom = &random;
    run->Init();

    // ------  Increase nb of step for CALO
    Int_t nSteps = -15000;
    TVirtualMC::GetMC()->SetMaxNStep(nSteps);

    // -----   Runtime database   ---------------------------------------------
    R3BFieldPar* fieldPar = (R3BFieldPar*)rtdb->getContainer("R3BFieldPar");
    if (NULL != magField)
    {
        fieldPar->SetParameters(magField);
        fieldPar->setChanged();
    }
    Bool_t kParameterMerged = kTRUE;
    FairParRootFileIo* parOut = new FairParRootFileIo(kParameterMerged);
    parOut->open(parFile.Data());
    rtdb->setOutput(parOut);
    rtdb->saveOutput();
    rtdb->print();

    // -----   Start run   ----------------------------------------------------
    if (nEvents > 0)
    {
        run->Run(nEvents);
    }

    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Output file is " << outFile << endl;
    cout << "Parameter file is " << parFile << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;

    cout << " Test passed" << endl;
    cout << " All ok " << endl;

}
