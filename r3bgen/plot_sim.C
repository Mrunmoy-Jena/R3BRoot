void plot_sim()
{
	const double PI = 3.14159265358979323846;
	bool IS_PN = false; //Set for (p,pn) = true
	TChain * ch = new TChain("evt");
	ch->Add("sim.root");
	TClonesArray * MCTrack = new TClonesArray("R3BMCTrack");
	ch->SetBranchAddress("MCTrack", &MCTrack);


	R3BMCTrack * proton[2];
	R3BMCTrack * candidate;

	TH2F * Theta_Correlation = new TH2F("Theta_Correlation","Theta Correlations; #theta_{1} [Deg]; #theta_{2} [Deg]",2000,0,90,2000,0,90);
	TH2F * OpangvsTheta = new TH2F("OpangvsTheta","Opening Angle vs Theta; #theta_{1} [Deg]; #theta_{op} [Deg]",2000,0,90,2000,0,120);
	TH2F * Theta_Correlation_califa_cut = new TH2F("Theta_Correlation_califa_cut","Theta Correlations_califa_cut; #theta_{1} [Deg]; #theta_{2} [Deg]",2000,0,90,2000,0,90);
	TVector3 proton_v3[2];

	TH2F * Phi_Correlation = new TH2F("Phi_Correlation","Phi Correlations; #phi_{1} [Deg]; #phi_{2} [Deg]",2000,-180,180,2000,-180,180);
	TH2F * Phi_Correlation_califa_cut = new TH2F("Phi_Correlation_califa_cut","Phi Correlations_califa_cut; #phi_{1} [Deg]; #phi_{2} [Deg]",2000,-180,180,2000,-180,180);
	TVector3 PhiProton_v3[2];
	
	TH2F * Kinetic_Energies_Correlation = new TH2F("Kinetic_Energies_Correlation", "Kinetic Energies Correlation [MeV]; #E_{kin,1}; #E_{kin,2} [MeV]", 200, 0, 2000, 200, 0, 2000);
	
	TH1F * phiproton1 = new TH1F("phiproton1","Phi 1; #phi_{1} [Deg];Events",2000,-180,180);
	TH1F * phiproton2 = new TH1F("phiproton2","Phi 2; #phi_{2} [Deg]; Events",2000,-180,180);
	TH1F * Theta_1 = new TH1F("Theta_1","Theta 1; #theta_{1} [Deg]; Events",500,-180,180);
	TH1F * Theta_2 = new TH1F("Theta_2","Theta 2; #theta_{2} [Deg]; Events",500,-180,180);
	TH1F * Relative_Theta = new TH1F("Opening_angle","Opening angle; #theta_{op} [Deg]; Events",2000,0,120);
	TH1F * Relative_Theta_califa_cut = new TH1F("Relative_Theta_califa_cut","Relative Theta_califa_cut; #theta_{op,cut} [Deg]; Theta 2[Deg]",2000,0,120);
	TH1F * Difference_Phi = new TH1F("Difference_Phi", "Difference between Phi(s); #Delta #phi [Deg]; Events",2000,0,360);
	TH1F * Distribution_Corrected_Gamma = new TH1F("Distribution_Corrected_Gamma", "Distribution of Corrected Gamma Energy; Energy [MeV]; Events", 100, 0, 30);
	TH1F * Gamma_Energy_Distribution = new TH1F("Gamma_Energy_Distribution", "Gamma Energy Distribution; Energy [MeV]; Events", 200, 0, 30);
	
//	TH1F * Angular_Energy_Distribution = new TH1F("Angular_Energy_Distribution", "Angular-Energy Distribution of Particle 1; Energy [MeV]; Theta [Deg]", 200, 0,180);	
//	TH1F * Angular_Energy_Distribution_2 = new TH1F("Angular_Energy_Distribution", "Angular-Energy Distribution of Particle 2 ; Energy [MeV]; Theta [Deg]", 200, 0,180);
	//Angular Distribution of Gamma Rays in Lab frme
	TH2F * Gamma_Energy_and_Theta = new TH2F("Gamma_Energy_and_Theta","Gamma Energy and Theta; #theta_{gamma} [Deg];Events",200,0,180,200,0,30);
	//eventloop
	

	TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
	
	TLorentzVector fragmentVector;
	TVector3 fragmentVector3;
	TVector3 boostVector;
	Double_t betaFragment;


	TRandom3 randomGen;
	for(int ev=0; ev < ch->GetEntries(); ev++)
	{
		ch->GetEntry(ev);

		proton[0] = nullptr;
		proton[1] = nullptr;

		for(int i=0; i<MCTrack->GetEntriesFast(); i++)
		{
			candidate = (R3BMCTrack*) MCTrack->At(i);

			if(candidate->GetMotherId()!=-1) continue ;
			//Fragment 2212 = proton $$ 2112 = neutron
			if(candidate->GetPdgCode()!=2112 && candidate->GetPdgCode() != 2212  && candidate->GetPdgCode()!=22 )
			{
				Double_t fragEnergy = candidate->GetEnergy();
				Double_t fragpx = candidate ->GetPx();
				Double_t fragpy = candidate ->GetPy();
				Double_t fragpz = candidate ->GetPz();
				fragmentVector.SetPxPyPzE(fragpx, fragpy, fragpz, fragEnergy);
				fragmentVector3.SetXYZ(fragpx, fragpy, fragpz);
				boostVector = fragmentVector.BoostVector();
			//	fragmentVector3.SetXYZ(boostVector.Px(), boostVector.Py(), boostVector.Pz()); //In lab frame
				betaFragment = boostVector.Mag();
				continue;
			}	
			//Gamma
			if(candidate->GetPdgCode()==22)
			{
				Double_t energy = candidate->GetEnergy()*1000.;
				Double_t px = candidate->GetPx();
				Double_t py = candidate->GetPy();
				Double_t pz = candidate->GetPz();
				TVector3 gammaVector(px, py, pz);
				
				//Relative Angle
				Double_t theta_relative = gammaVector.Angle(fragmentVector3);

				//Doppler Correction
				Double_t energyCorrected = energy* (1-betaFragment*TMath::Cos(theta_relative))/TMath::Sqrt(1-betaFragment*betaFragment);
				Gamma_Energy_and_Theta->Fill(theta_relative*180./TMath::Pi(), energyCorrected);
				Gamma_Energy_and_Theta->Fill(theta_relative*180./TMath::Pi(), energy);
				Distribution_Corrected_Gamma->Fill(energyCorrected);
				
				Gamma_Energy_Distribution->Fill(energy);
				continue;
			}
			//Protons
			if(IS_PN)
			{
				cout<<"\nThis is (p,np) scattering" <<endl;
			}
			else
			{
				cout<<"\nThis is (p,2p) scattering" <<endl;
			}
		//	if(candidate->GetPdgCode()==2212 ||candidate->GetPdgCode()==2112) 
		//	if(candidate->GetPdgCode()==2212)
		//	{
				cout << "\n Proton Px = " << candidate->GetPx() << "\t PdG Candidate: " << candidate->GetPdgCode() << endl;
				if(!proton[0])
					proton[0] = candidate;
				else
					proton[1] = candidate;
		//		continue;
		//	}
			//Randomly Swap the protons
//			if(randomGen.Integer(2) == 1)
//			{
//				std::swap(proton[0], proton[1]);
//			}	
		}
		if(!proton[0] || !proton[1])
		{
			cout << "\nERROR! No protons found!\n";
			continue;
		}
		int randomIndex = randomGen.Integer(2); //random btw 0 or 1
	//	R3BMCTrack * selectedProton = proton[randomIndex];
		int random2 = 1-randomIndex;
	//	for(auto k=0; k<2; k++)
	//	{

		proton_v3[0].SetXYZ(proton[randomIndex]->GetPx(), 
				proton[randomIndex]->GetPy(), 
				proton[randomIndex]->GetPz());

		proton_v3[1].SetXYZ(proton[random2]->GetPx(), proton[random2]->GetPy(), proton[random2]->GetPz());
	
		PhiProton_v3[0].SetXYZ(proton[randomIndex]->GetPx(), proton[randomIndex]->GetPy(), proton[randomIndex]->GetPz());

		PhiProton_v3[1].SetXYZ(proton[random2]->GetPx(), proton[random2]->GetPy(),proton[random2]->GetPz());
		//	}		
		Double_t T1 =TMath::Sqrt(pow(proton_v3[0].Mag()*1000.,2) + pow(938,2))-938;
		Double_t T2 =TMath::Sqrt(pow(proton_v3[1].Mag()*1000.,2) + pow(938,2))-938;
		if( proton_v3[0].Theta()*180./TMath::Pi() > 20 && proton_v3[1].Theta()*180./TMath::Pi() > 20){
			Theta_Correlation_califa_cut->Fill(proton_v3[0].Theta()*180./TMath::Pi(), proton_v3[1].Theta()*180./TMath::Pi());
			Phi_Correlation_califa_cut->Fill(PhiProton_v3[0].Phi()*180./TMath::Pi(), PhiProton_v3[1].Phi()*180./TMath::Pi());
			Relative_Theta_califa_cut->Fill(proton_v3[0].Angle(proton_v3[1])*180./TMath::Pi());
		}
		Theta_Correlation->Fill(proton_v3[0].Theta()*180./TMath::Pi(), proton_v3[1].Theta()*180./TMath::Pi());
		Phi_Correlation->Fill(PhiProton_v3[0].Phi()*180./TMath::Pi(), PhiProton_v3[1].Phi()*180./TMath::Pi());
		phiproton1->Fill(PhiProton_v3[0].Phi()*180./TMath::Pi());
		phiproton2->Fill(PhiProton_v3[1].Phi()*180./TMath::Pi());
		Theta_1->Fill(proton_v3[0].Theta()*180./TMath::Pi());
		Theta_2->Fill(proton_v3[1].Theta()*180./TMath::Pi());
		Relative_Theta->Fill(proton_v3[0].Angle(proton_v3[1])*180./TMath::Pi());
		OpangvsTheta->Fill(proton_v3[0].Theta()*180./TMath::Pi(),proton_v3[0].Angle(proton_v3[1])*180./TMath::Pi());
		Difference_Phi->Fill(fabs((proton_v3[0].Phi()-proton_v3[1].Phi()))*180./TMath::Pi());
		Kinetic_Energies_Correlation->Fill(T1, T2);
//		Angular_Energy_Distribution->Fill(T1, proton_v3[0].Theta());
//		Angular_Energy_Distribution_2->Fill(T2, proton_v3[1].Theta());
	}	

	TFile* file = new TFile("histogram_manuel.root", "RECREATE");
	TCanvas * c1 = new TCanvas("c1","c1",1000,1000);
	c1->Divide(4,4);
	c1->cd(1);
	Theta_Correlation->Draw("colz");
	Theta_Correlation->Write();
	Theta_Correlation_califa_cut->Write();
	OpangvsTheta->Write();
//	Theta_Correlation->Write();	

	c1->cd(2);
	Phi_Correlation->Draw("colz");
	Phi_Correlation->Write();
	Phi_Correlation_califa_cut->Write();

	c1->cd(3);
	Relative_Theta->Draw("colz");
	Relative_Theta->Write();
	Relative_Theta_califa_cut->Write();

	c1->cd(4);
	Theta_1->Draw("colz");
        Theta_1->Write();
        
	c1->cd(5);
	Theta_2->Draw("colz");	
        Theta_2->Write();
	c1->cd(6);
	phiproton1->Draw("colz");

	c1->cd(7);
	phiproton2->Draw("colz");

	c1->cd(8);
	Difference_Phi->Draw("colz");

	c1->cd(9);
	Gamma_Energy_and_Theta->Draw("colz");

	c1->cd(10);
	Distribution_Corrected_Gamma->Draw("colz");

	c1->cd(11);
	Gamma_Energy_Distribution->Draw("hist");

	c1->cd(12);
	Kinetic_Energies_Correlation->Draw("colz");
	
//	Angular_Energy_Distribution->SetLineColor(kRed);
//	Angular_Energy_Distribution_2->SetLineColor(kBlue);
//	c1->cd(13);
//	Angular_Energy_Distribution->Draw();
//	Angular_Energy_Distribution_2->Draw("same");
//	legend->AddEntry(Angular_Energy_Distribution, "First Particle", "l");
//	legend->AddEntry(Angular_Energy_Distribution_2, "Second Particle", "l");
	
	file->Close();	

}
