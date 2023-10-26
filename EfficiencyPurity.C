#include "TFile.h"
#include "TTree.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "ROOT/TSeq.hxx"
#include <vector>
#include "info.h"
#include "treeIO.h"
#include "histStyler.h"
#include "TEfficiency.h"
#include "TPad.h"

void EfficiencyPurity(const int pid){
	// ROOT::EnableImplicitMT(4);
	const int HLTRIGGER = 12;
	gStyle->SetOptStat(0);
	gStyle->SetCanvasPreferGL(true);
	string tree, inputFile_pp, inputFile_fw;
	string classStudy = "";
	if (pid==0) {
		tree = "trackAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_DF_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_DF_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "DF_trak";
	}
	else if (pid==1) {
		tree = "muonAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_mu_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_mu_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "QED_mu";
	}
	else if (pid==2) {
		tree = "elecAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_el_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_el_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "QED_el";
	}
	else if (pid==3) {
		tree = "lowPtElecAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_el_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/STARLIGHT_5p36TeV_2023Run3/ParticleAnalyzer_QED_el_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "QED_LowpTel";
	}
	else if (pid==4) {
		tree = "phoAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/SUPERCHIC_5p36TeV_2023Run3/ParticleAnalyzer_LbL_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/SUPERCHIC_5p36TeV_2023Run3/ParticleAnalyzer_LbL_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "photon";
	}
	else if (pid==5) {
		tree = "convAna/ParticleTree";
		inputFile_pp = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/SUPERCHIC_5p36TeV_2023Run3/ParticleAnalyzer_ChiC_ppRec_MC_HIRun2023_2023_10_24.root";
		inputFile_fw = "/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/ParticleAnalyzer/MC/2023_10_24/SUPERCHIC_5p36TeV_2023Run3/ParticleAnalyzer_ChiC_fwRec_MC_HIRun2023_2023_10_24.root";
		classStudy = "conversion";
	}

	// ParticleTreeReader r1("singleMuInclTrk_ana_ppRECO_MC.root", "trackAna/ParticleTree");
	// ParticleTreeReader r2("singleMuInclTrk_ana_HBS_LbyL_MC.root", "trackAna/ParticleTree");
	ParticleTreeReader r1(inputFile_pp.c_str(), tree.c_str());
	ParticleTreeReader r2(inputFile_fw.c_str(), tree.c_str());
	//ParticleTreeReader r1("diMu_ana_ppRECO_MC.root", "diMuAna/ParticleTree");
	//ParticleTreeReader r2("diMu_ana_HBS_LbyL_MC.root", "diMuAna/ParticleTree");
	TCanvas* c1 = new TCanvas("c1", "", 1400, 1200);
	c1->Divide(2,2);
	TH1D* hdau1 = new TH1D("hdau1","",400,0,20);
	TH1D* hdau2 = new TH1D("hdau2","",400,0,20);
	TH1D* h1Pass = new TH1D("h1pass","",100,0,10);
	TH1D* h2Pass = new TH1D("h2pass","",160, -4.0, 4.0);
	TH1D* h3Pass = new TH1D("h3pass","",628, -3.14, 3.14);

	TH2D* h12Pass = new TH2D("h12pass","",80, -4.0, 4.0,100,0,10);

	TH1D* h1Fail = (TH1D*) h1Pass->Clone("h1fail");
	TH1D* h2Fail = (TH1D*) h2Pass->Clone("h2fail");
	TH1D* h3Fail = (TH1D*) h3Pass->Clone("h3fail");
	TH2D* h12Fail = (TH2D*) h12Pass->Clone("h12fail");

	TH1D* h1Total = (TH1D*) h1Pass->Clone("h1total");
	TH1D* h2Total = (TH1D*) h2Pass->Clone("h2total");
	TH1D* h3Total = (TH1D*) h3Pass->Clone("h3total");
	TH2D* h12Total = (TH2D*) h12Pass->Clone("h12total");

	TH1D* hID = new TH1D("hID","",1000,0,1000);

	unsigned long count = 0;
	const float dedxThres1 = 2.6;
	const float dedxThres2 = 2.6;

	auto  getPID = [&](float p, float dedx){
		return std::vector<int>{e(p, dedx), pi(p, dedx) ,ka(p, dedx), pr(p, dedx)};
	};
	auto theJob = [&](ParticleTreeReader* r, int id){
		TH1D* _hdau1 = (TH1D*) hdau1->Clone(Form("_hdau1_%i", id));
		TH1D* _hdau2 = (TH1D*) hdau2->Clone(Form("_hdau2_%i", id));
		TH1D* _h1Pass = (TH1D*) h1Pass->Clone(Form("hpass_%i", id));
		TH1D* _h2Pass = (TH1D*) h2Pass->Clone(Form("hpass_%i", id));
		TH1D* _h3Pass = (TH1D*) h3Pass->Clone(Form("hpass_%i", id));
		TH2D* _h12Pass = (TH2D*) h12Pass->Clone(Form("h12pass_%i", id));
		TH1D* _h1Fail = (TH1D*) h1Pass->Clone(Form("h1fail_%i", id));
		TH1D* _h2Fail = (TH1D*) h2Pass->Clone(Form("h2fail_%i", id));
		TH1D* _h3Fail = (TH1D*) h3Pass->Clone(Form("h3fail_%i", id));
		TH2D* _h12Fail = (TH2D*) h12Pass->Clone(Form("h12fail_%i", id));
		TH1D* _h1Total = (TH1D*) h1Pass->Clone(Form("h1total_%i", id));
		TH1D* _h2Total = (TH1D*) h2Pass->Clone(Form("h2total_%i", id));
		TH1D* _h3Total = (TH1D*) h3Pass->Clone(Form("h3total_%i", id));
		TH2D* _h12Total = (TH2D*) h12Pass->Clone(Form("h12total_%i", id));
		TH1D* _hID = (TH1D*) hID->Clone(Form("hID_%i", id));
		r->theReader->Restart();
		while(r->theReader->Next()){
			int mother= -1;
			int med1 = -1;
			int med2 = -1;
			// if(!(**(r->evtSel))[0]) continue;
			// if(!(**(r->evtSel))[1]) continue;
			// if(!(**(r->passHLT))[HLTRIGGER]) continue;

			//std::vector<std::vector<unsigned int> > dau= {};

			//for( auto idx : ROOT::TSeqI((**(r->cand_mass)).size())){
			//	if( (**(r->cand_pdgId))[idx] == 443){
			//		mother = idx;
			//	}
			//	_hID->Fill( (**(r->cand_pdgId))[idx]);
			//}
			//if(mother == -1) continue;
			if( (**(r->gen_status)).size() != 2) std::cout << "Warning: Gen level muon size != 1 " << std::endl;
			for( auto genIdx : ROOT::TSeqI((**(r->gen_status)).size())){
				if( (**(r->gen_status))[genIdx] != (unsigned int) 1) continue;
				_h1Total->Fill((**(r->gen_pT))[genIdx]);
				_h2Total->Fill((**(r->gen_eta))[genIdx]);
				_h3Total->Fill((**(r->gen_phi))[genIdx]);
				_h12Total->Fill((**(r->gen_eta))[genIdx], (**(r->gen_pT))[genIdx]);
				if( (**(r->gen_candIdx))[genIdx].size() < 1)continue;
				if( (**(r->gen_candIdx))[genIdx][0] > (unsigned int) 10)continue;
				int recIdx = 1000;
				for( auto idx2 : ROOT::TSeqI((**(r->cand_pdgId)).size())){
					if( !( (**(r->cand_status))[idx2] == 1 &&(**(r->cand_pdgId))[idx2] == abs((**(r->gen_pdgId))[genIdx]) && (**(r->cand_charge))[idx2] == (((**(r->gen_pdgId))[genIdx]) < 0 ? 1 : -1)) ) continue;
					recIdx = idx2;
					break;
				}
				if(recIdx < 100){
					// std::cout << Form("\rcand pT Eta Phi : %.3f, %.3f, %.3f, pdgId : %d", (**(r->cand_pT))[recIdx],  (**(r->cand_eta))[recIdx], (**(r->cand_phi))[recIdx], (**(r->cand_pdgId))[recIdx] ) << std::endl;
					// std::cout << Form("gen pT Eta Phi : %.3f, %.3f, %.3f", (**(r->gen_pT))[genIdx],  (**(r->gen_eta))[genIdx], (**(r->gen_phi))[genIdx] ) << std::endl;
				}
				bool pass = (
							   (recIdx < 100 && (**(r->cand_genIdx))[recIdx]==genIdx)
							// && (**(r->trk_isHP))[(**(r->cand_trkIdx))[recIdx]] 
							// && (**(r->muon_isTracker))[(**(r->cand_srcIdx))[recIdx]]
							);
				if( !pass){ 
					_h1Fail->Fill((**(r->gen_pT))[genIdx]);
					_h2Fail->Fill((**(r->gen_eta))[genIdx]);
					_h3Fail->Fill((**(r->gen_phi))[genIdx]);
					_h12Fail->Fill((**(r->gen_eta))[genIdx], (**(r->gen_pT))[genIdx]);
					continue;
				}
				_h1Pass->Fill((**(r->gen_pT))[genIdx]);
				_h2Pass->Fill((**(r->gen_eta))[genIdx]);
				_h3Pass->Fill((**(r->gen_phi))[genIdx]);
				_h12Pass->Fill((**(r->gen_eta))[genIdx], (**(r->gen_pT))[genIdx]);
			}
		};
		for( auto idx :  ROOT::TSeqI(_h1Total->GetNbinsX()) ){
			if( _h1Total->GetBinContent(idx+1)==0 )  _h1Total->SetBinContent(idx+1, 0.000001);
		}
		for( auto idx :  ROOT::TSeqI(_h2Total->GetNbinsX()) ){
			if( _h2Total->GetBinContent(idx+1)==0 )  _h2Total->SetBinContent(idx+1, 0.000001);
		}
		for( auto idx :  ROOT::TSeqI(_h3Total->GetNbinsX()) ){
			if( _h3Total->GetBinContent(idx+1)==0 )  _h3Total->SetBinContent(idx+1, 0.000001);
		}
		for( auto idx1 :  ROOT::TSeqI(_h12Total->GetNbinsX()) ){
			for( auto idx2 :  ROOT::TSeqI(_h12Total->GetNbinsY()) ){
				if( _h12Total->GetBinContent(idx1+1,idx2+1)==0 )  _h12Total->SetBinContent(idx1+1, idx2+1, 0.000001);
			}
		}
		return std::map<std::string, std::tuple< std::vector<TH1D*>, std::vector<TH2D*> > >{
			{"pass"  ,{{_h1Pass, _h2Pass, _h3Pass}, {_h12Pass}}},
			{"fail"  ,{{_h1Fail, _h2Fail, _h3Fail}, {_h12Fail}}},
			{"total" ,{{_h1Total, _h2Total, _h3Total}, {_h12Total}}},
		};
	};
	auto res1vec = theJob(&r1, 1);

	TEfficiency* effpT1 = new TEfficiency(*(std::get<0>(res1vec["pass"])[0]), *(std::get<0>(res1vec["total"])[0]));
	TEfficiency* effeta1 = new TEfficiency(*(std::get<0>(res1vec["pass"])[1]), *(std::get<0>(res1vec["total"])[1]));
	TEfficiency* effphi1 = new TEfficiency(*(std::get<0>(res1vec["pass"])[2]), *(std::get<0>(res1vec["total"])[2]));

	TEfficiency* effpTeta1 = new TEfficiency(*(std::get<1>(res1vec["pass"])[0]), *(std::get<1>(res1vec["total"])[0]));
	effpT1->SetLineColor(kBlue);
	effeta1->SetLineColor(kBlue);
	effphi1->SetLineColor(kBlue);

	auto res2vec = theJob(&r2, 2);

	TEfficiency* effpT2 = new TEfficiency(*(std::get<0>(res2vec["pass"])[0]), *(std::get<0>(res2vec["total"])[0]));
	TEfficiency* effeta2 = new TEfficiency(*(std::get<0>(res2vec["pass"])[1]), *(std::get<0>(res2vec["total"])[1]));
	TEfficiency* effphi2 = new TEfficiency(*(std::get<0>(res2vec["pass"])[2]), *(std::get<0>(res2vec["total"])[2]));

	TEfficiency* effpTeta2 = new TEfficiency(*(std::get<1>(res2vec["pass"])[0]), *(std::get<1>(res2vec["total"])[0]));
	effpT2->SetLineColor(kRed);
	effeta2->SetLineColor(kRed);
	effphi2->SetLineColor(kRed);

	auto getRatioOfEff = []( TH1D* uu, TH1D* ud, TH1D* du, TH1D* dd  ){
		TH1D* x = ((TH1D*) uu->Clone(Form("_tmp_raio1_%s", uu->GetName())));x->Divide(ud);
		TH1D* y = ((TH1D*) du->Clone(Form("_tmp_raio1_%s", du->GetName())));y->Divide(dd);
		x->Divide(y);
		return x;
	};
	TH1D* ratiopT = getRatioOfEff((std::get<0>(res1vec["pass"])[0]), (std::get<0>(res1vec["total"])[0]), (std::get<0>(res2vec["pass"])[0]), (std::get<0>(res2vec["total"])[0]));
	TH1D* ratioeta = getRatioOfEff((std::get<0>(res1vec["pass"])[1]), (std::get<0>(res1vec["total"])[1]), (std::get<0>(res2vec["pass"])[1]), (std::get<0>(res2vec["total"])[1]));
	TH1D* ratiophi = getRatioOfEff((std::get<0>(res1vec["pass"])[2]), (std::get<0>(res1vec["total"])[2]), (std::get<0>(res2vec["pass"])[2]), (std::get<0>(res2vec["total"])[2]));


	c1->cd(1);
	TPad* p1 = new TPad("p1","", 0., 0.25, 1.0, 1.0);
	TPad* p1sub = new TPad("p1sub","", 0., 0.0, 1.0, 0.25);
	p1->SetBottomMargin(0);
	p1sub->SetTopMargin(0);
	p1->Draw();p1sub->Draw();
	p1->cd();
	setScaleToPad<TEfficiency>(effpT1, 1./0.75);
	effpT1->Draw("ape");
	effpT2->Draw("pe same");
	p1sub->cd();
	setScaleToPad<TH1D>(ratiopT, 1./0.25);
	ratiopT->Draw();

	c1->cd(2);
	TPad* p2 = new TPad("p2","", 0., 0.25, 1.0, 1.0);
	TPad* p2sub = new TPad("p2sub","", 0., 0.0, 1.0, 0.25);
	p2->SetBottomMargin(0);
	p2sub->SetTopMargin(0);
	p2->Draw();p2sub->Draw();
	p2->cd();
	setScaleToPad<TEfficiency>(effeta1, 1./0.75);
	effeta1->Draw("ape");
	effeta2->Draw("pe same");
	p2sub->cd();
	setScaleToPad<TH1D>(ratioeta, 1./0.25);
	ratioeta->GetYaxis()->SetRangeUser(0.8,1.2);
	ratioeta->Draw();

	c1->cd(3);
	TPad* p3 = new TPad("p3","", 0., 0.25, 1.0, 1.0);
	TPad* p3sub = new TPad("p3sub","", 0., 0.0, 1.0, 0.25);
	p3->SetBottomMargin(0);
	p3sub->SetTopMargin(0);
	p3->Draw();p3sub->Draw();
	p3->cd();
	setScaleToPad<TEfficiency>(effphi1, 1./0.75);
	effphi1->Draw("ape");
	effphi2->Draw("pe same");
	p3sub->cd();
	setScaleToPad<TH1D>(ratiophi, 1./0.25);
	ratiophi->GetYaxis()->SetRangeUser(0.9970, 1.003);
	ratiophi->Draw();

	c1->cd(4);
	effpTeta2->Draw("colz");
	
	TFile* output = new TFile(Form("outputEff_%s.root",classStudy.c_str()) ,"recreate");
	output->cd();
	c1->Write();
	effpT1->Write();
	effeta1->Write();
	effphi1->Write();
	effpT2->Write();
	effeta2->Write();
	effphi2->Write();
	effpTeta1->Write();
	effpTeta2->Write();
	ratiopT->Write();
	ratioeta->Write();
	ratiophi->Write();

	c1->SaveAs(Form("plot_%s.pdf", classStudy.c_str()));
	output->Close();

};
