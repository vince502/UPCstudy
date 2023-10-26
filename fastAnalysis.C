#include "TFile.h"
#include "TTree.h"
#include "TTreeReaderValue.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "ROOT/TSeq.hxx"
#include <vector>
#include "../info.h"
#include "treeIO.h"

void fastAnalysis(){

	ParticleTreeReader r1("diMu_ana_ppRECO_MC.root", "diMuAna/ParticleTree");
	ParticleTreeReader r2("diMu_ana_HBS_LbyL_MC.root", "diMuAna/ParticleTree");
	TCanvas* c1 = new TCanvas("c1", "", 1400, 1200);
	c1->Divide(2,2);
	TH1D* hdau1 = new TH1D("hdau1","",400,0,20);
	TH1D* hdau2 = new TH1D("hdau2","",400,0,20);
	TH1D* h = new TH1D("h","",400,0,20);
	TH1D* h2 = new TH1D("h2","",160, -4.0, 4.0);
	TH1D* h3 = new TH1D("h3","",628, -3.14, 3.14);
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
		TH1D* _h = (TH1D*) h->Clone(Form("_h_%i", id));
		TH1D* _h2 = (TH1D*) h2->Clone(Form("_h2_%i", id));
		TH1D* _h3 = (TH1D*) h3->Clone(Form("_h3_%i", id));
		TH1D* _hID = (TH1D*) hID->Clone(Form("_hID_%i", id));
		r->theReader->Restart();
		while(r->theReader->Next()){
			int mother= -1;
			int med1 = -1;
			int med2 = -1;

			std::vector<std::vector<unsigned int> > dau= {};

			for( auto idx : ROOT::TSeqI((**(r->cand_mass)).size())){
				if( (**(r->cand_pdgId))[idx] == 443){
					mother = idx;
				}
				_hID->Fill( (**(r->cand_pdgId))[idx]);
			}
			if(mother == -1) continue;
			_h->Fill((**(r->cand_pTDau))[0][0]);
			_h->Fill((**(r->cand_pTDau))[0][1]);
			_h2->Fill((**(r->cand_etaDau))[0][0]);
			_h2->Fill((**(r->cand_etaDau))[0][1]);
			_h3->Fill((**(r->cand_phiDau))[0][0]);
			_h3->Fill((**(r->cand_phiDau))[0][1]);
		};
		return std::vector<TH1D*>{_h, _h2, _h3};
	};
//	hdau->SetLineColor(kRed);
	auto res1vec = theJob(&r1, 1);
	for( auto res1 : res1vec){
		res1->Scale(1./res1->GetSum());
		res1->SetLineColor(kBlue);
		res1->SetLineWidth(5);
	}

	auto res2vec = theJob(&r2, 2);
	for( auto res2 : res2vec){
		res2->Scale(1./res2->GetSum());
		res2->SetLineColor(kRed);
	}


//	hdau->Draw("pe");
//	h->Draw("same");
	c1->cd(1);
	res1vec[0]->Draw();
	res2vec[0]->Draw("same");
	c1->cd(2);
	res1vec[1]->Draw();
	res2vec[1]->Draw("same");
	c1->cd(3);
	res1vec[2]->Draw();
	res2vec[2]->Draw("same");
	
	// h->Draw("");
	// c1->cd(2);
	// hdau1->Draw("");
//	hID->Draw();


};
