#ifndef __HIST_STYLE__
#define __HIST_STYLE__
#include "TEfficiency.h"

template<typename T>
void setScaleToPad(T* hist, double scale){
    if constexpr (!(std::is_same<T, TEfficiency>::value)){
        hist->GetXaxis()->SetTitleSize(hist->GetXaxis()->GetTitleSize()*scale);
        hist->GetXaxis()->SetLabelSize(hist->GetXaxis()->GetLabelSize()*scale);

        hist->GetYaxis()->SetTitleSize(hist->GetYaxis()->GetTitleSize()*scale);
        hist->GetYaxis()->SetLabelSize(hist->GetYaxis()->GetLabelSize()*scale);
    }
    if constexpr (std::is_same<T, TEfficiency>::value){
        hist->Draw();
        gPad->Update();
        hist->GetPaintedGraph()->GetXaxis()->SetTitleSize(hist->GetPaintedGraph()->GetXaxis()->GetTitleSize()*scale);
        hist->GetPaintedGraph()->GetXaxis()->SetLabelSize(hist->GetPaintedGraph()->GetXaxis()->GetLabelSize()*scale);

        hist->GetPaintedGraph()->GetYaxis()->SetTitleSize(hist->GetPaintedGraph()->GetYaxis()->GetTitleSize()*scale);
        hist->GetPaintedGraph()->GetYaxis()->SetLabelSize(hist->GetPaintedGraph()->GetYaxis()->GetLabelSize()*scale);
    }
};


#endif