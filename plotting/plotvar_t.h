#ifndef _PLOTVAR_T_H
#define _PLOTVAR_T_H

#if 0
PlotVar_t(char *inpv,double inmaxr,double inminr,int innbin,int inslog,char *inxl,char *inoutf,char *inoutf2,double inamax,double inamin,int inanb,int indl) :
  plotvar(inpv),
  MAXRange(inmaxr),MINRange(inminr),NBINS(innbin),
  slog(inslog),xlabel(inxl),outfile(inoutf),outfile2(inoutf2),
  AMAXRange(inamax),AMINRange(inamin),ANBINS(inanb),
  drawleg(indl) {}
#endif
struct plotVar_t {
  TString plotvar;
  double MINRange;
  double MAXRange;
  int    NBINS;
  int    slog;
  TString xlabel;
  TString outfile;
  double AMINRange;
  double AMAXRange;
  int    ANBINS;
  int  drawleg;
};

#endif // _PLOTVAR_T_H
