//Para generar archivo .root flat Field
{
  gROOT->Reset();
gStyle->SetOptStat(0);
//Para colorear en escala de gris
const Int_t NRGBs = 5;
const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
   Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
//Para cargar archivos .root y asignarlos a funciones TH2F
TFile *f=new TFile("30mAs_FF.root");//histograma OB
TString hnam="RTPSPDoseHistos: Dose XY_merged";
TH2F *h = (TH2F *)f->Get(hnam);
int N=h->GetSize(),overflow=2;
int xlow=-7,xmax=7,ylow=-7,ymax=7;
float maxbins=0;
maxbins=h->GetBinContent(h->GetMaximumBin());;
TH2F *h1 = new TH2F("","Simulacion Fantoma de Acreditacion, reg7, M=2X, di=30cm",sqrt(N)-overflow,xlow,xmax,sqrt(N)-overflow,ylow,ymax);
//Normalizar dividiendo entre el maximo
float temp=0,cont=0;
for(int i=0;i<N;i++){
temp=h->GetBinContent(i+1)/maxbins;
cont+=temp;
h1->SetBinContent(i+1,temp);
}
cout<<cont<<"\t"<<maxbins<<endl;
//End for
TCanvas *hola;
hola =new TCanvas("t4_t3","t4_t3 ",600,600);
h1->Draw("colz");
}//End root
