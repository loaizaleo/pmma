{
gROOT->Reset();
gStyle->SetOptStat(0);
//gStyle->SetPalette(52,0);//Faltaba poner el parametro (52,0)
//TFile *f0=new TFile("reg11_flat_normal_median.root");

const Int_t NRGBs = 5;
const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
   Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);


TFile *f=new TFile("30mAs_FF.root");
TString hnam="RTPSPDoseHistos: Dose XY_merged";
TH2F *h = (TH2F *)f->Get(hnam);
int d=10;//Ancho del strip
const int N=256,R=250;
float max=0;
TH2F *h2 = new TH2F("h2","Simulacion Fantoma Acreditacion, reg07, M=2X, di= 30cm",R,-7.0,7.0,d,0.,d);
max=h->GetBinContent(h->GetMaximumBin());
float binval[N][N];
for(int i=0;i<N;i++){
 for(int j=0;j<N;j++){
  binval[i][j]=h->GetBinContent(i,j)/max;
}
  }
//Definir desde cual biny observar
int biny=179;
for(int k=0;k<R;k++){
 for(int l=biny;l<(biny+d);l++){
  h2->SetBinContent(k+1,l-biny+1,binval[k+6][l]/d); 
}}

 TH1D *px2 = h2->ProjectionX("px2",0,d);

int cuenta=0;
float noise[R],threshold=0.529975;
float SNR[R],maxbins=0,minbins=0;

maxbins=px2->GetBinContent(px2->GetMaximumBin());
minbins=px2->GetBinContent(px2->GetMinimumBin());
TH1F* hh=new TH1F("hh","signal to noise",R,minbins,maxbins);
for(int m=0;m<R;m++){
SNR[m]=px2->GetBinContent(m);
if (SNR[m]<threshold){
cuenta+=1;
}
hh->Fill(SNR[m]);
}
TH1F* ni=new TH1F("nn","NOISE",R,minbins,threshold);
TH1F* sig=new TH1F("sig","SIGNAL",R-cuenta,threshold,maxbins);

//Primero toca ordenar los datos de SNR
float *ordenado = SNR;
for (int o = R-1; o > 0; --o) {
        for (int p = 0; p < o; ++p) {
            if (ordenado[p] > ordenado[p+1]) {
                float dTemp = ordenado[p];
                ordenado[p] = ordenado[p+1];
                ordenado[p+1] = dTemp;
            }        

}}
for(int v=0;v<=cuenta;v++){
noise[v]=SNR[v];
ni->Fill(noise[v]);
}

for(int vv=cuenta;vv<R;vv++){
noise[vv]=SNR[vv];
sig->Fill(noise[vv]);
}
float snr;
snr=(sig->GetMean()-ni->GetMean())/(3*ni->GetRMS());
cout<<"Signal to Noise Ratio = "<<snr<<endl;
cout<<"Mean signal= "<<sig->GetMean()<<endl;
cout<<"Mean noise= "<<ni->GetMean()<<endl;
//Escribir resultado en archivo .root
//TFile hh("strip00.root","new");
//px2->Write();


TCanvas *cc;
cc =new TCanvas("t4_t3","t4_t3 ",600,400);
//cc->Divide(2,2);

//cc->cd(1);
// h->Draw("COLZ");
//cc->cd(2); 
px2->Draw();
px2->GetXaxis()->SetTitle("Posicion X[mm]");
px2->GetYaxis()->SetTitle("Conteos");
//TText *latex = new TLatex(0.3,0.3,"snr= 1.37");
//latex->Draw(); 
//cc->cd(3);
// ni->Draw();
//cc->cd(4);
// sig->Draw();
}
