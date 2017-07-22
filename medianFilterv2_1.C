//En este programa se aplica Median Filter a una imagen
//raw para disminuir el ruido
{
//Parametros para entorno inicial 
gROOT->Reset();
//gROOT->SetStyle("Plain");
gStyle->SetOptStat(00);
//gStyle->SetPalette(52,0);//Faltaba poner el parametro (52,0)
const Int_t NRGBs = 5;
const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
   Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

//int NRGBs = 7, NCont = 999;
//gStyle->SetNumberContours(NCont);
//Double_t stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
//Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
//Double_t green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
//Double_t blue[NRGBs]  = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
//TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

//Aqui se define el archivo con el que se va a trabajar
TFile *f0=new TFile("50.6mAs_FF.root");
TString hnam="RTPSPDoseHistos: Dose XY_merged";
TH2F *h0 = (TH2F *)f0->Get(hnam);
TH2F *h1 = new TH2F(*h0);//Aqui se almacenara la imagen filtrada y se le heredan atributos de *h0

const int R=sqrt(h1->GetSize());
const int N=R-0;//Menos over flow
float binval[N][N];

//Se recorre toda la imagen para obtener el contenido de cada bin

for(int i=0;i<N;i++){
 for(int j=0;j<N;j++){
  binval[i][j]=h0->GetBinContent(i,j); 
  }
  }

//Parametros para calcular la mediana


//Se necesita rellenar lo demas con la mediana de los vecinos
for(int k=2;k<N-2;k++){
 for(int l=2;l<N-2;l++){

//Codigo para calcular la mediana 
int iSize=25;
float vecinos[iSize]={binval[k][l],binval[k][l+1],binval[k][l-1],binval[k+1][l],binval[k-1][l],binval[k+1][l+1],binval[k-1][l+1],binval[k-1][l-1],binval[k+1][l-1],binval[k+1][l+2],binval[k+1][l-2],binval[k+0][l+2],binval[k+0][l-2],binval[k-1][l+2],binval[k-1][l-2],binval[k-2][l+2],binval[k-2][l-1],binval[k-2][l],binval[k-2][l-1],binval[k-2][l-2],binval[k+2][l+2],binval[k+2][l+1],binval[k+2][l],binval[k+2][l-1],binval[k+2][l-2]};
float *ordenado = vecinos;

for (int o = iSize - 1; o > 0; --o) {
        for (int p = 0; p < o; ++p) {
            if (ordenado[p] > ordenado[p+1]) {
                float dTemp = ordenado[p];
                ordenado[p] = ordenado[p+1];
                ordenado[p+1] = dTemp;
            }
        }}

 h1->SetBinContent(k,l,ordenado[12]);
 }
  }
//Codigo para arreglar los bordes

float binsval[N][N];
float ccv=0;
for(int oi=0;oi<N;oi++){
for (int io=0;io<N;io++){
  binsval[oi][io]=h1->GetBinContent(oi,io);
ccv+=binsval[oi][io];
}}//End

for(int m=0;m<N;m++){
for (int n=0;n<N;n++){
if (m<2 || n<2 || m>N-3 || n>N-3){
binsval[m][n]=ccv/pow(N,2);
h1->SetBinContent(m,n,binsval[m][n]);
}
}}
//ENd Codigo para arreglar los bordes

//Escribir resultado en archivo .root
TFile h11("50.6mAs_M.root","new");
h1->Write();

//cout<<N<<endl;
TCanvas *cc;
cc =new TCanvas("t4_t3","hola",900,500);
//cc->SetGrayscale();
cc->Divide(2,1);

cc->cd(1);
 h0->Draw("COLZ");
cc->cd(2);
 h1->Draw("colz");


}



