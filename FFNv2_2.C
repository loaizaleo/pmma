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
TFile *f=new TFile("OB_30mAs_60cm.root");//histograma OB
TFile *f1=new TFile("30mAs.root");//histograma raw
TFile *f2=new TFile("imagen_DI.root");//histograma DI
TString hnam="RTPSPDoseHistos: Dose XY_merged";
TH2F *h = (TH2F *)f->Get(hnam);
TH2F *h1 = (TH2F *)f1->Get(hnam);
TH2F *h2 = (TH2F *)f2->Get(hnam);
//Declaracion de variables
float maxraw=0,sumraw=0,meanraw=0;//Para datos RAW
float maxOB=0,sumOB=0,meanOB=0;//Para datos OB
const float maxumbral=0.95,minumbral=0.5;
//Declaracion de matrices
int R=sqrt(h1->GetSize());
int N=R-2, contadoraw=0,contadorOB=0;
float raw[N][N],normalraw[N][N],OB[N][N],normalOB[N][N],DI[N][N];
float ctrlRAW, ctrlOB;
maxOB=h->GetBinContent(h->GetMaximumBin());
maxraw=h1->GetBinContent(h1->GetMaximumBin());
//Primero calcular el promedio de las matrices RAW y OB
for(int ii=0;ii<N;ii++){
 for(int jj=0;jj<N;jj++){
OB[ii][jj]=h->GetBinContent(ii,jj);//Se asignan los valores del .root OB a una matriz
h->SetBinContent(ii,jj,OB[ii][jj]/maxOB);//Se normaliza h 
raw[ii][jj]=h1->GetBinContent(ii,jj);//Se asignan los valores del .root RAW a una matriz
h1->SetBinContent(ii,jj,raw[ii][jj]/maxraw);//Se normaliza h1 
DI[ii][jj]=h2->GetBinContent(ii,jj);//Se asignan los valores del .root DI a una matriz
normalOB[ii][jj]=h->GetBinContent(ii,jj);
normalraw[ii][jj]=h1->GetBinContent(ii,jj);
sumOB+=h->GetBinContent(ii,jj,normalOB[ii][jj]);
sumraw+=h1->GetBinContent(ii,jj,normalraw[ii][jj]);
}}
meanraw=sumraw/pow(N,2);
meanOB=sumOB/pow(N,2);
//Identificar bins muertos y muriendo de matriz RAW
for(int i=0;i<N;i++){
 for(int j=0;j<N;j++){
if (normalraw[i][j]>=0 && normalraw[i][j]<=minumbral || normalraw[i][j]>=maxumbral && normalraw[i][j]<=1){//identifica pixel que estan muriendo.
contadoraw+=1;
h1->SetBinContent(i,j,meanraw);
}//end del if
else{
h1->SetBinContent(i,j,normalraw[i][j]);
}//end del else
normalraw[i][j]=h1->GetBinContent(i,j);
}}//end del for
//cout<<"promedio = "<<meanraw<<"\tbin pichos = "<<contador<<endl;

//Asignar promedio de ocho vecinos
float mini[3][3];
for(int k=0;k<N;k++){
 for(int l=0;l<N;l++){
ctrlRAW=normalraw[k][l];
  if (k-1<0 || l-1<0 || k+1>255 || l+1>255){
h1->SetBinContent(k,l,normalraw[k][l]);
}//end del if
 else{

//Pero sin contar los pichos
if (ctrlRAW==meanraw){
float promini=0;
int cc=0;
for(int m=0;m<3;m++){
 for(int n=0;n<3;n++){

mini[m][n]=normalraw[k+(m-1)][l+(n-1)];
if(mini[m][n]==meanraw){
cc+=1;
mini[m][n]=0;
}
promini+=mini[m][n];

}}
if(cc==9){
normalraw[k][l]=meanraw;
}
else{normalraw[k][l]=promini/(9-cc);}
//cout<<"mini ="<<mini[0][1]<<" \t"<<"promini ="<<promini/(9-cc)<<" \t"<<cc<<" \t"<<k<<" \t"<<l<<endl;
}
h1->SetBinContent(k,l,normalraw[k][l]);
else{
h1->SetBinContent(k,l,normalraw[k][l]);
normalraw[k][l]=h1->GetBinContent(k,l);//Por ultimo se vuelve a poner todo en forma de matriz
}
}
    }}

//Identificar bins muertos y muriendo de matriz OB
for(int o=0;o<N;o++){
 for(int p=0;p<N;p++){
if (normalOB[o][p]>=0 && normalOB[o][p]<=minumbral || normalOB[o][p]>=maxumbral && normalOB[o][p]<=1){//identifica pixel que estan muriendo.
contadorOB+=1;
h->SetBinContent(o,p,meanOB);
}//end del if
else{
h->SetBinContent(o,p,normalOB[o][p]);
}//end del else
normalOB[o][p]=h->GetBinContent(o,p);
}}//end del for
//Asignar promedio de ocho vecinos
float miniOB[3][3];
for(int kk=0;kk<N;kk++){
 for(int ll=0;ll<N;ll++){
ctrlOB=normalOB[kk][ll];
  if (kk-1<0 || ll-1<0 || kk+1>255 || ll+1>255){
h->SetBinContent(kk,ll,normalOB[kk][ll]);
}//end del if
 else{

//Pero sin contar los pichos
if (ctrlOB==meanOB){
float prominiOB=0;
int ccOB=0;
for(int mm=0;mm<3;mm++){
 for(int nn=0;nn<3;nn++){

miniOB[mm][nn]=normalOB[kk+(mm-1)][ll+(nn-1)];
if(miniOB[mm][nn]==meanOB){
ccOB+=1;
miniOB[mm][nn]=0;
}
prominiOB+=miniOB[mm][nn];
}}
if(ccOB==9){
normalOB[kk][ll]=meanOB;
}
else{normalOB[kk][ll]=prominiOB/(9-ccOB);}
//cout<<"mini ="<<miniOB[0][1]<<" \t"<<"prominiOB ="<<prominiOB/(9-ccOB)<<" \t"<<ccOB<<" \t"<<kk<<" \t"<<ll<<endl;
}
h->SetBinContent(kk,ll,normalOB[kk][ll]);
else{
h->SetBinContent(kk,ll,normalOB[kk][ll]);
normalOB[kk][ll]=h->GetBinContent(kk,ll);//Por ultimo se vuelve a poner todo en forma de matriz
}
}
 }}
//Aplicar correccion Flat Field FF
float FF[N][N],max=0,neg=0;
for(int oo=0;oo<N;oo++){
 for(int pp=0;pp<N;pp++){
FF[oo][pp]=(normalraw[oo][pp]-DI[oo][pp])/(normalOB[oo][pp]-DI[oo][pp]);
h2->SetBinContent(oo,pp,FF[oo][pp]);
}}
//Para obtener el negativo de la imagen
max=h2->GetBinContent(h2->GetMaximumBin());
int nbins=h2->GetSize()-2; // Minus Underflow and Overflow
for (int w=0;w<nbins;w++){
neg=max-(h2->GetBinContent(w));
h2->SetBinContent(w,neg);
}
float maxx=0;
maxx=h2->GetBinContent(h2->GetMaximumBin());
float binval[N][N];
for(int ff=0;ff<N;ff++){
 for(int g=0;g<N;g++){
  binval[ff][g]=h2->GetBinContent(ff,g)/maxx;
h2->SetBinContent(ff,g,binval[ff][g]);
}
  }

//Codigo para arreglar los bordes, primero necesito contabilizar el promedio de los buenos
float noise[N],threshold=0.3901;
const int RR=sqrt(h2->GetSize());
float binsval[RR][RR], maxbins, minbins;
float ccv=0;
for(int oi=0;oi<RR;oi++){
for (int io=0;io<RR;io++){
  binsval[oi][io]=h1->GetBinContent(oi+1,io+1);
ccv+=binsval[oi][io];
}}//End
for(int mh=0;mh<RR;mh++){
for (int nh=0;nh<RR;nh++){
if (mh<1 || nh<1 || mh>RR-4 || nh>RR-4){
binsval[mh][nh]=ccv/(2*pow(RR,2));
h2->SetBinContent(mh+1,nh+1,binsval[mh][nh]);
}
}}

//End
//for (int oooo=0;oooo<RR;oooo++){
//if (binsval[oooo]>0.91){
//binsval[oooo]=ccv/RR;
//h2->SetBinContent(oooo,binsval[oooo]);
//}
//}//ENd Codigo para arreglar los bordes
int cuenta=0;
//TH1F* ni=new TH1F("ni","RUIDO",N,minbins,threshold);
//TH1F* sig=new TH1F("sig","SEÑAL",N-cuenta,threshold,maxbins);


//Para obtener el strip o perfil

int d=5;//Ancho del strip
TH2F *h3 = new TH2F("h2","strip de 256 por 10 pixels",256,-7.0,7.0,d,0.,d);
int biny=137;//Definir desde cual biny observar
for(int a=0;a<N;a++){
 for(int s=biny;s<(biny+d);s++){
  h3->SetBinContent(a+1,s-biny+1,binval[a][s]/d); 
}}
 TH1D *px2 = h3->ProjectionX("px2",0,d);

int nbinss=px2->GetSize()-2;
float SNR[nbinss],maxbins=0,minbins=0;
maxbins=px2->GetBinContent(px2->GetMaximumBin());
minbins=px2->GetBinContent(px2->GetMinimumBin());
TH1F* hh=new TH1F("hh","signal to noise",nbinss,minbins,maxbins);
for(int m=0;m<N;m++){
SNR[m]=px2->GetBinContent(m);
if (SNR[m]<threshold){
cuenta+=1;
}
hh->Fill(SNR[m]);
}

TH1F* ni=new TH1F("nn","NOISE",N,minbins,threshold);
TH1F* sig=new TH1F("sig","SIGNAL",N-cuenta,threshold,maxbins);
//Primero toca ordenar los datos de SNR
float *ordenado = SNR;
for (int o = N-1; o > 0; --o) {
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

for(int vv=cuenta;vv<N;vv++){
noise[vv]=SNR[vv];
sig->Fill(noise[vv]);
}
float snr;
snr=(sig->GetMean()-ni->GetMean())/(3*ni->GetRMS());
cout<<"Signal to Noise Ratio = "<<snr<<endl;



TCanvas *hola;
hola =new TCanvas("t4_t3","t4_t3 ",650,650);
hola->Divide(2,2);

TFile ppt("30mAs_FF.root","new");
h2->Write();


hola->cd(1);
 h2->Draw("COLZ");
hola->cd(2);
 px2->Draw();
hola->cd(3);
 ni->Draw();
hola->cd(4);
 sig->Draw();

}//Aqui se termina root

