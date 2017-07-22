{
//Configuracion inicial
gROOT->Reset();
//gROOT->SetStyle("Plain");
gStyle->SetOptStat(0);
//gStyle->SetPalette(52,0);//Faltaba poner el parametro (52,0)

const Int_t NRGBs = 5;
const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
   Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);


//Archivo GAMOS para trabajar
TFile *f=new TFile("dato04.root");
TString hnam="RTPSPDoseHistos: Dose XY_merged";
TH2F *h0 = (TH2F *)f->Get(hnam);
TH2F *h1 = new TH2F(*h0);//Aqui se almacenara la imagen raw de laboratorio

const int R=sqrt(h1->GetSize());//El valor de N es 256
const int N=R-2;//Menos over flow

//Codigo C++ para obtener el contenido de los datos experimentales y ponerlos en un array
    string line;
    int row,col;
    float my_array[N][N];
    ifstream pFile ("30mAs.txt");//Nombre del archivo experimental
    if (pFile.is_open())
    {
        row=0;
        while(!pFile.eof())
        {
            getline(pFile, line);
            stringstream ss(line);
            col=0;
            while(ss >> my_array[row][col])
            {
                col++;
            }
            row++;
        } 
        pFile.close();
    }
    else 
        cout << "Verificar que el archivo este donde necesitamos"<<endl;
//Se recorre toda la imagen para fijar el contenido del array en cada bin

for(int k=0;k<N;k++){
 for(int l=0;l<N;l++){
  h1->SetBinContent(k,l,my_array[k][l]); 
  }
  }
//Escribir resultado en archivo .root
TFile h11("30mAs.root","new");
TString hnam="RTPSPDoseHistos: Dose XY_merged";
h1->Write();

h1->Draw("COLZ");

}
