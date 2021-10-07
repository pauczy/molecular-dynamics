void hist()
{
  TH1D *hx = new TH1D("hx", "hx", 200, -10., 10.);
  TH1D *hy = new TH1D("hy", "hy", 200, -10., 10.);
  TH1D *hz = new TH1D("hz", "hz", 200, -10., 10.);

  ifstream fin;

  int n = 0;
  fin.open("params.txt");
  fin >> n;
  fin.close();
  int N = n*n*n;
  

  fin.open("out_mom.txt");
  double px, py, pz;
  for(int i = 0; i < N; i++)
  {
    fin >> px >> py >> pz;
    hx->Fill(px);
    hy->Fill(py);
    hz->Fill(pz);
  }

  fin.close();
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  hx->Draw();
  c1->cd(2);
  hy->Draw();
  c1->cd(3);
  hz->Draw();

}
