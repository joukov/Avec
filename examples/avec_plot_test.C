void avec_plot_test()
{
  Avec v1=grand(Avec(1000));
  Avec v2=grand(Avec(1000));

  TCanvas * cv=new TCanvas("avec_plot_test","avec_plot_test");
  cv->Divide(2,2);
  cv->cd(1);
  avec_plot(v1,30);
  cv->cd(2);
  avec_plot(v1,v2,30,30);
  cv->cd(3);
  avec_plot(v2*v1, 30);
}
