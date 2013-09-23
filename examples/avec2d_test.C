void avec2d_test()
{
  Avec y1=exp(-0.3*pow(Avec(20,-10,10),2));
  Avec2D z2=vtrans(Avec2D(20,y1))*y1;
  TCanvas* cv=new TCanvas("avec2d_test","avec2d_test");
  avec_draw(z2);
}

