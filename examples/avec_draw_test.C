void avec_draw_test_tmultigraph()
{
  //Create some Avecs for drawing
  Avec x1(100,0,7); Avec y1=sin(x1); Avec y2=cos(x1);
  x1.SetName("x1");
  y1.SetName("sin(x)");
  y2.SetName("cos(x)");
  Avec x2(100,1.5, 8.5); Avec y3=cos(x2);
  Avec dy4=grand(Avec(10),0.2);
  Avec dy5=grand(Avec(10),0.1);
  Avec x3(10,0,7); Avec y4=sin(x3)+dy4; Avec y5=cos(x3)+dy5;
  Avec err4(10,0.2);
  Avec err5(10,0.1);

  std::cout << "Testing drawing of single Avec" << endl;
  // Test for single Avec drawing
  TCanvas* c1=new TCanvas("cv_sing","avec_draw test single");
  c1->Divide(2,2);

  // Single Avec
  c1->cd(1);
  //avec_draw(y1,"Avec","x","y","AL",1,"avec_draw_test_1");
  avec_draw(y1,"Avec","x");

  // Single Avec against Avec
  c1->cd(2);
  //avec_draw(x1, y1,"Avec \% Avec","x","y","AL",1,"avec_draw_test_2"); // All parameters passed
  avec_draw(x1, y1); //Shortest possible version 

  // Single Avec with Errors against Avec
  c1->cd(3);
  avec_draw(x3, y4, y5, err4, "Avec \% Avec with Errors","x","y","AP");


  std::cout << "Testing TMultiGraphs" << endl;
  // Test for TMultiGraph
  TCanvas* c2=new TCanvas("cv_mgr","multigraph test");
  c2->Divide(2,2);

  // Draw set of Avecs against one Avec
  c2->cd(1);
  avec_draw(x1,(y1 | y2),"Set of Avecs \% Avec","XTitle","YTitle");

  // Draw set of Avecs against set of Avecs
  c2->cd(2);
  avec_draw((x1 | x2),(y1 | y3),"Set of Avecs \% Set of Avecs","XTitle","YTitle","AL");

  // Draw set of Avecs with errors against set of Avecs
  c2->cd(3);
  TMultiGraph* mgr = avec_draw(x3, (y4 | y5),(err4 | err5),"Set of Avecs with errors \% Avec","XTitle","YTitle","AP");

  // Draw columns in Avec2D against one Avec
  c2->cd(4);
  Avec2D y2D;
  y2D.push_back(y1);
  y2D.push_back(y2);
  avec_draw(x1,y2D,"Avec2D \% Avec","XTitle","YTitle");
  
  
  TCanvas* c2=new TCanvas("cv_mgr_2","Test MultiGraph Palette");

  Avec px(100, 0, TMath::TwoPi());
  std::vector<Avec> vy;
  for(double coeff = 1.0; coeff <= 3.0; coeff += 0.1) vy.push_back( coeff * sin(px)); 
  Avec py1 = sin(px);
  Avec py2 = 1.1 * sin(px);
  Avec py3 = 1.2 * sin(px);
  Avec py4 = 1.3 * sin(px);
  Avec py5 = 1.4 * sin(px);
  Avec py6 = 1.5 * sin(px);
  Avec py7 = 1.6 * sin(px);
  Avec py8 = 1.7 * sin(px);
  Avec py9 = 1.8 * sin(px);
  Avec py10 = 1.9 * sin(px);
  //avec_draw(px, (py1|py2|py3|py4|py5|py6|py7|py8|py9|py10), "", "", "", "ALP");
  avec_draw(px, vy, "", "x", "y", "ALP");


/*   TLegend* leg = new TLegend(0.4,0.5,0.89,0.89); */
/*   //leg->SetTextFont(72); */
/*   leg->SetTextSize(0.08); */
/*   leg->SetHeader("Avec Legend"); */
/*   leg->AddEntry(mgr->GetListOfGraphs()->At(0),"Avec 1","P"); */
/*   leg->AddEntry(mgr->GetListOfGraphs()->At(1),"Avec 2","P"); */
/*   leg->Draw(); */

}
