{
   gSystem -> Load ( "Avec.so" );

   gStyle -> SetCanvasBorderMode ( 0 );
   gStyle -> SetPadBorderMode ( 0 );
   gStyle -> SetCanvasColor ( 0 );
   gStyle->SetTitleXSize(0.05);
   gStyle->SetTitleYSize(0.05);
   gStyle->SetLabelSize(0.05);
   gStyle -> SetOptFit ( 1011 );

   cout << "rootlogon loaded" << endl;

   Avec::VERBOSE_FLAG = 1;

   cout << "Testing avec_draw" << endl;
   gROOT -> LoadMacro ( "examples/avec_draw_test.C" );
   avec_draw_test_tmultigraph();

   cout << "Testing TObject inheritance" << endl;
   gROOT -> LoadMacro ( "examples/avec_tobject_test.C" );
   avec_tobject_test();
   avec2d_tobject_test();

   cout << "Testing Avec2D" << endl;
   gROOT -> LoadMacro ( "examples/avec2d_test.C" );
   avec2d_test();

   cout << "Testing avec_plot" << endl;
   gROOT -> LoadMacro ( "examples/avec_plot_test.C" );
   avec_plot_test();

}
