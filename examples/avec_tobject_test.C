void avec_tobject_test()
{
   Avec y1(10,1,100);
   Avec y2=y1*y1;
   y1.SetName("y1");
   y2.SetName("y2");
   
   cout << "ClassName: " << y1.ClassName() << endl;
   cout << "Dump: " << y1.Dump() << endl; // lists all data members and their current valsue
   //cout << "Inspect: " << y1.Inspect() << endl; // opens a window to browser data members at all levels
   //cout << "DrawClass: " << y1.DrawClass() << endl; // Draws the class inheritance tree

   //TClass *cl = y1->IsA();
   TFile f("output.root","RECREATE");
   y1.Write("y1");
   y2.Write("y2");
   f.GetListOfKeys()->Print();
   f.Close();
   TFile f2("output.root","READ");
   f2.GetListOfKeys()->Print();
   Avec* v1=(Avec*)f2.Get("y2");
   TCanvas* c2=new TCanvas("c2","c2");
   avec_draw(*v1);
   f2.Close();
}

void avec2d_tobject_test()
{
   Avec2D x2d(20,Avec(20,0,7));

  	Avec2D y2d=sin(x2d)*cos(vtrans(x2d));
	y2d.SetName("y2d");

   cout << "ClassName: " << y2d.ClassName() << endl;
   cout << "Dump: " << y2d.Dump() << endl; // lists all data members and their current valsue
   //cout << "Inspect: " << y2d.Inspect() << endl; // opens a window to browser data members at all levels
   //cout << "DrawClass: " << y2d.DrawClass() << endl; // Draws the class inheritance tree

   //TClass *cl = y2d->IsA();
   TFile f("output.root","RECREATE");
   y2d.Write();
   f.GetListOfKeys()->Print();
   f.Close();
   TFile f2("output.root","READ");
   f2.GetListOfKeys()->Print();
   Avec2D* v2d=(Avec2D*)f2.Get("y2d");
   cout << "v2d->size(): " << v2d->size() << endl;
   TCanvas* c3=new TCanvas("c3","c3");
   //c3->cd(1);
   avec_draw(*v2d);
   f2.Close();
}
