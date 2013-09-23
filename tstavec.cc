//--------------------------------------------------
#include "TROOT.h"
#include "TRint.h"


int main(int argc, char **argv)
{
  // Create an interactive ROOT interface
  TRint *theApp = new TRint("ROOT session", &argc, argv, NULL, 0);

  // Run interactive interface
  theApp->Run();

    return(0);
}
//--------------------------------------------------
