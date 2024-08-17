#include "RhocpApp.h"
#include "MooseMain.h"

// #include "MooseInit.h"
// #include "Moose.h"
// #include "MooseApp.h"
// #include "AppFactory.h"

int
main(int argc, char * argv[])
{
//   // Initialize MPI, solvers and MOOSE
//   MooseInit init(argc, argv);
//
//   // Register this application's MooseApp and any it depends on
//   RhocpApp::registerApps();
//
//   // Create an instance of the application and store it in a smart pointer for easy cleanup
//   std::shared_ptr<MooseApp> app = AppFactory::createAppShared("RhocpApp", argc, argv);
//
//   app->setErrorOverridden();
//
//   // Execute the application
//   app->run();

  Moose::main<RhocpApp>(argc, argv);

  return 0;
}
