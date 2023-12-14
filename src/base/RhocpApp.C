#include "RhocpApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
RhocpApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  return params;
}

registerKnownLabel("RhocpApp");

RhocpApp::RhocpApp(const InputParameters & parameters) : MooseApp(parameters)
{
  RhocpApp::registerAll(_factory, _action_factory, _syntax);
}

RhocpApp::~RhocpApp() {}

void
RhocpApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<RhocpApp>(f, af, s);
  Registry::registerObjectsTo(f, {"RhocpApp"});
  Registry::registerActionsTo(af, {"RhocpApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
RhocpApp::registerApps()
{
  registerApp(RhocpApp);
}

void
RhocpApp::registerObjects(Factory & factory)
{
  mooseDeprecated("use registerAll instead of registerObjects");
  Registry::registerObjectsTo(factory, {"RhocpApp"});
}

void
RhocpApp::registerExecFlags(Factory & /*factory*/)
{
  mooseDeprecated("use registerAll instead of registerExecFlags");
}

void
RhocpApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  mooseDeprecated("use registerAll instead of associateSyntax");
  Registry::registerActionsTo(action_factory, {"RhocpApp"});
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
RhocpApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RhocpApp::registerAll(f, af, s);
}
extern "C" void
RhocpApp__registerApps()
{
  RhocpApp::registerApps();
}
