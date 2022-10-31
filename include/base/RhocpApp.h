#ifndef RHOCPAPP_H
#define RHOCPAPP_H

#include "MooseApp.h"

class RhocpApp;

template<>
InputParameters validParams<RhocpApp>();

class RhocpApp : public MooseApp
{
public:
  RhocpApp(const InputParameters & parameters);
  virtual ~RhocpApp();

  //static void registerApps();
  //static void registerObjects(Factory & factory);
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
  static void registerExecFlags(Factory & factory);
};

#endif /* CPFEAPP_H */
