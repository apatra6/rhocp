#ifndef EBSDMESHREADERPOINTDATAAUX_H
#define EBSDMESHREADERPOINTDATAAUX_H

#include "AuxKernel.h"
#include "EBSDMeshReader.h"

class EBSDMeshReaderPointDataAux : public AuxKernel, EBSDAccessFunctors
{
public:
  static InputParameters validParams();

  EBSDMeshReaderPointDataAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  virtual void precalculateValue();

  const EBSDMeshReader & _ebsd_reader;

  /// MooseEnum that stores the type of data this AuxKernel extracts.
  MooseEnum _data_name;

  /// Accessor functor to fetch the selected data field form the EBSD data point
  MooseSharedPointer<EBSDPointDataFunctor> _val;

  /// precalculated element value
  Real _value;
};

#endif // EBSDMESHREADERPOINTDATAAUX_H
