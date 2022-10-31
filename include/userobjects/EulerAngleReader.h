#ifndef EULERANGLEREADER_H
#define EULERANGLEREADER_H

#include "GeneralUserObject.h"

// Forward declaration
class EulerAngleReader;

template <>
InputParameters validParams<EulerAngleReader>();

class EulerAngleReader : public GeneralUserObject
{
public:
  EulerAngleReader(const InputParameters & parameters);
  
  virtual ~EulerAngleReader() {}

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}
  
  Real getData(unsigned int, unsigned int) const;

protected:
  FileName _file_name;
  std::vector<Real> _data;
  unsigned int _ngrain;
  
  void readFile();
};

#endif // EULERANGLEREADER_H