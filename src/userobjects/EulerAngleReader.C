#include "EulerAngleReader.h"

#include <fstream>

registerMooseObject("RhocpApp", EulerAngleReader);

InputParameters
EulerAngleReader::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Read Euler angles from a file and provide it to other objects.");
  params.addRequiredParam<FileName>("file_name", "Euler angle data file name");
  return params;
}

EulerAngleReader::EulerAngleReader(const InputParameters & params)
  : GeneralUserObject(params), _file_name(getParam<FileName>("file_name"))
{
  readFile();
}

void
EulerAngleReader::readFile()
{
  // Read in Euler angles from _file_name
  MooseUtils::checkFileReadable(_file_name);
  std::ifstream file_prop;
  file_prop.open(_file_name.c_str());

  // Skip first x lines (needed for VPSC texture files)
  /* for (unsigned int i = 0; i < x; ++i)
    file_prop.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); */

  // Read number of grains
  file_prop >> _ngrain;
  _data.resize(3 * _ngrain);

  for (unsigned int i = 0; i < _ngrain; i++)
    for (unsigned int j = 0; j < 3; j++)
      if (!(file_prop >> _data[i * 3 + j]))
        mooseError("Error EulerAngleReader: Premature end of file");

  file_prop.close();
}

Real
EulerAngleReader::getData(unsigned int igrain, unsigned int ieuler) const
{
  return _data[(igrain) * 3 + ieuler];
}
