#include "AssignSubdomainIDfromPhase.h"
#include "MooseMesh.h"
#include "CastUniquePointer.h"

#include "libmesh/elem.h"

#include <fstream>

registerMooseObject("RhocpApp", AssignSubdomainIDfromPhase);

InputParameters
AssignSubdomainIDfromPhase::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<FileName>("EBSDFilename","Name of the EBSD mesh file");
  return params;
}

AssignSubdomainIDfromPhase::AssignSubdomainIDfromPhase(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _EBSDFilename(getParam<FileName>("EBSDFilename"))
{
  readEBSDfile();
}

std::unique_ptr<MeshBase>
AssignSubdomainIDfromPhase::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  for (auto & elem : mesh->element_ptr_range()){
      const EBSDAccessFunctors::EBSDPointData & d = getData(elem->vertex_average());
      _subdomain_id = d._phase;
      elem->subdomain_id() = _subdomain_id;
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}

void
AssignSubdomainIDfromPhase::readEBSDfile()
{
  std::ifstream stream_in(_EBSDFilename);
  if (!stream_in)
    mooseError("Can't open EBSD file: ", _EBSDFilename);

  // Labels to look for in the header
  std::vector<std::string> labels = {
      "x_step", "x_dim", "y_step", "y_dim", "z_step", "z_dim", "x_min", "y_min", "z_min"};

  // Dimension variables to store once they are found in the header
  // X_step, X_Dim, Y_step, Y_Dim, Z_step, Z_Dim
  // We use Reals even though the Dim values should all be integers...
  std::vector<Real> label_vals(labels.size(), 0.0);

  unsigned int header_flag = 0;

  std::string line;
  while (std::getline(stream_in, line))
  {
    if (line.find("#") == 0)
    {
      // Process lines that start with a comment character (comments and meta data)
      std::transform(line.begin(), line.end(), line.begin(), ::tolower);

      for (unsigned i = 0; i < labels.size(); ++i)
        if (line.find(labels[i]) != std::string::npos)
        {
          std::string dummy;
          std::istringstream iss(line);
          iss >> dummy >> dummy >> label_vals[i];

          // One label per line, break out of for loop over labels
          break;
        }
    }
    else if (line.find("#") != 0)
    {
      // first non comment line marks the end of the header
      if (header_flag == 0){
        // Copy stuff out of the label_vars array into class variables
        _dx = label_vals[0];
        _nx = label_vals[1];
        _minx = label_vals[6];
        _maxx = _minx + _dx * _nx;

        _dy = label_vals[2];
        _ny = label_vals[3];
        _miny = label_vals[7];
        _maxy = _miny + _dy * _ny;

        _dz = label_vals[4];
        _nz = label_vals[5];
        _minz = label_vals[8];
        _maxz = _minz + _dz * _nz;

        // Resize the _data array
        unsigned total_size = _nz < 1 ? _nx * _ny : _nx * _ny * _nz;
        _data.resize(total_size);

        _mesh_dimension = _nz < 1 ? 2 : 3;

        header_flag = 1;
      }
      // Temporary variables to read in on each line
      EBSDAccessFunctors::EBSDPointData d;
      Real x, y, z;

      std::istringstream iss(line);
      iss >> d._phi1 >> d._Phi >> d._phi2 >> x >> y >> z >> d._feature_id >> d._phase >>
          d._symmetry;

      // Transform angles to degrees
      d._phi1 *= 180.0 / libMesh::pi;
      d._Phi *= 180.0 / libMesh::pi;
      d._phi2 *= 180.0 / libMesh::pi;

      if (x < _minx || y < _miny || x > _maxx || y > _maxy ||
          (_mesh_dimension == 3 && (z < _minz || z > _maxz)))
        mooseError("EBSD Data ouside of the domain declared in the header ([",
                   _minx,
                   ':',
                   _maxx,
                   "], [",
                   _miny,
                   ':',
                   _maxy,
                   "], [",
                   _minz,
                   ':',
                   _maxz,
                   "]) dim=",
                   _mesh_dimension,
                   "\n",
                   line);

      d._p = Point(x, y, z);

      unsigned int global_index = indexFromPoint(Point(x, y, z));
      _data[global_index] = d;
    }
  }
  stream_in.close();
}

const EBSDAccessFunctors::EBSDPointData &
AssignSubdomainIDfromPhase::getData(const Point & p) const
{
  return _data[indexFromPoint(p)];
}

unsigned int
AssignSubdomainIDfromPhase::indexFromPoint(const Point & p) const
{
  // Don't assume an ordering on the input data, use the (x, y,
  // z) values of this centroid to determine the index.
  unsigned int x_index, y_index, z_index, global_index;

  x_index = (unsigned int)((p(0) - _minx) / _dx);
  y_index = (unsigned int)((p(1) - _miny) / _dy);

  if (_mesh_dimension == 3)
  {
    z_index = (unsigned int)((p(2) - _minz) / _dz);
    global_index = z_index * _ny;
  }
  else
    global_index = 0;

  // Compute the global index into the _data array.  This stores points
  // in a [z][y][x] ordering.
  global_index = (global_index + y_index) * _nx + x_index;

  // Don't access out of range!
  mooseAssert(global_index < _data.size(),
              "global_index " << global_index << " points out of _data range: " << _data.size());

  return global_index;
}
