// This code is an adapted and shorter version of moose/modules/phase_field/src/mesh_generators/EBSDMeshGenerator.C
// Only parts related to reading in EBSD data and mapping to elements are retained here
// Known issues: While reading pixel coordinates with decimal values, an error regarding some coordinates not within bounds of the EBSD domain might appear
// Fix: Reduce the number of decimal points or convert the step size and input coordinates to integers

#include "EBSDMeshReader.h"
#include "EBSDMesh.h"
#include "MooseMesh.h"
#include "Conversion.h"
#include "NonlinearSystem.h"
#include <Eigen/Dense>

#include <fstream>

registerMooseObject("RhocpApp", EBSDMeshReader);

InputParameters
EBSDMeshReader::validParams()
{
  InputParameters params = EulerAngleProvider::validParams();
  params.addClassDescription("Load and manage DREAM.3D EBSD data files for running simulations on "
                             "reconstructed microstructures.");
  params.addParam<unsigned int>(
      "custom_columns", 0, "Number of additional custom data columns to read from the EBSD file");
  params.addParam<unsigned int>("bins", 20, "Number of bins to segregate quaternions");
  params.addParam<Real>("L_norm", 1, "Specifies the type of average the user intends to perform");
  params.addParam<std::string>("ebsd_meshgenerator",
                               "Specify the name of the EBSDMeshGenerator. The EBSDReader can "
                               "autodetect this, if only one such MeshGenerator exists.");
  return params;
}

EBSDMeshReader::EBSDMeshReader(const InputParameters & params)
  : EulerAngleProvider(params),
    _mesh(_fe_problem.mesh()),
    _nl(_fe_problem.getNonlinearSystemBase(_sys.number())),
    _grain_num(0),
    _custom_columns(getParam<unsigned int>("custom_columns")),
    _time_step(_fe_problem.timeStep()),
    _mesh_dimension(_mesh.dimension()),
    _bins(getParam<unsigned int>("bins")),
    _L_norm(getParam<Real>("L_norm")),
    _nx(0),
    _ny(0),
    _nz(0),
    _dx(0.),
    _dy(0.),
    _dz(0.)
{
  readFile();
  // throws an error for zero bins
  if (_bins == 0)
    mooseError("One cannot have zero bins");
}

void
EBSDMeshReader::readFile()
{
  // No need to re-read data upon recovery
  if (_app.isRecovering())
    return;

  std::string ebsd_filename;
  EBSDMeshGenerator::Geometry geometry;

  // Fetch and check mesh or meshgenerators
  EBSDMesh * mesh = dynamic_cast<EBSDMesh *>(&_mesh);
  if (mesh != NULL)
  {
    ebsd_filename = mesh->getEBSDFilename();
    geometry = mesh->getEBSDGeometry();
  }
  else
  {
    std::string ebsd_meshgenerator_name;
    if (isParamValid("ebsd_meshgenerator"))
      ebsd_meshgenerator_name = getParam<std::string>("ebsd_meshgenerator");
    else
    {
      auto meshgenerator_names = _app.getMeshGeneratorNames();
      for (auto & mgn : meshgenerator_names)
      {
        const EBSDMeshGenerator * emg =
            dynamic_cast<const EBSDMeshGenerator *>(&_app.getMeshGenerator(mgn));
        if (emg)
        {
          if (!ebsd_meshgenerator_name.empty())
            mooseError("Found multiple EBSDMeshGenerator objects (",
                       ebsd_meshgenerator_name,
                       " and ",
                       mgn,
                       "). Use the 'ebsd_meshgenerator' parameter to specify which one to use.");
          ebsd_meshgenerator_name = mgn;
        }
      }

      if (ebsd_meshgenerator_name.empty())
        mooseError("Failed to autodetect an EBSDMeshGenerator (or a deprecated EBSDMesh object).");
    }

    // get the selected or detected mesh generator
    const EBSDMeshGenerator * emg =
        dynamic_cast<const EBSDMeshGenerator *>(&_app.getMeshGenerator(ebsd_meshgenerator_name));
    if (!emg)
      paramError("ebsd_meshgenerator", "No valid EBSDMeshGenerator object found.");

    ebsd_filename = emg->getEBSDFilename();
    geometry = emg->getEBSDGeometry();
  }

  std::ifstream stream_in(ebsd_filename.c_str());
  if (!stream_in)
    mooseError("Can't open EBSD file: ", ebsd_filename);

  // Copy file header data from the EBSDMesh
  _dx = geometry.d[0];
  _nx = geometry.n[0];
  _minx = geometry.min[0];
  _maxx = _minx + _dx * _nx;

  _dy = geometry.d[1];
  _ny = geometry.n[1];
  _miny = geometry.min[1];
  _maxy = _miny + _dy * _ny;

  _dz = geometry.d[2];
  _nz = geometry.n[2];
  _minz = geometry.min[2];
  _maxz = _minz + _dz * _nz;

  // Resize the _data array
  unsigned total_size = geometry.dim < 3 ? _nx * _ny : _nx * _ny * _nz;
  _data.resize(total_size);

  std::string line;
  while (std::getline(stream_in, line))
  {
    if (line.find("#") != 0)
    {
      // Temporary variables to read in on each line
      EBSDPointData d;
      Real x, y, z;

      std::istringstream iss(line);
      iss >> d._phi1 >> d._Phi >> d._phi2 >> x >> y >> z >> d._feature_id >> d._phase >>
          d._symmetry;

//       // Transform angles to degrees
//       d._phi1 *= 180.0 / libMesh::pi;
//       d._Phi *= 180.0 / libMesh::pi;
//       d._phi2 *= 180.0 / libMesh::pi;

      // Custom columns
      d._custom.resize(_custom_columns);
      for (unsigned int i = 0; i < _custom_columns; ++i)
        if (!(iss >> d._custom[i]))
          mooseError("Unable to read in EBSD custom data column #", i);

      if (x < _minx || y < _miny || x > _maxx || y > _maxy ||
          (geometry.dim == 3 && (z < _minz || z > _maxz)))
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
                   geometry.dim,
                   "\n",
                   line);

      d._p = Point(x, y, z);

      // determine number of grains in the dataset
      if (_global_id_map.find(d._feature_id) == _global_id_map.end())
        _global_id_map[d._feature_id] = _grain_num++;

      unsigned int global_index = indexFromPoint(Point(x, y, z));
      _data[global_index] = d;
    }
  }
  stream_in.close();
}

EBSDMeshReader::~EBSDMeshReader() {}

const EBSDMeshReader::EBSDPointData &
EBSDMeshReader::getData(const Point & p) const
{
  return _data[indexFromPoint(p)];
}

const EBSDMeshReader::EBSDAvgData &
EBSDMeshReader::getAvgData(unsigned int var) const
{
  return _avg_data[indexFromIndex(var)];
}

const EulerAngles &
EBSDMeshReader::getEulerAngles(unsigned int var) const
{
  return _avg_angles[indexFromIndex(var)];
}

const EBSDMeshReader::EBSDAvgData &
EBSDMeshReader::getAvgData(unsigned int phase, unsigned int local_id) const
{
  return _avg_data[indexFromIndex(_global_id[phase][local_id])];
}

unsigned int
EBSDMeshReader::getGrainNum() const
{
  return _grain_num;
}

unsigned int
EBSDMeshReader::getGrainNum(unsigned int phase) const
{
  return _global_id[phase].size();
}

unsigned int
EBSDMeshReader::indexFromPoint(const Point & p) const
{
  // Don't assume an ordering on the input data, use the (x, y,
  // z) values of this centroid to determine the index.
  unsigned int x_index, y_index, z_index, global_index;

  x_index = (unsigned int)((p(0) - _minx) / _dx);
  y_index = (unsigned int)((p(1) - _miny) / _dy);
  if ((p(0) < _minx || p(0) > _maxx || p(1) < _miny || p(1) > _maxy) && _mesh_dimension ==2)
    mooseError("Data points must be on the interior of the mesh elements with mesh dimensions 2. In EBSDMeshReader ", name(), ", coordinates: ", p(0), ", ", p(1), ", ", p(2), ", min_range: ", _minx, ", ", _miny, ", ", _minz, ", max_range: ", _maxx, ", ", _maxy, ", ", _maxz);

  if (_mesh_dimension == 3)
  {
    z_index = (unsigned int)((p(2) - _minz) / _dz);
    global_index = z_index * _ny;
    if (p(2) < _minz || p(2) > _maxz)
      mooseError("Data points must be on the interior of the mesh elements with mesh dimensions 3. In EBSDMeshReader ",
                 name(), ", coordinates: ", p(0), ", ", p(1), ", ", p(2), ", min_range: ", _minx, ", ", _miny, ", ", _minz, ", max_range: ", _maxx, ", ", _maxy, ", ", _maxz);
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

unsigned int
EBSDMeshReader::indexFromIndex(unsigned int var) const
{

  // Transfer the index into the _avg_data array.
  unsigned avg_index = var;

  // Don't access out of range!
  if (avg_index >= _avg_data.size())
    mooseError("Error! Index out of range in EBSDMeshReader::indexFromIndex(), index: ",
               avg_index,
               " size: ",
               _avg_data.size());

  return avg_index;
}

const std::map<dof_id_type, std::vector<Real>> &
EBSDMeshReader::getNodeToGrainWeightMap() const
{
  return _node_to_grain_weight_map;
}

const std::map<dof_id_type, std::vector<Real>> &
EBSDMeshReader::getNodeToPhaseWeightMap() const
{
  return _node_to_phase_weight_map;
}

unsigned int
EBSDMeshReader::getGlobalID(unsigned int feature_id) const
{
  auto it = _global_id_map.find(feature_id);
  if (it == _global_id_map.end())
    mooseError("Invalid Feature ID");
  return it->second;
}

void
EBSDMeshReader::meshChanged()
{
  // maps are only rebuild for use in initial conditions, which happens in time step zero
  if (_time_step == 0)
    buildNodeWeightMaps();
}

void
EBSDMeshReader::buildNodeWeightMaps()
{
  // Import nodeToElemMap from MooseMesh for current node
  // This map consists of the node index followed by a vector of element indices that are associated
  // with that node
  const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map =
      _mesh.nodeToActiveSemilocalElemMap();
  libMesh::MeshBase & mesh = _mesh.getMesh();

  // Loop through each node in mesh and calculate eta values for each grain associated with the node
  for (const auto & node : as_range(mesh.active_nodes_begin(), mesh.active_nodes_end()))
  {
    // Get node_id
    const dof_id_type node_id = node->id();

    // Initialize map entries for current node
    _node_to_grain_weight_map[node_id].assign(getGrainNum(), 0.0);
    _node_to_phase_weight_map[node_id].assign(getPhaseNum(), 0.0);

    // Loop through element indices associated with the current node and record weighted eta value
    // in new map
    const auto & node_to_elem_pair = node_to_elem_map.find(node_id);
    if (node_to_elem_pair != node_to_elem_map.end())
    {
      unsigned int n_elems =
          node_to_elem_pair->second
              .size(); // n_elems can range from 1 to 4 for 2D and 1 to 8 for 3D problems

      for (unsigned int ne = 0; ne < n_elems; ++ne)
      {
        // Current element index
        unsigned int elem_id = (node_to_elem_pair->second)[ne];

        // Retrieve EBSD grain number for the current element index
        const Elem * elem = mesh.elem_ptr(elem_id);
        const EBSDMeshReader::EBSDPointData & d = getData(elem->vertex_average());

        // get the (global) grain ID for the EBSD feature ID
        const unsigned int global_id = getGlobalID(d._feature_id);

        // Calculate eta value and add to map
        _node_to_grain_weight_map[node_id][global_id] += 1.0 / n_elems;
        _node_to_phase_weight_map[node_id][d._phase] += 1.0 / n_elems;
      }
    }
  }
}

MooseSharedPointer<EBSDAccessFunctors::EBSDPointDataFunctor>
EBSDMeshReader::getPointDataAccessFunctor(const MooseEnum & field_name) const
{
  EBSDPointDataFunctor * ret_val = NULL;

  switch (field_name)
  {
    case 0: // phi1
      ret_val = new EBSDPointDataPhi1();
      break;
    case 1: // phi
      ret_val = new EBSDPointDataPhi();
      break;
    case 2: // phi2
      ret_val = new EBSDPointDataPhi2();
      break;
    case 3: // grain
      ret_val = new EBSDPointDataFeatureID();
      break;
    case 4: // phase
      ret_val = new EBSDPointDataPhase();
      break;
    case 5: // symmetry
      ret_val = new EBSDPointDataSymmetry();
      break;
    default:
    {
      // check for custom columns
      for (unsigned int i = 0; i < _custom_columns; ++i)
        if (field_name == "CUSTOM" + Moose::stringify(i))
        {
          ret_val = new EBSDPointDataCustom(i);
          break;
        }
    }
  }

  // If ret_val was not set by any of the above cases, throw an error.
  if (!ret_val)
    mooseError("Error:  Please input supported EBSD_param");

  // If we made it here, wrap up the the ret_val in a
  // MooseSharedPointer and ship it off.
  return MooseSharedPointer<EBSDPointDataFunctor>(ret_val);
}

MooseSharedPointer<EBSDAccessFunctors::EBSDAvgDataFunctor>
EBSDMeshReader::getAvgDataAccessFunctor(const MooseEnum & field_name) const
{
  EBSDAvgDataFunctor * ret_val = NULL;

  switch (field_name)
  {
    case 0: // phi1
      ret_val = new EBSDAvgDataPhi1();
      break;
    case 1: // phi
      ret_val = new EBSDAvgDataPhi();
      break;
    case 2: // phi2
      ret_val = new EBSDAvgDataPhi2();
      break;
    case 3: // phase
      ret_val = new EBSDAvgDataPhase();
      break;
    case 4: // symmetry
      ret_val = new EBSDAvgDataSymmetry();
      break;
    case 5: // local_id
      ret_val = new EBSDAvgDataLocalID();
      break;
    case 6: // feature_id
      ret_val = new EBSDAvgDataFeatureID();
      break;
    default:
    {
      // check for custom columns
      for (unsigned int i = 0; i < _custom_columns; ++i)
        if (field_name == "CUSTOM" + Moose::stringify(i))
        {
          ret_val = new EBSDAvgDataCustom(i);
          break;
        }
    }
  }

  // If ret_val was not set by any of the above cases, throw an error.
  if (!ret_val)
    mooseError("Error:  Please input supported EBSD_param");

  // If we made it here, wrap up the the ret_val in a
  // MooseSharedPointer and ship it off.
  return MooseSharedPointer<EBSDAvgDataFunctor>(ret_val);
}
