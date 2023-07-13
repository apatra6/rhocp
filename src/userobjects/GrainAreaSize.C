#include "GrainAreaSize.h"
#include "MooseMesh.h"
#include "MooseRandom.h"
#include "libmesh/mesh_base.h"

registerMooseObject("RhocpApp", GrainAreaSize);

InputParameters
GrainAreaSize::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("UO for calculating grain/block radius by assuming a 3D volume as a sphere and a 2D volume as a circle");
  params.addParam<UserObjectName>("EBSDFileReader", "Name of the EBSDReader UO");

  return params;
}

GrainAreaSize::GrainAreaSize(const InputParameters & params)
  : GeneralUserObject(params),
  _mesh(_subproblem.mesh()),
  _EBSDFileReader(isParamValid("EBSDFileReader")
                               ? &getUserObject<EBSDMeshReader>("EBSDFileReader")
                               : NULL)
{
    updateArea();
}

void
GrainAreaSize::updateArea()
{
  //Generate neighbor list for each element

  auto & mesh = _mesh.getMesh();
  const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map = _mesh.nodeToElemMap();
  _total_elems = mesh.n_elem();

  _grain_area.clear();
  _grain_size.clear();

  _grain_area.resize(_total_elems);
  _grain_size.resize(_total_elems);

  unsigned int max_grains = 0;

  for (const auto & elem : mesh.active_element_ptr_range()){

      unsigned int elemid = elem->id();

      if (_EBSDFileReader){
          // current element
          Point p = elem->centroid();
          EBSDAccessFunctors::EBSDPointData data = _EBSDFileReader->getData(p);
          _grain_area[data._feature_id] += elem->volume();

          if (data._feature_id > max_grains) {
              max_grains = data._feature_id;
          }
      }
      else {
          // current element
          _grain_area[elem->subdomain_id()] += elem->volume();

          if (elem->subdomain_id() > max_grains) {
              max_grains = elem->subdomain_id();
          }
      }
  }

  if (_mesh.dimension() == 2){
      for (unsigned int igr = 0; igr<=max_grains; igr++) {
          _grain_size[igr] = sqrt(_grain_area[igr]/(22.0/7.0));
          }
  }

  if (_mesh.dimension() == 3){
      for (unsigned int igr = 0; igr<=max_grains; igr++) {
          _grain_size[igr] = cbrt((3.0/4.0)*_grain_area[igr]/(22.0/7.0));
          }
  }

}

Real
GrainAreaSize::getGrainArea(unsigned int grainid) const
{
  return _grain_area[grainid]; //double check: grainid or (grainid - 1)?
}

Real
GrainAreaSize::getGrainSize(unsigned int grainid) const
{
  return _grain_size[grainid]; //double check: grainid or (grainid - 1)?
}
