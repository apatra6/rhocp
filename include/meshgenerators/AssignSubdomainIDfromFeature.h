#ifndef ASSIGNSUBDOMAINIDFROMFEATURE_H
#define ASSIGNSUBDOMAINIDFROMFEATURE_H

// MOOSE includes
#include "MeshGenerator.h"
#include "EBSDAccessFunctors.h"

/**
 * MeshGenerator for assigning a subdomain ID to all elements
 */
class AssignSubdomainIDfromFeature : public MeshGenerator, public EBSDAccessFunctors
{
public:
  static InputParameters validParams();

  AssignSubdomainIDfromFeature(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> & _input;

  virtual void readEBSDfile();

  FileName _EBSDFilename;

  /// Logically three-dimensional data indexed by geometric points in a 1D vector
  std::vector<EBSDPointData> _data;

  /// Dimension of the problem domain
  unsigned int _mesh_dimension;

  /// The number of values in the x, y and z directions.
  unsigned _nx, _ny, _nz;

  /// The spacing of the values in x, y and z directions.
  Real _dx, _dy, _dz;

  /// Grid origin
  Real _minx, _miny, _minz;

  /// Maximum grid extent
  Real _maxx, _maxy, _maxz;

  /// The subdomain ID to assign to every elemennt
  SubdomainID _subdomain_id;

    /**
   * Get the requested type of data at the point p.
   */
  const EBSDPointData & getData(const Point & p) const;

  /// Computes a global index in the _data array given an input *centroid* point
  unsigned indexFromPoint(const Point & p) const;
};

#endif // ASSIGNSUBDOMAINIDFROMPHASE_H
