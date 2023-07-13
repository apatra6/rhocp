#ifndef GRAINAREASIZE_H
#define GRAINAREASIZE_H

#pragma once

#include "GeneralUserObject.h"
#include "EulerAngleReader.h"
#include "EBSDMeshReader.h"

class GrainAreaSize : public GeneralUserObject
{
public:
  GrainAreaSize(const InputParameters & params);

  static InputParameters validParams();

  virtual ~GrainAreaSize() {}

  virtual void updateArea();

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}

  /// Reference to the current simulation mesh
  MooseMesh & _mesh;

  const EBSDMeshReader * _EBSDFileReader;

  Real getGrainArea(unsigned int) const;
  Real getGrainSize(unsigned int) const;

protected:
  unsigned int _total_elems;
  std::vector<Real> _grain_area;
  std::vector<Real> _grain_size;
};

#endif // GRAINAREASIZE_H
