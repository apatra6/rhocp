import pathlib

import numpy as np

import meshio

import netCDF4

this_dir = pathlib.Path(__file__).resolve().parent

inpfile = this_dir / "810grains_51840elements.inp"
outfile = this_dir / "810grains_51840elements_meshio.e"

mesh = meshio.abaqus.read(inpfile)

mesh_ref = meshio.Mesh(mesh.points, mesh.cells, cell_sets=mesh.cell_sets)

print("Nodes:",np.sum(mesh_ref.points))

print("Elements:",sum(len(cells.data) for cells in mesh_ref.cells))

print("Elsets:",len(mesh_ref.cell_sets))

meshio.exodus.write(outfile,mesh_ref)
