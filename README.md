# HBI
Multi-dimensional earthquake cycle simulation code based on Boundary Element Method with H-matrices.

Features

- 2D and 3D planar/nonplanar faults in a full/half space

- Rectangular and triangular meshes for 3D problems (require .stl mesh file for triangular mesh)

- Rate-State friction law

- Quasi-dynamic approximation using radiation-damping term

- Hybrid MPI/open-MP parallelization

See Documentation_for_HBI.pdf for more information.

This software is freely available under the MIT license.
If you write a paper using this code, please cite the following paper:

So Ozawa, Akihiro Ida, Tetsuya Hoshino, Ryosuke Ando (2023),
"Large-scale earthquake sequence simulations of 3D geometrically complex faults using the boundary element method accelerated by lattice H-matrices", Geophysical Journal International,232 (3), 1471-1481 https://doi.org/10.1093/gji/ggac386
