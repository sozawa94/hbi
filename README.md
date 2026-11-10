# HBI
Multi-dimensional earthquake cycle simulation code based on Boundary Element Method with H-matrices.

Features

- 2D and 3D planar/nonplanar faults in a full/half space

- Rectangular and triangular meshes for 3D problems (require .stl mesh file for triangular mesh)

- Rate and state friction law (aging law, slip law, cut-off velocity model, flash heating, and modified CNS model)

- Fault-zone linear and nonlinear viscous flow

- Quasi-dynamic approximation using radiation-damping term

- Along-fault fluid pressure diffusion with permeability evolution

- Many examples for input and parameter files

- Hybrid MPI/open-MP parallelization

See [Documentation](https://drive.google.com/file/d/142j0Ga_sKLtk8SxxUmfIjNq2blFGf-l2/view?usp=sharing) for more information.

This software is freely available under the MIT license.
If you write a paper using this code, please cite the following paper:

So Ozawa, Akihiro Ida, Tetsuya Hoshino, Ryosuke Ando (2023),
"Large-scale earthquake sequence simulations of 3D geometrically complex faults using the boundary element method accelerated by lattice H-matrices", Geophysical Journal International,232 (3), 1471-1481 https://doi.org/10.1093/gji/ggac386

If you write a paper using this code with fluid pressure evolution, please cite the following paper:

So Ozawa, Yuyun Yang, Eric M. Dunham (2024), "Fault Valve Instability: A mechanism for slow slip events", Journal of Geophysical Research: Solid Earth, 129,  https://doi.org/10.1029/2024JB029165,

Papers from other groups using HBI:

Sun, H., & Becker, T. W. (2026). Single fault earthquake cycle complexity. Geophysical Journal International, 247(1), ggag316.

Liu, Y., Xu, B., Yang, Y., & Huang, F. (2026). Aseismic fault-valve leakage and hydraulic memory during underground gas storage. Geoenergy Science and Engineering, 214795.

https://github.com/user-attachments/assets/d1331cee-73f1-4077-9a82-a148b9e94dba

https://github.com/user-attachments/assets/709b4ca6-b1ca-49b2-9375-bb2ad3c86319



