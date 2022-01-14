LAMMPS fix for calculating many-body forces using [image method](https://aip.scitation.org/doi/full/10.1063/1.4962832).
-----------------------------------------------------------------

This code is released as a part of [MICCoM](http://miccom-center.org/index) (Midwest Integrated Center for Computational Materials) to calculate three-body forces using `fix` command in LAMMPS. MICCoM develops and disseminates interoperable computational tools - open source software, data, simulation templates, and validation procedures - that enable simulations and predictions of properties of functional materials for energy conversion and of solid-state materials for quantum information science.


src
--------
- `fix_colloidImage.h` and `fix_colloidImage.cpp` : LAMMPS fix for computing three-body forces using image method (Eq. 29 in [Qin et al.](https://aip.scitation.org/doi/full/10.1063/1.4962832))

Installation
------------
- Copy `src/fix_colloidImage.cpp` and `src/fix_colloidImage.h` into lammps/src/ directory.
- Build [LAMMPS](https://docs.lammps.org/Build_cmake.html).

Syntax
-------
- `fix ID group-ID colloid/image ion-type-start IIII einner XXXX`
	- IIII: start ID of ion type.
	- XXXX: dielectric constant of particle.
- Eg: `fix fix_colloid all colloid/image ion-type-start 11 einner 100`
	- See example LAMMPS input script in `examples` for performing ten particle simulations. 

Contributors
------------
- Original author: Jiyuan Li ([dePablo group](https://pme.uchicago.edu/group/de-pablo-group), PME, UChicago)
- Contributors: Siva Dasetty ([Ferguson lab](https://www.ferglab.com), PME, UChicago)

Acknowledgements
----------------
This work was supported by [MICCoM](http://miccom-center.org/index) (Midwest Integrated Center for Computational Materials), as part of the Computational Materials Science Program funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Materials Sciences and Engineering Division.
