# Vortices in the magnetically-induced current density in LiH molecule studied through the lens of the Omega function


| ![Omega_bz.png](screenshots/LiH_MICD/repImageGray.png){width=500} | 
|:-:|
|<div style="width:500px"><b>An automatic approach based on Topological Data Analysis extracts axial (blue) and toroidal (green) vortices in magnetically-induced current density as specific sub-sets of the separatrices (gray curves) of the Morse-Smale complex of the Omega index.</b></div>|

## Pipeline description

This example illustrates (1) the calculation of the magnetically-induced current density (MICD) tensor in the LiH molecule in the `DIRAC` software, followed by (2) the calculation of the Omega index with the `qcten` script and (3) its subsequent topological analysis in the `TTK` software.



The first step involving quantum chemistry calculations aims to compute the MICD tensor and its gradient and export them on a 3D grid.



The purpose of the second step is a pointwise derivation of a scalar function from these tensor fields. In this case, the studied scalar field represents the so-called Omega index, used as an indicator of vortices in the first-order current density field. This step also involves translating data exported from `DIRAC` in TXT to the VTI format favored by the `TTK` code. Simultaneously, it applies the resampling filter ("ResampleToImage") without changing the number of grid points or grid bounds.



The final step involves analyzing the Omega scalar field in the `TTK` software. It starts with extracting all critical point pairs, determining a persistence threshold for the salient pairs, and using this threshold to simplify the topology of the Omega scalar field. Then, the computation of the Morse-Smale complex of such simplified field results, among others, in the extraction of its one-dimensional separatrices. A subset of these separatrices connecting 2-saddles and maxima well captures the center lines of vortices. These can then be associated with axial and toroidal vortices in the MICD of the LiH molecule. Notable `TTK` filters employed in this analysis are the `PersistenceDiagram`, `TopologicalSimplification`, `MorseSmaleComplex`, and `PersistentGenerators`.



For more details on this analysis, please see the publication [on arXiv](https://arxiv.org/abs/2212.08690). All data files generated in this analysis are on [zenodo](https://zenodo.org/record/7446735#.Y8BlkNKE4XU). All links are [at the bottom of this page](#resources-and-additional-information).

Below, we describe these three steps in more detail.



## Quantum chemistry calculations

### Setup

The experimental geometry of LiH molecule from the [NIST database](https://cccbdb.nist.gov/exp2x.asp?casno=7580678&charge=0#NISTdiatomic and https://www.nist.gov/pml/diatomic-spectral-database) was used (R(Li-H) = 1.595 Angstrom).

The MICD tensor and its gradient were calculated analytically in the development version of the `DIRAC` software (commit hash `2330f11`) with the Dirac-Coulomb Hamiltonian, the B3LYP exchange-correlation functional, and the def-TZVP basis set applied for both atoms. London atomic orbitals and the simple magnetic balance scheme were applied in response calculations. The densities were exported on the cube grid of 128 points in each Cartesian direction using the default visualization options in `DIRAC`.


### `DIRAC` inputs

* Molecular geometry of LiH molecule in XYZ format (in Angstrom): [LiH.xyz](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/LiH.xyz)
* Input for a wave function optimization: [scf.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/scf.inp)
* Input for calculations of the magnetic-field response (uses NMR shielding calculations): [prp.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/prp.inp)
* Inputs for calculations of the components of the MICD tensor, composed of the elements of the current density vector induced by the magnetic field applied in 
    the "x"-direction ([jbx.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/jbx.inp)), 
    the "y"-direction ([jby.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/jby.inp)), and 
    the "z"-direction ([jbz.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/jbz.inp))
* Inputs for calculations of the components of the gradient of the MICD tensor, composed of the elements of the gradient of the current density vector induced by the magnetic field applied in 
    the "x"-direction ([gradjbx.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/gradjbx.inp)), 
    the "y"-direction ([gradjby.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/gradjby.inp)), and 
    the "z"-direction ([gradjbz.inp](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/inputs/visgrid_cube_128/gradjbz.inp))

### `DIRAC` outputs

* Files with exported elements of the MICD tensor and its gradient on a grid in TXT format; these are also available on [zenodo](https://zenodo.org/record/7446735#.Y8BlkNKE4XU).
* `DIRAC` text output files, available in [the repository](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/dirac/dc_b3lyp_def2tzvp/outputs) and on [zenodo](https://zenodo.org/record/7446735#.Y8BlkNKE4XU).

### Execution

* Below, we assume that the `pam` script of `DIRAC` is available in `$PATH`.

* Step 1. Wave function optimization:

```
mol=LiH.xyz
inp_scf=scf.inp
pam --inp=$inp_scf --mol=$mol --outcmo
```

* Step 2. Response calculations:

```
mol=LiH.xyz
inp_prp=prp.inp
pam --inp=$inp_prp --mol=$mol --incmo --get="DFCOEF=DFCOEF.smb TBMO PAMXVC"
```


* Step 3. Calculations and export of real-space densities,

    * involving the elements of the MICD tensor, on an example of the `jbx.inp` file:

      ```
      mol=LiH.xyz
      vis=jbx
      inp_vis=$vis.inp
      pam --inp=$inp_vis --mol=$mol --put="DFCOEF.smb=DFCOEF TBMO PAMXVC" --get="plot.3d.vector=$vis.txt"
      ```

    * involving the elements of the gradient of the MICD tensor, on an example of the `gradjbx.inp` file:

      ```
      mol=LiH.xyz
      vis=gradjbx
      inp_vis=$vis.inp
      pam --inp=$inp_vis --mol=$mol --put="DFCOEF.smb=DFCOEF TBMO PAMXVC" --get="plot.3d.tensor=$vis.txt"
      ```

    * analogous computations should be done with the `jby.inp` file (please change of `vis=jbx` to `vis=jby`) and the `jbz.inp` file (please change `vis=jbx` to `vis=jby`),
    * analogous computations should be done with the `gradjby.inp` file (please change of `vis=gradjbx` to `vis=gradjby`) and the `gradjbz.inp` file (please change `vis=gradjbx` to `vis=gradjby`).


## Calculation of scalar functions for the topological data analysis

### Inputs

* Data exported from `DIRAC`:
    * The "x"/"y"/"z" elements of the current density vector field induced by the magnetic field applied in the "x" direction are in the 4th/5th/6th-column of the `jbx.txt` file.
    * The "x"/"y"/"z" elements of the current density vector field induced by the magnetic field applied in the "y" direction are in the 4th/5th/6th-column of the `jby.txt` file.
    * The "x"/"y"/"z" elements of the current density vector field induced by the magnetic field applied in the "z" direction are in the 4th/5th/6th-column of the `jbz.txt` file.
    * The "x"/"y"/"z" elements of the gradient of the current density vector field induced by the magnetic field applied in the "z" direction are in the `gradjbx.txt` file, starting from the 4th column).
    * In all TXT files exported from `DIRAC`, first three columns refer to the "x"/"y"/"z"-coordinates of grid vertices.

### Outputs
### Run script


## Topological Data Analysis

### Inputs

* It may help to mark the positions of the Li and H nuclei on the plots; for this purpose, the molecular geometry of the LiH molecule in CSV format (in atomic units) is available in the [LiH.csv](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/LiH.csv) file.

* MICD-related data in VTI format:

    * Omega function derived from the magnetically-induced current density vector corresponding to the perturbation of the magnetic field applied perpendicularly to the Li-H bond: `start_data_omega_bz.vti` file in [the repository](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/vti/start_data_omega_bz.vti) and on [zenodo](https://zenodo.org/record/7446735#.Y8E2dtKZNhF); data description:
        * `omega_bz` - corresponds to Omega function calculated for the magnetic field applied perpendicularly to the Li-H bond ("bz");
        * `bz_wz` - corresponds to the "z"-component of the curl of the current density vector induced by the magnetic field applied perpendicularly to the Li-H bond ("bz"); it is a zz-component of the vorticity tensor.

    * The elements of the magnetically-induced current density vector corresponding to the perturbation of the magnetic field applied perpendicularly to the Li-H bond are on the `start_data_bz.vti` file in [the repository](https://github.com/tda-qchem/tda-qchem-examples/tree/main/data/LiH_MICD/vti/start_data_bz.vti) and on [zenodo](https://zenodo.org/record/7446735#.Y8E2dtKZNhF); the "x"/"y"/"z" elements of this vector are marked as  `bz_jx`, `bz_jy`, `bz_jz`, respectively.

    * Additionally, on [zenodo](https://zenodo.org/record/7446735#.Y8E2dtKZNhF), we also share the VTI file which contains all the elements of the full MICD tensor (`startdatajbtensor.vti` file).
    

### Outputs

* The `Paraview` state files and all other files enabling the reproduction of the above screenshot and all images attached to the companion publication are in [the repository](https://github.com/tda-qchem/tda-qchem-examples/tree/main/pvsm).


### ParaView

To reproduce the images and to explore the TDA pipeline, go to the root directory of [this repository](https://github.com/tda-qchem/tda-qchem-examples) and enter the following command:

``` bash
paraview --state=pvsm/lih.pvsm
```

## Resources and additional information

* [Prerequisites](https://tda-qchem.github.io/tda-qchem-examples/)
* [DIRAC](http://www.diracprogram.org/)
* [TTK](https://topology-tool-kit.github.io/)
* [qcten]()

* The calculations and export of the magnetically-induced current density are also discusssed in [the official DIRAC tutorial](http://www.diracprogram.org/doc/release-22/tutorials/visual/general/tutorial.html#densities-and-currents-induced-by-a-magnetic-field).

* Related data: the publication on [arXiv](https://arxiv.org/abs/2212.08690) and its [1-page summary](), data files generated in this analysis on [zenodo](https://zenodo.org/record/7446735#.Y8BlkNKE4XU).

* To fully reproduce the results reported in the publication, please check [this link]().
* To fully reproduce the publication, please check [this link]().
