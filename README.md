![ScreenShot](/logo.png)

# Privateer
Privateer is a tool for carbohydrate structure validation, re-refinement and graphical analysis. It carries out automatic assignments of ring conformation (IUPAC nomenclature), anomeric form, absolute configuration and comparison to reference values for validation. It computes omit mFo-DFc maps and calculates a correlation coefficient between model and electron density. For structure refinement, it is able to generate chemical dictionaries with unimodal torsion restraints which will help keep the lowest energy conformation. In terms of graphical analysis, it will produce vector diagrams in SNFG nomenclature (SVG format), which are annotated using the validation information (ring conformation, anomeric form, etc).

Privateer is implemented in C++11 and Python3 (via pybind11), and produces Scheme and Python scripts for use with Coot (https://github.com/pemsley/coot).

Skip to sections: [Build it yourself!](#build-it-yourself) – [Step by step build](#step-by-step-build) - [Validating a model programmatically](#validating-a-model-programmatically) 

## Authors
Jon Agirre (@glycojones), with contributions from Haroldas Bagdonas (@GABRAH), Jordan Dialpuri (@Dialpuri) and Lucy Schofield (@lcs551), York Structural Biology Laboratory (YSBL), Department of Chemistry, The University of York, UK. We have implemented great suggestions and feedback from Gideon Davies FRS, Wendy Offen, Christian Roth, Marcin Wojdyr, Robbie Joosten, Paul Emsley, Garib Murshudov, Eleanor Dodson FRS, Saioa Urresti, Jacob Sorensen and Stuart McNicholas. 

## Funding 
* 2013-2015 BBSRC grant awarded to Keith Wilson (York)
* 2015-2017 BBSRC grant awarded to CCP4/Kevin Cowtan (York)
* 2017-2025 Royal Society University Research Fellowship awarded to Jon Agirre (York)
* 2018-2021 EPSRC DTP studentship allocated to Jon Agirre & awarded to Mihaela Atanasova (York)
* 2019-2023 Royal Society Research Grant awarded to Jon Agirre (York), covering 4-year studentship awarded to Haroldas Bagdonas

## Citations
If you find Privateer useful, please cite:

- Agirre, J., Iglesias-Fernández, J., Rovira, C., Davies, G.J., Wilson, K.S. and Cowtan, K.D., 2015. Privateer: software for the conformational validation of carbohydrate structures. Nature Structural and Molecular Biology, 22(11), p.833.
- Bagdonas, H., Ungar, D. and Agirre, J., 2020. Leveraging glycomics data in glycoprotein 3D structure validation with Privateer. Beilstein Journal of Organic Chemistry, 16(1), 2523-2533.
- Agirre, J., 2017. Strategies for carbohydrate model building, refinement and validation. Acta Crystallographica Section D, 73(2), pp.171-186.



## **Build it yourself!**
Stable versions of Privateer are routinely released by CCP4 (https://www.ccp4.ac.uk) and CCP-EM (https://www.ccpem.ac.uk). While we recommend that users download and install those suites for X-ray crystallography and Electron cryo-microscopy work respectively, we recognise that they do not always come with the latest functionality. If you would like to test the latest versions (we'd love you to), then you'll need to build Privateer yourself. It is not as hard as it sounds, thanks in part to the work Haroldas Bagdonas (@GABRAH) put into a build script that gathers and builds all dependencies.  

#### Operating systems supported
MacOS (tested up to Ventura) and GNU/Linux (tested on Ubuntu 22.04). We do not test our builds on Windows, at least not routinely. 

#### Prerequisites

`bzr` (breezy on Mac homebrew)
> This is CCP4's version control system, and it will be required until CCP4 starts using git.

`cmake` (minimum version required 3.12)
> Privateer itself requires cmake.

`python3 virtualenv coreutils wget` 
> Python3 and virtualenv will allow you to have you own little Python with access to the Privateer module (import privateer). We need coreutils and wget so the build scripts can function.

`gcc g++ gfortran m4` 
> These are the compilers, for C/C++/Fortran. On Debian/Ubuntu: `apt-get install build-essential m4 gfortran`. On MacOS (via Homebrew): `brew install gcc gfortran m4`


#### Step by step build 
First, select your preferred Python (up to 3.10 suppported) using `pyenv install [version]` and then `pyenv global [version]`
This will give you that particular version of Python in your virtualenv.

```
git clone https://github.com/glycojones/privateer.git
```
> Creates a copy of the repository on your machine and switches it to the *master* branch.

```
cd privateer
```
> Goes into the directory with the repository, which has been created by the previous command.

```
git submodule update --init --recursive
```
> This will fetch those dependencies that are managed by git.

```
virtualenv privateerpython
```
> It creates your own Python environment within the repository.

```
source ccp4.envsetup-sh
```
> CCP4 environment variables, required for Privateer to work.

```
source privateerpython/bin/activate
```
> Will setup the Python environment.

```
pip install -r requirements.txt
```
> Installs the list of Python dependencies.

```
python setup.py install
```
> This is the one command that actually builds something. It will take a while.


## Validating a model programmatically

Once built, the executable can be found in `build/executable/privateer`

Before executing Privateer it is important to have ccp4.envsetup-sh properly sourced

```
.././privateer -pdbin ../../../tests/test_data/5fjj.pdb -mtzin ../../../tests/test_data/5fjj.mtz -glytoucan -databasein ../../../src/privateer/privateer_database.json -debug_output -all_permutations`
```

Running the binary is not the only way of accessing Privateer. Our fully featured Python module, built on top of the C++ code, can also be accessed via the python virtualenv you should have created as part of the build process:

    from privateer import privateer_core as pvt
    print(dir(pvt))
    wurcs = pvt.print_wurcs("/Users/haroldas/Dev/privateer_dev_noccp4/tests/test_data/2h6o_carbremediation.pdb")
    print(wurcs)
