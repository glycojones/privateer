![ScreenShot](/logo.png)

# Privateer

Privateer is a tool for carbohydrate structure validation, re-refinement and graphical analysis. It carries out automatic assignments of ring conformation (IUPAC nomenclature), anomeric form, absolute configuration and comparison to reference values for validation. It computes omit mFo-DFc maps and calculates a correlation coefficient between model and electron density. For structure refinement, it is able to generate chemical dictionaries with unimodal torsion restraints which will help keep the lowest energy conformation. In terms of graphical analysis, it will produce vector diagrams in SNFG nomenclature (SVG format), which are annotated using the validation information (ring conformation, anomeric form, etc).

This bundle allows for validation of carbohydrate-containing macromolecular structures within the molecular visualisation program ChimeraX. Core Privateer functionality is implemented in C++11, with wrapper functions in Python3 (via pybind11) and a Qt based user interface for use within ChimeraX.

## **Get Privateer for ChimeraX**
The easiest way to get the Privateer bundle for ChimeraX is to install it from within ChimeraX itself. To do this, open the "Tools" menu and select "More Tools..." to open ChimeraX's Toolshed. From here, search for "Privateer" and install the relevant version for your operating system (currently available on macOS and Linux). You can also perform this from the command line by typing `toolshed install ChimeraX_Privateer`.

Alternatively, you can navigate to ChimeraX's Toolshed in an external browser and download the bundle from there, then install from within ChimeraX, as outlined in the section [Installing the Bundle](#installing-the-bundle).

It is also possible to build the bundle from the source code, if desired, as outlined below.

### Operating systems supported
Only MacOS and Linux are supported. There is no bundle or build system for Windows at present.

### Prerequisites

The bundle builds through ChimeraX's bundle builder, so an installation of ChimeraX is required. The Privateer bundle is compatible with ChimeraX version 1.11.0 or higher.

The bundle also relies on precompiled C++ librairies within the ChimeraX-Clipper bundle, so it is necessary to install this prior to compilation. This can be done from within ChimeraX by opening the "Tools" menu and selecting "More Tools..." to open ChimeraX's Toolshed, and installing the relevation version of ChimeraX-Clipper from there. The Privateer bundle is compatible with ChimeraX-Clipper version 0.26 or higher.


### Building the Bundle
First, clone the repository and switch to the ChimeraX bundle branch:

```
git clone --single-branch --branch chimeraX_bundle https://github.com/glycojones/privateer.git chimeraX_bundle
```
Next, navigate to the directory with the repository, created by the previous command:

```
cd chimeraX_bundle
```
Finally, build the bundle:

```
make
```
If you have previously built the bundle and want to clean up compilation products, you can type:
```
make clean
```
This will produce a python .whl file, created in the directory `dist`.

### Installing the Bundle

To install the bundle within ChimeraX from a .whl file rather than from within ChimeraX's Toolshed, utilise the `toolshed install` command within the ChimeraX command line. This involves typing `toolshed install` followed by the absolute path to the .whl file. If a previous version of the bundle has been installed, you will need to first type `toolshed uninstall ChimeraX_Privateer`.
