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

## **Compilation instructions**

**Requirements** 

* A Bourne-like shell
* bzr
* git
* curl
* patch
* emsdk (Steps 1 and 2 below)
* cmake
* A *native* C++ compiler. (This is required for part of the `boost` build system).
* `autoconf`,`autotools`
* `libtool`

1. Install emscripten (following  [https://emscripten.org/docs/getting_started/downloads.html](https://emscripten.org/docs/getting_started/downloads.html)):  
`git clone https://github.com/emscripten-core/emsdk.git`  
`cd emsdk`  
`git pull`  
`./emsdk install latest`  
`./emsdk activate latest`

2. Each time you want to use emscripten:  
`source ./emsdk_env.sh`

3. Get the sources:  
`git clone https://github.com/glycojones/privateer.git privateer_wasm`  
`cd privateer_wasm`  
`git checkout wasm`
`./get_sources`

4. Build required libraries and privateer WASM module using:  
`emcmake cmake .`  
`emmake make .`

5. To run the Privateer application use: 
`cd webserver`
`npm run dev`

## Releases
- V0.1 16.8.23 - Full UI Redesign with ID correction - Initial Release
