![ScreenShot](/logo.png)

# Privateer Webserver

![Privateer Webserver Backend](https://github.com/glycojones/privateer/actions/workflows/backend.yml/badge.svg) ![Privateer Webserver Frontend](https://github.com/glycojones/privateer/actions/workflows/frontend.yml/badge.svg)

Privateer is a tool for carbohydrate structure validation, re-refinement and graphical analysis. It carries out automatic assignments of ring conformation (IUPAC nomenclature), anomeric form, absolute configuration and comparison to reference values for validation. It computes omit mFo-DFc maps and calculates a correlation coefficient between model and electron density. For structure refinement, it is able to generate chemical dictionaries with unimodal torsion restraints which will help keep the lowest energy conformation. In terms of graphical analysis, it will produce vector diagrams in SNFG nomenclature (SVG format), which are annotated using the validation information (ring conformation, anomeric form, etc).

This web server allows for online validation of carbohydrate-containing macromolecular structures in PDB or mmCIF format. All processing happens clientside using C++ compiled to WebAssembly using emscripten. 

## **Frontend Development Instructions** ##
To just work on the front end of the Privateer web app do the following:
- Clone this repository
- `cd webapp`
- `npm i`
- `npm run dev`

Before ANY commit to this repository the following code must be run:
`npm run prettier`
AND the following tests must be performed: 
- Upload PDB file only
- Upload PDB file and MTZ file
- Fetch 1234 from the PDB
- Fetch 5FJI from the PDB
- Click on table
- Drag GlycanDetailView around
- Check Moorhen loads map and model


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

## **Ubuntu Compilation instructions**
To compile on Ubuntu with precompiled dependency libraries, use the following commands 
```
  git clone https://github.com/emscripten-core/emsdk.git
  cd emsdk
  git pull
  ./emsdk install latest
  ./emsdk activate latest
  cd ..
  git clone https://github.com/Dialpuri/privateer_webserver_dependencies.git 
  cd privateer_webserver_dependencies
  cp -r * ../
  cd .. 
  source ./emsdk/emsdk_env.sh
  emcmake cmake . -D MODE=TESTING
  emmake make -j 
```
