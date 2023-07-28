source /Applications/ccp4-8.0/bin/ccp4.setup-sh 
source ~/Development/privateer/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j 

mv privateer.js web/privateer/privateer/wasm/privateer.js
mv privateer.wasm web/privateer/privateer/wasm/privateer.wasm
mv privateer.data web/privateer/privateer/wasm/privateer.data
