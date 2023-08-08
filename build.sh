source /Applications/ccp4-8.0/bin/ccp4.setup-sh 
source ~/Development/privateer/emsdk/emsdk_env.sh
source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh 
source ~/dev/privateer_wasm/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j 

mv privateer.js web/privateer/privateer/src/wasm/privateer.js
mv privateer.wasm web/privateer/privateer/src/wasm/privateer.wasm
mv privateer.data web/privateer/privateer/src/wasm/privateer.data
