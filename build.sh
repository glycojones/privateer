#source /Applications/ccp4-8.0/bin/ccp4.setup-sh
#source ~/Development/privateer/emsdk/emsdk_env.sh
 source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh
 source ~/dev/privateer_wasm/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j

mv privateer.js webapp/src/wasm/privateer.js
mv privateer.wasm webapp/src/wasm/privateer.wasm
cp privateer.data webapp/public/privateer.data
