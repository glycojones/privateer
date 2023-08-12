source /Applications/ccp4-8.0/bin/ccp4.setup-sh 
source ~/Development/privateer/emsdk/emsdk_env.sh
# source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh 
# source ~/dev/privateer_wasm/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j 



mv privateer.js webserver/src/wasm/privateer.js
mv privateer.wasm webserver/src/wasm/privateer.wasm
cp privateer.data webserver/public/privateer.data


# mv privateer.js web/privateer/privateer/src/wasm/privateer.js
# mv privateer.wasm web/privateer/privateer/src/wasm/privateer.wasm
# mv privateer.data web/privateer/privateer/src/wasm/privateer.data

# cp web/privateer/privateer/src/wasm/privateer.data web/privateer/privateer/public/privateer.data
# cp web/privateer/privateer/src/wasm/privateer.wasm web/privateer/privateer/public/privateer.wasm
# cp web/privateer/privateer/src/wasm/privateer.js web/privateer/privateer/public/privateer.js
