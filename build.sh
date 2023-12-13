source /Applications/ccp4-8.0/bin/ccp4.setup-sh 
source ~/Development/privateer/emsdk/emsdk_env.sh
# source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh 
# source ~/dev/privateer_wasm/emsdk/emsdk_env.sh

emcmake cmake .
emmake make -j 

# cp privateer.js webapp/api/privateer.js
# cp privateer.wasm webapp/api/privateer.wasm
# cp privateer.data webapp/api/privateer.data

mv privateer.js webapp/src/wasm/privateer.js
mv privateer.wasm webapp/src/wasm/privateer.wasm
cp privateer.data webapp/public/privateer.data


# mv privateer.js web/privateer/privateer/src/wasm/privateer.js
# mv privateer.wasm web/privateer/privateer/src/wasm/privateer.wasm
# mv privateer.data web/privateer/privateer/src/wasm/privateer.data

# cp web/privateer/privateer/src/wasm/privateer.data web/privateer/privateer/public/privateer.data
# cp web/privateer/privateer/src/wasm/privateer.wasm web/privateer/privateer/public/privateer.wasm
# cp web/privateer/privateer/src/wasm/privateer.js web/privateer/privateer/public/privateer.js
