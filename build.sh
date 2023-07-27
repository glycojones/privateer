source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh 
source /home/jordan/dev/privateer_wasm/emsdk/emsdk_env.sh

emcmake cmake . 
emmake make -j 

mv privateer.js web/privateer_wasm/wasm/privateer.js
mv privateer.wasm web/privateer_wasm/wasm/privateer.wasm