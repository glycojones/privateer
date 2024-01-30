let CCP4Module,RSRModule,currentTaskName="";function rsrPrint(e){postMessage({messageTag:"output",result:e,taskName:currentTaskName})}let molecules_container=null;importScripts("./web_example.js"),importScripts("./pako.js"),createCCP4Module({print:rsrPrint,printErr:rsrPrint}).then((function(e){CCP4Module=e})).catch((e=>{console.log("CCP4 problem :("),console.log(e)})),importScripts("./moorhen.js");const Lib={locateFile:e=>e,onRuntimeInitialized:()=>{console.log("onRuntimeInitialized")},mainScriptUrlOrBlob:"./moorhen.js",print:rsrPrint,printErr:rsrPrint};createRSRModule(Lib).then((function(e){RSRModule=e,molecules_container=new RSRModule.molecules_container_js,molecules_container.geometry_init_standard(),molecules_container.fill_rotamer_probability_tables(),console.log("##################################################"),console.log(molecules_container),console.log("##################################################")})).catch((e=>{console.log("RSR problem :("),console.log(e)}));let dataObjects={pdbFiles:{},mtzFiles:{},cifFiles:{}},dataObjectsNames={pdbFiles:{},mtzFiles:{},cifFiles:{},ramaInfo:{},bvalInfo:{},mol_cont_idx:{},map_cont_idx:{},densityFitInfo:{},rotamersInfo:{}};function guid(){var e=Date.now();return"xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g,(function(a){var t=(e+16*Math.random())%16|0;return e=Math.floor(e/16),("x"===a?t:3&t|8).toString(16)}))}function updateDataObjectsNames(){postMessage({messageTag:"dataObjectsNames",result:dataObjectsNames})}function downLoadFiles(e){let a=[],t=[],s=[],o=[],n=[];for(let i=0;i<e.length;i++)e[i].endsWith(".mtz")?n.push(e[i]):e[i].endsWith(".pdb")||e[i].endsWith(".ent")?s.push(e[i]):e[i].endsWith(".cif")?a.push(e[i]):e[i].endsWith(".pdb.gz")||e[i].endsWith(".ent.gz")?o.push(e[i]):e[i].endsWith(".cif.gz")?t.push(e[i]):console.log("Unknown file suffix");console.log(t),Promise.all(t.map((function(e){return fetch(e).then((e=>e.blob())).then((a=>[e,a]))}))).then((e=>{for(ires=0;ires<e.length;ires++){const a=e[ires];console.log(a),guid()}}))}function loadFiles(e){let a=[],t=[];for(let s=0;s<e.length;s++)e[s].name.endsWith(".mtz")?t.push(e[s]):(e[s].name.endsWith(".pdb")||e[s].name.endsWith(".ent"))&&a.push(e[s]);Promise.all(a.map((function(e){return e.text().then((a=>[e.name,a]))}))).then((function(e){for(ires=0;ires<e.length;ires++){const a=e[ires],t=guid();CCP4Module.FS_createDataFile(".",t+".pdb",a[1],!0,!0),RSRModule.FS_createDataFile(".",t+".pdb",a[1],!0,!0),dataObjects.pdbFiles[t]={fileName:t+".pdb",originalFileName:a[0],contents:a[1]},dataObjectsNames.pdbFiles[t]={fileName:t+".pdb",originalFileName:a[0]};const s=molecules_container.read_pdb(t+".pdb");dataObjectsNames.mol_cont_idx[t]=s}updateDataObjectsNames()})).catch((function(e){console.log(e)})),Promise.all(t.map((function(e){return e.arrayBuffer().then((a=>[e.name,a]))}))).then((function(e){for(ires=0;ires<e.length;ires++){const a=e[ires],t=guid(),s=new Uint8Array(a[1]);CCP4Module.FS_createDataFile(".",t+".mtz",s,!0,!0),RSRModule.FS_createDataFile(".",t+".mtz",s,!0,!0),dataObjects.mtzFiles[t]={fileName:t+".mtz",originalFileName:a[0],contents:s},dataObjectsNames.mtzFiles[t]={fileName:t+".mtz",originalFileName:a[0]};const o="FC",n="PHIC",i="",r=!1,d=!1,l=molecules_container.read_mtz(t+".mtz",o,n,i,r,d);dataObjectsNames.map_cont_idx[t]=l}updateDataObjectsNames()})).catch((function(e){console.log(e)}))}function getDensityFit(e){console.log(e.data),e.data.jobId,dataObjects.pdbFiles[e.data.pdbinKey].fileName;const a=e.data.chainId,t=(dataObjects.mtzFiles[e.data.hklinKey].fileName,dataObjectsNames.mol_cont_idx[e.data.pdbinKey]),s=dataObjectsNames.map_cont_idx[e.data.hklinKey];console.log(t,s);const o=molecules_container.density_fit_analysis(t,s);console.log(o);const n=o.get_index_for_chain(a);console.log(n),console.log(o.cviv.get(n));const i=o.cviv.get(n).rviv;console.log(i),console.log(i.size());let r=[];for(let e=0;e<i.size();e++){const a=i.get(e).residue_spec.chain_id,t=i.get(e).residue_spec.res_no,s=i.get(e).distortion,o={chainId:a,insCode:i.get(e).residue_spec.ins_code,seqNum:t,restype:"UNK",density_fit:1/s};r.push(o)}postMessage({messageId:e.data.messageId,messageTag:"result",result:r,taskName:currentTaskName})}function getBVals(e){console.log(e.data),e.data.jobId;const a=dataObjects.pdbFiles[e.data.pdbinKey].fileName,t=e.data.chainId,s=RSRModule.getBVals(a,t);let o=[];for(let e=0;e<s.size();e++){const a=s.get(e),t={chainId:a.chainId,insCode:a.insCode,seqNum:a.seqNum,restype:a.restype,bval:a.property};o.push(t)}postMessage({messageId:e.data.messageId,messageTag:"result",result:o,taskName:currentTaskName})}function getRama(e){console.log(e.data),e.data.jobId;const a=dataObjects.pdbFiles[e.data.pdbinKey].fileName,t=e.data.chainId,s=RSRModule.getRamachandranData(a,t);console.log(s);let o=[];for(let e=0;e<s.size();e++){const a=s.get(e),t={chainId:a.chainId,insCode:a.insCode,seqNum:a.seqNum,restype:a.restype,phi:a.phi,psi:a.psi,isOutlier:a.isOutlier,is_pre_pro:a.is_pre_pro};o.push(t)}postMessage({messageId:e.data.messageId,messageTag:"result",result:o,taskName:currentTaskName})}function drawMesh(e,a){molecules_container.count_simple_mesh_vertices(e);const t=e.vertices,s=(t.size(),e.triangles);s.size();let o=[],n=[],i=[],r=[];for(let e=0;e<s.size();e++){const a=s.get(e).point_id;o.push(...a)}for(let e=0;e<t.size();e++){const a=t.get(e);n.push(...a.pos),i.push(...a.normal),r.push(...a.color)}const d={prim_types:[["TRIANGLES"]],idx_tri:[[o]],vert_tri:[[n]],norm_tri:[[i]],col_tri:[[r]]};postMessage({messageId:a.data.messageId,messageTag:"result",result:d,taskName:currentTaskName})}function drawDodos(e){const a=dataObjectsNames.mol_cont_idx[e.data.pdbinKey];drawMesh(molecules_container.get_rotamer_dodecs(a),e)}function drawRamaBalls(e){const a=dataObjectsNames.mol_cont_idx[e.data.pdbinKey];drawMesh(molecules_container.ramachandran_validation_markup_mesh(a),e)}function drawCube(e){drawMesh(molecules_container.test_origin_cube(),e)}function getRotamers(e){e.data.jobId;const a=dataObjects.pdbFiles[e.data.pdbinKey].fileName,t=e.data.chainId,s=RSRModule.getRotamersMap();let o=[];const n=RSRModule.getResidueListForChain(a,t),i=RSRModule.getResidueSpecListForChain(a,t);for(let e=0;e<n.size();e++){const a=s.get(n.get(e));if(a){let s=0,r=[];for(s=0;s<a.size();s++){const e=a.get(s),t=e.get_chi(1),o=e.get_chi(2),n=e.get_chi(3),i=e.get_chi(4);r.push([t,o,n,i])}o.push({chainId:t,seqNum:i.get(e).res_no,insCode:i.get(e).ins_code,restype:n.get(e),data:r})}else o.push({chainId:t,seqNum:i.get(e).res_no,insCode:i.get(e).ins_code,restype:n.get(e),data:[]})}postMessage({messageId:e.data.messageId,messageTag:"result",result:o,taskName:currentTaskName})}function autoFitRotamer(e){const a=e.data.jobId,t=(dataObjects.pdbFiles[e.data.pdbinKey].fileName,e.data.chainId),s=e.data.resno,o=a+"out.pdb",n=dataObjectsNames.mol_cont_idx[e.data.pdbinKey],i=dataObjectsNames.map_cont_idx[e.data.hklinKey],r=molecules_container.auto_fit_rotamer(n,t,s,"","",i),d=molecules_container.writePDBASCII(dataObjectsNames.mol_cont_idx[e.data.pdbinKey],o);console.log("result of which is",r,d);var l=RSRModule.FS.readFile(o,{encoding:"utf8"});postMessage({messageId:e.data.messageId,messageTag:"result",result:r,taskName:currentTaskName}),postMessage({messageId:e.data.messageId,messageTag:"pdb_out",jobId:a,result:l,taskName:currentTaskName})}function flipPeptide(e){const a=e.data.jobId,t=(dataObjects.pdbFiles[e.data.pdbinKey].fileName,e.data.chainId),s=e.data.resnoFlip,o=a+"out.pdb",n=new RSRModule.residue_spec_t(t,s,"");console.log("Should in fact call molecules_container.flipPeptide_rs with",dataObjectsNames.mol_cont_idx[e.data.pdbinKey]);const i=molecules_container.flipPeptide_rs(dataObjectsNames.mol_cont_idx[e.data.pdbinKey],n,""),r=molecules_container.writePDBASCII(dataObjectsNames.mol_cont_idx[e.data.pdbinKey],o);console.log("result of which is",i,r);var d=RSRModule.FS.readFile(o,{encoding:"utf8"});postMessage({messageId:e.data.messageId,messageTag:"result",result:i,taskName:currentTaskName}),postMessage({messageId:e.data.messageId,messageTag:"pdb_out",jobId:a,result:d,taskName:currentTaskName})}function miniRSR(e){console.log(e.data);const a=e.data.jobId;var t=new RSRModule.VectorString;t.push_back("mini_rsr"),t.push_back("--pdbin"),t.push_back(dataObjects.pdbFiles[e.data.pdbinKey].fileName),t.push_back("--hklin"),t.push_back(dataObjects.mtzFiles[e.data.hklinKey].fileName),t.push_back("--pdbout"),t.push_back(a+"out.pdb"),t.push_back("--resno-start"),t.push_back(e.data.resnoStart.toString()),t.push_back("--resno-end"),t.push_back(e.data.resnoEnd.toString()),t.push_back("--chain-id"),t.push_back(e.data.chainId),t.push_back("--debug"),console.log("Calling");var s=RSRModule.mini_rsr(t),o=RSRModule.FS.readFile(a+"out.pdb",{encoding:"utf8"});postMessage({messageId:e.data.messageId,messageTag:"result",result:s,taskName:currentTaskName}),postMessage({messageId:e.data.messageId,messageTag:"pdb_out",jobId:a,result:o,taskName:currentTaskName})}onmessage=function(e){switch(e.data.method){case"loadFile":console.log("Load file(s)",e.data.files),loadFiles(e.data.files);break;case"loadUrl":console.log("Download file(s)",e.data.urls),downLoadFiles(e.data.urls);break;case"get_rama":currentTaskName="get_rama",getRama(e),currentTaskName="";break;case"flip_peptide":console.log("Do peptide-flip in cryst worker ..."),currentTaskName="flip_peptide",flipPeptide(e),currentTaskName="";break;case"mini_rsr":console.log("Do mini-rsr in cryst worker ..."),currentTaskName="mini_rsr",miniRSR(e),currentTaskName="";break;case"auto_fit_rotamer":console.log("Do auto-fit rotamer in cryst worker ..."),currentTaskName="mini_rsr",autoFitRotamer(e),currentTaskName="";break;case"get_bvals":currentTaskName="get_bvals",getBVals(e),currentTaskName="";break;case"density_fit":currentTaskName="density_fit",getDensityFit(e),currentTaskName="";break;case"get_xyz":currentTaskName="get_xyz";const a=/.pdb$/,t=Object.keys(dataObjectsNames.pdbFiles);let s;if(e.data.resInfo.molKey)s=e.data.resInfo.molKey;else for(let o=0;o<t.length;o++){const n=t[o];if(dataObjectsNames.pdbFiles[n].originalFileName.replace(a,"")===e.data.resInfo.molName){console.log("Use",n,dataObjects.pdbFiles[n]),s=n;break}}if(s){let a;a=e.data.resInfo.seqNum?CCP4Module.getXYZSeqNumInsCode(dataObjects.pdbFiles[s].fileName,e.data.resInfo.chain,e.data.resInfo.seqNum,""):CCP4Module.getXYZResNo(dataObjects.pdbFiles[s].fileName,e.data.resInfo.chain,e.data.resInfo.resNo),3===a.size()&&postMessage({messageId:e.data.messageId,messageTag:"result",result:[-a.get(0),-a.get(1),-a.get(2)],taskName:currentTaskName})}currentTaskName="";break;case"get_rotamers":currentTaskName="rotamers",getRotamers(e),currentTaskName="";break;case"draw_cube":currentTaskName="draw_cube",drawCube(e),currentTaskName="";break;case"rotamer_dodecs":currentTaskName="rotamer_dodecs",drawDodos(e),currentTaskName="";break;case"rama_balls":currentTaskName="rama_balls",drawRamaBalls(e),currentTaskName="";break;default:console.log("default, do nothing",e.data.method)}};