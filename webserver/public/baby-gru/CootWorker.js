var cootModule,molecules_container,ccp4Module;importScripts("./wasm/moorhen.js"),importScripts("./wasm/web_example.js");var guid=function(){var e=Date.now();return"xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g,(function(a){var t=(e+16*Math.random())%16|0;return e=Math.floor(e/16),("x"===a?t:3&t|8).toString(16)}))},print=function(e){console.log(e),postMessage({consoleMessage:JSON.stringify(e)})},instancedMeshToMeshData=function(e,a,t,o){void 0===t&&(t=!1),void 0===o&&(o=1e4);for(var s=[],r=[],n=[],i=[],l=[],c=[],d=[],_=[],m=[],u=e.geom,p=e.markup,g=u.size(),h=0;h<g;h++){var f=t,y=[],v=[],S=[],M=[],b=[],T=[],z=[],I=u.get(h);"spherical-atoms"===I.name&&(f=!0);for(var x=I.vertices,A=I.triangles,D=A.size(),F=0;F<D;F++){var k=A.get(F).point_id;a?(y.push(k[0]),y.push(k[2]),y.push(k[1])):(y.push(k[0]),y.push(k[1]),y.push(k[2]))}A.delete();for(var w=x.size(),N=0;N<w;N++){var J=x.get(N),C=J.pos;v.push(C[0]),v.push(C[1]),v.push(C[2]);var E=J.normal;S.push(E[0]),S.push(E[1]),S.push(E[2]),J.delete()}x.delete();var R=I.instancing_data_A,P=R.size();if(P>0)for(var U=0;U<P;U++){var O=(H=R.get(U)).position;T.push(O[0]),T.push(O[1]),T.push(O[2]);var B=H.colour;b.push(B[0]),b.push(B[1]),b.push(B[2]),b.push(B[3]);var j=H.size;M.push(j[0]),M.push(j[1]),M.push(j[2]),z.push.apply(z,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]),H.delete()}R.delete();var L=I.instancing_data_B,q=L.size();if(q>0)for(U=0;U<q;U++){var H;if(!((j=(H=L.get(U)).size)[2]>o)){M.push(j[0]),M.push(j[1]),M.push(j[2]),O=H.position,T.push(O[0]),T.push(O[1]),T.push(O[2]),B=H.colour,b.push(B[0]),b.push(B[1]),b.push(B[2]),b.push(B[3]);var G=H.orientation;z.push(G[0][0]),z.push(G[0][1]),z.push(G[0][2]),z.push(G[0][3]),z.push(G[1][0]),z.push(G[1][1]),z.push(G[1][2]),z.push(G[1][3]),z.push(G[2][0]),z.push(G[2][1]),z.push(G[2][2]),z.push(G[2][3]),z.push(G[3][0]),z.push(G[3][1]),z.push(G[3][2]),z.push(G[3][3]),H.delete()}}L.delete(),I.delete(),n.push(S),r.push(v),s.push(y),i.push(M),c.push(T),d.push(z),l.push(b),_.push(!0),f?m.push("PERFECT_SPHERES"):m.push("TRIANGLES")}u.delete();var K=simpleMeshToMeshData(p);return e.delete(),K.idx_tri.length>0&&K.idx_tri[0].length>0&&K.idx_tri[0][0].length>0?t?{prim_types:[m,K.prim_types[0]],idx_tri:[s,K.idx_tri[0]],vert_tri:[c,K.vert_tri[0]],norm_tri:[n,K.norm_tri[0]],col_tri:[l,K.col_tri[0]],instance_use_colors:[_,null],instance_sizes:[i,null],instance_origins:[c,null],instance_orientations:[d,null]}:{prim_types:[m,K.prim_types[0]],idx_tri:[s,K.idx_tri[0]],vert_tri:[r,K.vert_tri[0]],norm_tri:[n,K.norm_tri[0]],col_tri:[l,K.col_tri[0]],instance_use_colors:[_,null],instance_sizes:[i,null],instance_origins:[c,null],instance_orientations:[d,null]}:{prim_types:[m],idx_tri:[s],vert_tri:[r],norm_tri:[n],col_tri:[l],instance_use_colors:[_],instance_sizes:[i],instance_origins:[c],instance_orientations:[d]}},simpleMeshToMeshData=function(e,a){void 0===a&&(a=!1);for(var t=e.vertices,o=e.triangles,s=[],r=[],n=[],i=[],l=o.size(),c=0;c<l;c++){var d=o.get(c).point_id;a?s.push.apply(s,[d[0],d[2],d[1]]):s.push.apply(s,d)}o.delete();var _=t.size();for(c=0;c<_;c++){var m=t.get(c),u=m.pos,p=m.normal,g=m.color;r.push.apply(r,u),a?n.push.apply(n,[-p[0],-p[1],-p[2]]):n.push.apply(n,p),i.push.apply(i,g),m.delete()}return t.delete(),e.delete(),{prim_types:[["TRIANGLES"]],idx_tri:[[s]],vert_tri:[[r]],norm_tri:[[n]],col_tri:[[i]]}},SuperposeResultsToJSArray=function(e){for(var a=e.aligned_pairs,t=a.size(),o=[],s=0;s<t;s++){var r=a.get(s),n=r.first,i=n.residue_spec,l=r.second,c=l.residue_spec,d={reference:{chainId:i.chain_id,insCode:i.ins_code,seqNum:i.res_no,restype:"UNK",value:n.function_value,label:n.label},moving:{chainId:c.chain_id,insCode:c.ins_code,seqNum:c.res_no,restype:"UNK",value:l.function_value,label:l.label}};l.delete(),c.delete(),n.delete(),i.delete(),o.push(d)}return a.delete(),{referenceSequence:e.alignment.first,movingSequence:e.alignment.second,supperposeInfo:e.superpose_info,alignedPairsData:o}},colourRulesToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push(s)}return e.delete(),a},floatArrayToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push(s)}return e.delete(),a},mapMoleculeCentreInfoToJSObject=function(e){var a=e.updated_centre,t={updated_centre:[a.x(),a.y(),a.z()],success:e.success,suggested_contour_level:e.suggested_contour_level};return a.delete(),t},intArrayToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push(s)}return e.delete(),a},stringArrayToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push(s)}return e.delete(),a},symmetryToJSData=function(e){for(var a=[],t=e.first,o=e.second,s=t.cell,r=t.symm_trans,n=r.size(),i=0;i<n;i++){var l=r.get(i),c=l.first,d=l.second,_=o.get(i);a.push({x:c.x(),y:c.y(),z:c.z(),asString:c.symm_as_string,isym:c.isym(),us:d.us,ws:d.ws,vs:d.vs,matrix:_}),c.delete()}return s.delete(),r.delete(),o.delete(),t.delete(),a},mmrrccStatsToJSArray=function(e){var a=function(e){for(var a=[],t=e.keys(),o=t.size(),s=0;s<o;s++){var r=t.get(s),n=e.get(r);a.push({resNum:r.res_no,insCode:r.ins_code,modelNumber:r.model_number,chainId:r.chain_id,n:n.n,correlation:n.correlation()}),r.delete(),n.delete()}return t.delete(),a},t=e.first,o=e.second,s={"All atoms":a(t),"Side-chains":a(o)};return t.delete(),o.delete(),s},residueSpecToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push({resNum:s.res_no,insCode:s.ins_code,modelNumber:s.model_number,chainId:s.chain_id}),s.delete()}return e.delete(),a},validationDataToJSArray=function(e,a){void 0===a&&(a=null);for(var t=[],o=e.cviv,s=o.size(),r=0;r<s;r++){var n=o.get(r);if(null!==a&&n.chain_id!==a);else{for(var i=n.rviv,l=i.size(),c=0;c<l;c++){var d=i.get(c),_=d.residue_spec;t.push({chainId:_.chain_id,insCode:_.ins_code,seqNum:_.res_no,restype:"UNK",value:d.function_value}),d.delete(),_.delete()}i.delete()}n.delete()}return o.delete(),e.delete(),t},linesBoxToJSArray=function(e){for(var a=[],t=e.line_segments,o=t.size(),s=0;s<o;s++){for(var r=[],n=t.get(s),i=n.size(),l=0;l<i;l++){var c=n.get(l),d=c.getStart(),_=c.getFinish(),m=c.amplitude(),u={x:d.x(),y:d.y(),z:d.z()},p={x:_.x(),y:_.y(),z:_.z()};r.push({start:u,end:p,dist:m}),d.delete(),_.delete(),c.delete()}n.delete(),a.push(r)}return t.delete(),e.delete(),a},vectorHBondToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o);a.push({hb_hydrogen:s.hb_hydrogen,donor:s.donor,acceptor:s.acceptor,donor_neigh:s.donor_neigh,acceptor_neigh:s.acceptor_neigh,angle_1:s.angle_1,angle_2:s.angle_2,angle_3:s.angle_3,dist:s.dist,ligand_atom_is_donor:s.ligand_atom_is_donor,hydrogen_is_ligand_atom:s.hydrogen_is_ligand_atom,bond_has_hydrogen_flag:s.bond_has_hydrogen_flag})}return e.delete(),a},interestingPlaceDataToJSArray=function(e){for(var a=[],t=e.size(),o=0;o<t;o++){var s=e.get(o),r=s.residue_spec;a.push({modelNumber:r.model_number,chainId:r.chain_id,insCode:r.ins_code,resNum:r.res_no,featureType:s.feature_type,featureValue:s.feature_value,buttonLabel:s.button_label,badness:s.badness,coordX:s.x,coordY:s.y,coordZ:s.z}),s.delete(),r.delete()}return e.delete(),a},ramachandranDataToJSArray=function(e,a){for(var t=[],o=e.size(),s=0;s<o;s++){var r=e.get(s),n=r.phi_psi;n.chain_id===a&&t.push({chainId:n.chain_id,insCode:n.ins_code,seqNum:n.residue_number,restype:r.residue_name(),isOutlier:!r.is_allowed_flag,phi:n.phi(),psi:n.psi(),is_pre_pro:"PRO"===r.residue_name()}),r.delete(),n.delete()}return e.delete(),t},simpleMeshToLineMeshData=function(e,a){for(var t=e.vertices,o=e.triangles,s=[],r=[],n=[],i=[],l=o.size(),c=0;c<l;c++){var d=o.get(c).point_id;s.push.apply(s,[d[0],d[1],d[0],d[2],d[1],d[2]])}o.delete();var _=t.size();for(c=0;c<_;c++){var m=t.get(c);r.push.apply(r,m.pos),n.push.apply(n,m.normal),i.push.apply(i,m.color),m.delete()}return t.delete(),e.delete(),a?{prim_types:[["NORMALLINES"]],useIndices:[[!0]],idx_tri:[[s]],vert_tri:[[r]],additional_norm_tri:[[n]],norm_tri:[[n]],col_tri:[[i]]}:{prim_types:[["LINES"]],useIndices:[[!0]],idx_tri:[[s]],vert_tri:[[r]],norm_tri:[[n]],col_tri:[[i]]}},read_pdb=function(e,a){var t=guid();cootModule.FS_createDataFile(".","".concat(t,".pdb"),e,!0,!0);var o="./".concat(t,".pdb"),s=molecules_container.read_pdb(o);return cootModule.FS_unlink(o),s},auto_open_mtz=function(e){var a=guid(),t=new Uint8Array(e);cootModule.FS_createDataFile(".","".concat(a,".mtz"),t,!0,!0);var o="./".concat(a,".mtz"),s=molecules_container.auto_read_mtz(o);return cootModule.FS_unlink(o),s},read_dictionary=function(e,a){var t=guid();cootModule.FS_createDataFile(".","".concat(t,".cif"),e,!0,!0);var o="./".concat(t,".cif"),s=molecules_container.import_cif_dictionary(o,a);return cootModule.FS_unlink(o),s},replace_molecule_by_model_from_file=function(e,a){var t=guid(),o="./".concat(t,".pdb");cootModule.FS_createDataFile(".",o,a,!0,!0);var s=molecules_container.replace_molecule_by_model_from_file(e,o);return cootModule.FS_unlink(o),s},replace_map_by_mtz_from_file=function(e,a,t){var o=guid(),s="./".concat(o,".mtz"),r=new Uint8Array(a);cootModule.FS_createDataFile(".",s,r,!0,!0);var n=[e,s,t.F,t.PHI,"",!1],i=molecules_container.replace_map_by_mtz_from_file.apply(molecules_container,n);return cootModule.FS_unlink(s),i},new_positions_for_residue_atoms=function(e,a){var t=0,o=new cootModule.Vectormoved_residue_t;return a.forEach((function(e){if(e.length>0){var a=e[0].label.split("/"),t=a[3].split("."),s=t[0],r=t[1];r=r||"";var n=new cootModule.moved_residue_t(a[2],parseInt(s),r);e.forEach((function(e){var a=new cootModule.moved_atom_t(e.name,e.alt_loc,e.x,e.y,e.z,-1);n.add_atom(a),a.delete()})),o.push_back(n),n.delete()}})),t+=molecules_container.new_positions_for_atoms_in_residues(e,o),o.delete(),t},read_mtz=function(e,a,t){var o=guid(),s=new Uint8Array(e);cootModule.FS_createDataFile(".","".concat(o,".mtz"),s,!0,!0);var r="./".concat(o,".mtz"),n=[r,t.F,t.PHI,"",!1,t.isDifference],i=molecules_container.read_mtz.apply(molecules_container,n);return cootModule.FS_unlink(r),i},associate_data_mtz_file_with_map=function(e,a,t,o,s){var r=new Uint8Array(a.data);cootModule.FS_createDataFile(".","".concat(a.fileName,".mtz"),r,!0,!0);var n="./".concat(a.fileName,".mtz"),i=[e,n,t,o,s];return molecules_container.associate_data_mtz_file_with_map.apply(molecules_container,i),n},read_ccp4_map=function(e,a,t){var o=guid(),s=new Uint8Array(e);cootModule.FS_createDataFile(".","".concat(o,".map"),s,!0,!0);var r="./".concat(o,".map"),n=[r,t],i=molecules_container.read_ccp4_map.apply(molecules_container,n);return cootModule.FS_unlink(r),i},doColourTest=function(e){console.log("DEBUG: Start test...");var a={0:{cid:"//A/1-10/",rgb:[255,0,0]},1:{cid:"//A/11-20/",rgb:[0,255,0]},2:{cid:"//A/21-30/",rgb:[0,0,255]}},t=new cootModule.MapIntFloat3;for(var o in Object.keys(a))t[o]=a[o].rgb;var s=new cootModule.VectorStringUInt_pair;for(var o in Object.keys(a)){var r={first:a[o].cid,second:parseInt(o)};s.push_back(r)}console.log("DEBUG: Running molecules_container.set_user_defined_bond_colours"),molecules_container.set_user_defined_bond_colours(e,t),console.log("DEBUG: Running molecules_container.set_user_defined_atom_colour_by_residue"),molecules_container.set_user_defined_atom_colour_by_residue(e,s),s.delete(),t.delete()};onmessage=function(e){if("CootInitialize"===e.data.message)createRSRModule({locateFile:function(e){return"./wasm/".concat(e)},onRuntimeInitialized:function(){},mainScriptUrlOrBlob:"moorhen.js",print,printErr:print}).then((function(a){postMessage({consoleMessage:"Initialized molecules_container",message:e.data.message,messageId:e.data.messageId}),(molecules_container=new(cootModule=a).molecules_container_js(!1)).set_show_timings(!1),molecules_container.fill_rotamer_probability_tables(),molecules_container.set_map_sampling_rate(1.7),cootModule.FS.mkdir("COOT_BACKUP")})).catch((function(e){console.log(e),print(e)})),createCCP4Module({locateFile:function(e){return"./wasm/".concat(e)},onRuntimeInitialized:function(){},mainScriptUrlOrBlob:"web_example.js",print,printErr:print}).then((function(e){ccp4Module=e})).catch((function(e){console.log(e),print(e)}));else if("get_atoms"===e.data.message){var a=guid(),t="./".concat(a,".pdb");"pdb"===e.data.format?molecules_container.writePDBASCII(e.data.molNo,t):"mmcif"===e.data.format?molecules_container.writeCIFASCII(e.data.molNo,t):console.log("Unrecognised format... ".concat(e.data.format));var o=cootModule.FS.readFile(t,{encoding:"utf8"});cootModule.FS_unlink(t),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,consoleMessage:"Fetched coordinates of molecule ".concat(e.data.molNo),message:e.data.message,result:{molNo:e.data.molNo,pdbData:o}})}else if("get_mtz_data"===e.data.message){var s=cootModule.FS.readFile(e.data.fileName,{encoding:"binary"});postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,consoleMessage:"Fetched mtz data for map ".concat(e.data.molNo),message:e.data.message,result:{molNo:e.data.molNo,mtzData:s}})}else if("get_map"===e.data.message){a=guid(),t="./".concat(a,".map"),molecules_container.writeCCP4Map(e.data.molNo,t);var r=cootModule.FS.readFile(t,{encoding:"binary"});cootModule.FS_unlink(t),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,consoleMessage:"Fetched map of map ".concat(e.data.molNo),message:e.data.message,result:{molNo:e.data.molNo,mapData:r.buffer}})}else if("read_mtz"===e.data.message)try{a=guid(),cootModule.FS_createDataFile(".","".concat(a,".mtz"),e.data.data,!0,!0,!0),t="./".concat(a,".mtz");var n=molecules_container.read_mtz(t,"FWT","PHWT","",!1,!1);cootModule.FS_unlink(t),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,consoleMessage:"Read map MTZ as molecule ".concat(n),message:e.data.message,result:{molNo:n,name:e.data.name}})}catch(e){print(e)}else if("get_rama"===e.data.message){a=guid(),t="./".concat(a,".pdb"),molecules_container.writePDBASCII(e.data.molNo,t);var i=cootModule.getRamachandranData(t,e.data.chainId);cootModule.FS_unlink(t);for(var l=[],c=0;c<i.size();c++){var d=i.get(c),_={chainId:d.chainId,insCode:d.insCode,seqNum:d.seqNum,restype:d.restype,phi:d.phi,psi:d.psi,isOutlier:d.isOutlier,is_pre_pro:d.is_pre_pro};l.push(_)}postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,messageTag:"result",result:l})}else if("copy_fragment"===e.data.message){var m=molecules_container.copy_fragment_using_residue_range(e.data.molNo,e.data.chainId,e.data.res_no_start,e.data.res_no_end);postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,messageTag:"result",result:m})}else"delete"===e.data.message?(i=molecules_container.close_molecule(e.data.molNo),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,messageTag:"result",result:i})):"delete_file_name"===e.data.message&&(i=cootModule.FS_unlink(e.data.fileName),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,messageTag:"result",result:i}));if("coot_command"===e.data.message){var u=e.data,p=u.returnType,g=u.command,h=u.commandArgs,f=u.messageId;try{var y,v="Message from main thread to worker took ".concat(Date.now()-e.data.myTimeStamp," ms (").concat(g,") - (").concat(f.slice(0,5),")"),S=new Date;y="shim_read_pdb"===g?read_pdb.apply(void 0,h):"shim_new_positions_for_residue_atoms"===g?new_positions_for_residue_atoms.apply(void 0,h):"shim_read_mtz"===g?read_mtz.apply(void 0,h):"shim_auto_open_mtz"===g?auto_open_mtz.apply(void 0,h):"shim_read_ccp4_map"===g?read_ccp4_map.apply(void 0,h):"shim_read_dictionary"===g?read_dictionary.apply(void 0,h):"shim_associate_data_mtz_file_with_map"===g?associate_data_mtz_file_with_map.apply(void 0,h):"shim_replace_molecule_by_model_from_file"===g?replace_molecule_by_model_from_file.apply(void 0,h):"shim_replace_map_by_mtz_from_file"===g?replace_map_by_mtz_from_file.apply(void 0,h):"shim_do_colour_test"===g?doColourTest.apply(void 0,h):"shim_smiles_to_pdb"===g?cootModule.SmilesToPDB.apply(cootModule,h):molecules_container[g].apply(molecules_container,h);var M=new Date,b=M-S,T="libcootAPI command ".concat(g," took ").concat(b," ms  - (").concat(f.slice(0,5),")"),z=void 0;switch(S=new Date,p){case"instanced_mesh_perm":z=instancedMeshToMeshData(y,!0);break;case"symmetry":z=symmetryToJSData(y);break;case"mmrrcc_stats":z=mmrrccStatsToJSArray(y);break;case"colour_rules":z=colourRulesToJSArray(y);break;case"instanced_mesh_perfect_spheres":z=instancedMeshToMeshData(y,!1,!0);break;case"instanced_mesh":z=instancedMeshToMeshData(y,!1,!1,5);break;case"mesh_perm":z=simpleMeshToMeshData(y,!0);break;case"mesh":z=simpleMeshToMeshData(y);break;case"lit_lines_mesh":z=simpleMeshToLineMeshData(y,!0);break;case"lines_mesh":z=simpleMeshToLineMeshData(y,!1);break;case"float_array":z=floatArrayToJSArray(y);break;case"int_array":z=intArrayToJSArray(y);break;case"map_molecule_centre_info_t":z=mapMoleculeCentreInfoToJSObject(y);break;case"string_array":z=stringArrayToJSArray(y);break;case"residue_specs":z=residueSpecToJSArray(y);break;case"ramachandran_data":z=ramachandranDataToJSArray(y,e.data.chainID);break;case"validation_data":z=validationDataToJSArray(y,e.data.chainID);break;case"interesting_places_data":z=interestingPlaceDataToJSArray(y);break;case"superpose_results":z=SuperposeResultsToJSArray(y);break;case"generic_3d_lines_bonds_box":z=linesBoxToJSArray(y);break;case"vector_hbond":z=vectorHBondToJSArray(y);break;case"status_instanced_mesh_pair":z={status:y.first,mesh:instancedMeshToMeshData(y.second,!1,!1,5)};break;default:z=y}b=(M=new Date)-S;var I="conversion of output of ".concat(g," to JS data took ").concat(b," ms  - (").concat(f.slice(0,5),")");postMessage({timelibcootAPI:T,timeconvertingWASMJS:I,timeMainThreadToWorker:v,messageId:f,messageSendTime:Date.now(),consoleMessage:"Completed ".concat(g," in ").concat(Date.now()-e.data.myTimeStamp," ms"),result:{status:"Completed",result:z}})}catch(a){console.log(a),postMessage({messageId:e.data.messageId,myTimeStamp:e.data.myTimeStamp,message:e.data.message,consoleMessage:"EXCEPTION RAISED IN ".concat(g,", ").concat(a),result:{status:"Exception"}})}}};