var CCP4Module;importScripts("web_example.js"),importScripts("papaparse.min.js"),createCCP4Module({print(e){postMessage(["output",e])},printErr(e){postMessage(["output",e])}}).then((function(e){CCP4Module=e})).catch((e=>{console.log("CCP4 problem :("),console.log(e)})),onmessage=function(e){new CCP4Module.VectorString;try{CCP4Module.FS_createDataFile(".",e.data[1],e.data[0],!0,!0)}catch(e){}var t=CCP4Module.calculate_normal_modes(e.data[1],1);const s=t.GetBValues(),o=s.size();let l=[],a=[];const p=CCP4Module.get_CA_bvalues_from_file(e.data[1]);let r;if(o>0){r=0;for(let e=0;e<o;e++)r+=s.get(e)/p.get(e);r/=o,r=1/r}else r=1;const u=CCP4Module.GetCorrelations(t,r);for(let e=0;e<o;e++)l.push({x:e+1,y:s.get(e)*r}),a.push({x:e+1,y:p.get(e)});const n=u.get_rows(),c=u.get_columns();let g={x:[],y:[],z:[]};for(let e=0;e<n;e++)for(let t=0;t<c;t++)g.x.push(e),g.y.push(t),g.z.push(u.get(e,t));postMessage(["corrMat",g]),postMessage(["bvalues",[l,a]]),postMessage(["result",0]);let C=CCP4Module.calculate_normal_modes(e.data[1],0);const d=CCP4Module.GetEigen(C).get(1),m=d.get_rows();d.get_columns(),postMessage(["output","Normal modes"]);let i=[];for(let e=6;e<m;e++)i.push(d.get(e,0));postMessage(["energies",i]);const M=CCP4Module.GetDisplacements(C),_=M.nModes(),f=M.nSteps();let h={modes:[]};for(let e=0;e<_;e++){h.modes.push({steps:[]});for(let t=0;t<f;t++){h.modes[e].steps.push({xyz:[]});const s=M.getDisplacements(e,t),o=s.size();for(let l=0;l<o;l++)h.modes[e].steps[t].xyz.push(s.get(l))}}postMessage(["displacements",h])};