var CCP4Module;importScripts("web_example.js"),importScripts("papaparse.min.js"),createCCP4Module({print(t){postMessage(["output",t])},printErr(t){postMessage(["output",t])}}).then((function(t){CCP4Module=t})).catch((t=>{console.log("CCP4 problem :("),console.log(t)})),onmessage=function(t){var e=new CCP4Module.VectorString;e.push_back("gesamt");var a=new CCP4Module.VectorString;for(let s=0;s<t.data[0].length;s++){try{CCP4Module.FS_createDataFile(".",t.data[1][s],t.data[0][s],!0,!0)}catch(t){}e.push_back(t.data[1][s]),a.push_back("dummy")}e.push_back("-csv"),e.push_back("out.csv");var s=CCP4Module.gesamt(e),l=CCP4Module.FS.readFile("out.csv",{encoding:"utf8"});let r=Papa.parse(l).data;if(2==t.data[0].length){let e=!1,a=!1,s=[],l=[],p=-1;for(let t=0;t<r.length;t++){if(e&&t<p+4&&(l.push(parseFloat(r[t][0])),l.push(parseFloat(r[t][1])),l.push(parseFloat(r[t][2])),l.push(parseFloat(r[t][3]))),a){if(r[t][0].length>0){let e=parseFloat(r[t][0]),a=parseFloat(r[t][4].trim().split(" ")[r[t][4].trim().split(" ").length-1]);s.push({x:a,y:e})}if(r[t].length<5)break}4==r[t].length&&"Rx"===r[t][0].trim()&&"Ry"===r[t][1].trim()&&"Rz"===r[t][2].trim()&&"T"===r[t][3].trim()&&(e=!0,p=t),r[t].length>4&&"Dist [A]"===r[t][0].trim()&&"Query"===r[t][3].trim()&&"Target"===r[t][4].trim()&&(a=!0)}12==l.length&&(l.push(0),l.push(0),l.push(0),l.push(1));let o={};o.alignData=s,o.transformMatrices=[l],t.data.length>2&&(o.jobid=t.data[2]),postMessage(["csvResult",o])}else{let e=!1,a=!1,s=[],l=[],p=[],o=-1;for(let i=0;i<r.length;i++){if(e&&i<o+4&&(l.push(parseFloat(r[i][0])),l.push(parseFloat(r[i][1])),l.push(parseFloat(r[i][2])),l.push(parseFloat(r[i][3]))),12==l.length&&(l.push(0),l.push(0),l.push(0),l.push(1),p.push(l),l=[],e=!1),4==r[i].length&&"Rx"===r[i][0].trim()&&"Ry"===r[i][1].trim()&&"Rz"===r[i][2].trim()&&"T"===r[i][3].trim()&&(e=!0,o=i),a){if(r[i][0].length>0){let t=parseFloat(r[i][0]),e=parseFloat(r[i][1].trim().split(" ")[r[i][1].trim().split(" ").length-1]);s.push({x:e,y:t})}if(r[i].length<t.data[0].length+1)break}"Disp.[A]"===r[i][0].trim()&&i>2&&"RESIDUE ALIGNMENT"===r[i-2][0].trim()&&(a=!0)}let i={};i.transformMatrices=p,i.alignData=s,t.data.length>2&&(i.jobid=t.data[2]),postMessage(["csvResult",i])}postMessage(["result",s])};