import { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';

export default function TorsionPlot({linkage_type, sorted_torsion_list}) {

    const linkage_db = { 
      "ASN-1,2-NAG": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/ASN-NAG_reduced.json', 
      "NAG-1,4-NAG": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/NAG-NAG_reduced.json', 
      "NAG-1,4-MAN": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/NAG-MAN_reduced.json', 
      "MAN-1,3-MAN": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/MAN-MAN_reduced.json', 
      "MAN-1,6-MAN": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/MAN-MAN_reduced.json', 
      "MAN-1,3-BMA": 'https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/BMA-MAN_reduced.json', 

    }

    const bin_db = {
      "ASN-1,2-NAG": {start: 0, end: 360,size: 4},
      "NAG-1,4-NAG": {start: -180, end: 180, size: 4},
      "NAG-1,4-MAN": {start: -180, end: 180, size: 4},
      "MAN-1,3-MAN": {start: -180, end: 180, size: 4},
      "MAN-1,6-MAN": {start: -180, end: 180, size: 4},
      "MAN-1,3-BMA": {start: -180, end: 180, size: 4},
    }

    const [trace, setTrace] = useState({})
    const [overlay, setOverlay] = useState({})

    useEffect(() => {console.log("sorted_torsion_list updated", sorted_torsion_list)},[sorted_torsion_list])

    useEffect(() => {
        fetch(linkage_db[linkage_type])
        .then((response) => response.json())
        .then((responseJson) => {
            let xData = []
            let yData = []

            for (let i = 0; i < responseJson.length; i++) { 
                xData.push(parseFloat(responseJson[i].Phi))
                yData.push(parseFloat(responseJson[i].Psi))
            }

            setTrace({
                x: xData,
                y: yData,
                name: 'density',
                ncontours: 100,
                colorscale: 'Hot',
                reversescale: true,
                showscale: true,
                type: 'histogram2d',
                dragmode: false, 
                
                xbins: {
                    start: -180,
                    end: 180,
                    size: 4
                  },
                  autobiny: false,
                  ybins: bin_db[linkage_type]
              })

        })
        .catch((error) => {
            console.error(error);
            console.log(linkage_type, " is not in the DB most likely ")
        });
        
        let overlay_phi = []
        let overlay_psi = []

        for (let i = 0; i < sorted_torsion_list[linkage_type].length; i++) { 
          overlay_phi.push(sorted_torsion_list[linkage_type][i].phi)
          overlay_psi.push(sorted_torsion_list[linkage_type][i].psi)

        }

        setOverlay({
          x: overlay_phi, 
          y: overlay_psi, 
          mode: 'markers',
          type: 'scatter',
          marker: {
            size: 8,
            color: 'blue',
            symbol: ['x']
        },
        })

    }, [linkage_type])

    return (
      <Plot
        data={[
            trace, overlay 
        ]}
        layout={ {width: 500, height: 500, title: linkage_type, 
        yaxis: {fixedrange: true},
        xaxis : {fixedrange: true},} }
      />
    );
  }

