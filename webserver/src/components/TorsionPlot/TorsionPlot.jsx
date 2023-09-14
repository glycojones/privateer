import { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';

export default function TorsionPlot({}) {

    const [trace, setTrace] = useState({})

    useEffect(() => {
        fetch('https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/ASN-NAG_reduced.json')
        .then((response) => response.json())
        .then((responseJson) => {
            console.log(responseJson)

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
                  ybins: {
                    start: 0,
                    end: 360,
                    size: 4
                  },
              })

        })
        .catch((error) => {
            console.error(error);
        });
    }, [])

    return (
      <Plot
        data={[
            trace
        ]}
        layout={ {width: 500, height: 500, title: 'ASN-NAG ', 
        yaxis: {fixedrange: true},
        xaxis : {fixedrange: true},} }
      />
    );
  }

