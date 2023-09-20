import { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { linkage_db, bin_db} from '../../data/Constants';

export default function TorsionPlot({linkage_type, sorted_torsion_list}) {

 
    const [trace, setTrace] = useState({})
    const [overlay, setOverlay] = useState({})
    const [linkageFound, setLinkageFound] = useState({})

    useEffect(() => {
        fetch(linkage_db[linkage_type])
        .then((response) => response.json())
        .then((responseJson) => {
            let xData = []
            let yData = []

            for (let i = 0; i < responseJson.length; i++) { 
                xData.push(parseFloat(responseJson[i].phi))
                yData.push(parseFloat(responseJson[i].psi))
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
              setLinkageFound(true)
        })
        .catch((error) => {
            console.error(error);
            console.log(linkage_type, " is not in the DB most likely ")
            setLinkageFound(false)
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
      !linkageFound? (
        <h3>{linkage_type} does not have enough datapoints to generate a torsion plot.</h3>
      ): 
      <Plot
        data={[
            trace, overlay 
        ]}
        layout={ {width: 500, height: 500, title: linkage_type, 
        yaxis: {
          fixedrange: true, 
          range: (linkage_type in bin_db)? [bin_db[linkage_type].start, bin_db[linkage_type].end] : [-180,180],
          showgrid:false
        },
        xaxis : {
          fixedrange: true, 
          range:[-180,180], 
          showgrid:false
        },} }
      />
    );
  }

