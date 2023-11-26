import { useEffect, useState, lazy } from 'react';

const Plot = lazy(() => import('react-plotly.js'));

function calculate_points(data) {
    let glycans = data.data.glycans

    let x_axis = []
    let y_axis = []
    let text = []

    for (const key in glycans) {
        let glycan_type = glycans[key]
        for (let i = 0; i < glycan_type.length; i++) {
            let sugars = glycan_type[i].Sugars
            for (let j = 0; j < sugars.length; j++) {
                x_axis.push(sugars[j].BFactor)
                y_axis.push(sugars[j].RSCC)
                text.push(sugars[j]["Sugar ID"])
            }
        }
    }


    return [x_axis, y_axis, text]
}

export default function BFactorVsRSCC(props) {

    const [trace, setTrace] = useState({})

    useEffect(() => {
        const [x_axis, y_axis, text] = calculate_points(props)

        setTrace({
            x: x_axis,
            y: y_axis,
            text: text,
            hoverinfo: "text",
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 8,
                color: 'green',
                symbol: ['o']
            },
        })
    }, [])


    return (
        <div className='flex flex-col mx-auto'>
            <span className='text-xl'>BFactor vs RSCC</span>
            <Plot
                data={[
                    trace
                ]}
                layout={{
                    width: 500, height: 400, title: "",
                    plot_bgcolor: "#D6D9E5",
                    paper_bgcolor: "#D6D9E5",
                    margin: {
                        l: 50,
                        r: 50,
                        b: 50,
                        t: 10,
                        pad: 4
                    },
                    yaxis: {
                        title: {
                            text: "RSCC"
                        },
                        fixedrange: true,
                        range: [0, 1],
                        showgrid: true
                    },
                    xaxis: {
                        title: {
                            text: "B Factor"
                        },
                        fixedrange: true,
                        range: [0, 100],
                        showgrid: true
                    },
                }}
            />
        </div>
    )



}