import { useEffect, useState, lazy } from 'react';

const Plot = lazy(() => import('react-plotly.js'));

function calculate_points(data) {
    let glycans = data.data.glycans

    let x_axis = []
    let y_axis = []
    let text = []

    let error_x_axis = []
    let error_y_axis = []
    let error_text = []

    for (const key in glycans) {
        let glycan_type = glycans[key]
        for (let i = 0; i < glycan_type.length; i++) {
            let sugars = glycan_type[i].Sugars
            for (let j = 0; j < sugars.length; j++) {
                if (sugars[j].Diagnostic != "yes") { 
                    error_x_axis.push(sugars[j].Phi)
                    error_y_axis.push(sugars[j].Theta)
                    error_text.push(sugars[j]["Sugar ID"])
                }
                else { 
                    x_axis.push(sugars[j].Phi)
                    y_axis.push(sugars[j].Theta)
                    text.push(sugars[j]["Sugar ID"])
                }
                console.log(sugars[j])
               
            }
        }
    }


    return [x_axis, y_axis, text, error_x_axis, error_y_axis, error_text]
}

export default function CremerPopleGraph(props) {

    const [trace, setTrace] = useState({})
    const [badTrace, setBadTrace] = useState({})

    useEffect(() => {
        const [x_axis, y_axis, text, error_x_axis, error_y_axis, error_text] = calculate_points(props)

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
            name: "No Issues",
        })
        setBadTrace({
            x: error_x_axis,
            y: error_y_axis,
            text: error_text,
            hoverinfo: "text",
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 8,
                color: 'red',
                symbol: ['o']
            },
            name: "Issues",

        })

    }, [])


    return (
        <div className='flex flex-col mx-auto'>
            <span className='text-xl'>Conformational landscape for pyranoses</span>

            <Plot
                data={[
                    trace, badTrace
                ]}
                layout={{
                    showlegend: false,
                    width: 500, height: 400, title: "",
                    plot_bgcolor: "#FFFFFF",
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
                            text: "Theta / °"
                        },
                        fixedrange: true,
                        range: [180, 0],
                        showgrid: true
                    },
                    xaxis: {
                        title: {
                            text: "Phi / °"
                        },
                        fixedrange: true,
                        range: [360, 0],
                        showgrid: true
                    },
                }}
            />

        </div>

    )



}