import React, { useEffect, useState, lazy } from 'react';

const Plot = lazy(async () => await import('react-plotly.js'));

function calculatePoints(data) {
    const glycans = data.data.glycans;

    const xAxis: number[] = [];
    const yAxis: number[] = [];
    const text: string[] = [];

    const errorXAxis: number[] = [];
    const errorYAxis: number[] = [];
    const errorText: string[] = [];

    for (const key in glycans) {
        const glycanType = glycans[key];
        for (let i = 0; i < glycanType.length; i++) {
            const sugars = glycanType[i].sugars;
            for (let j = 0; j < sugars.length; j++) {
                if (sugars[j].diagnostic !== 'yes') {
                    errorXAxis.push(sugars[j].phi as number);
                    errorYAxis.push(sugars[j].theta as number);
                    errorText.push(sugars[j].sugarId as string);
                } else {
                    xAxis.push(sugars[j].phi as number);
                    yAxis.push(sugars[j].theta as number);
                    text.push(sugars[j].sugarId as string);
                }
                // console.log(sugars[j])
            }
        }
    }

    return [xAxis, yAxis, text, errorXAxis, errorYAxis, errorText];
}

export default function CremerPopleGraph(props: any) {
    const [trace, setTrace] = useState({});
    const [badTrace, setBadTrace] = useState({});

    useEffect(() => {
        const [xAxis, yAxis, text, errorXAxis, errorYAxis, errorText] =
            calculatePoints(props);

        setTrace({
            x: xAxis,
            y: yAxis,
            text,
            hoverinfo: 'text',
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 8,
                color: 'blue',
                symbol: ['o'],
            },
            name: 'No Issues',
        });
        setBadTrace({
            x: errorXAxis,
            y: errorYAxis,
            text: errorText,
            hoverinfo: 'text',
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 8,
                color: 'red',
                symbol: ['o'],
            },
            name: 'Issues',
        });
    }, [props]);

    return (
        <div className="flex flex-col mx-auto">
            <span className="text-xl">
                Conformational landscape for pyranoses
            </span>

            <Plot
                data={[trace, badTrace]}
                layout={{
                    showlegend: true,
                    legend: {
                        x: 1,
                        xanchor: 'right',
                        y: 0.5,
                        bgcolor: 'rgba(0,0,0,0)',
                        borderwidth: 0.2,
                        bordercolor: "gray"


                    },
                    width: 500,
                    height: 400,
                    title: '',
                    plot_bgcolor: '#FFFFFF',
                    paper_bgcolor: '#D6D9E5',
                    margin: {
                        l: 50,
                        r: 50,
                        b: 50,
                        t: 10,
                        pad: 4,
                    },
                    yaxis: {
                        title: {
                            text: 'Theta / °',
                        },
                        fixedrange: true,
                        range: [180, 0],
                        showgrid: true,
                    },
                    xaxis: {
                        title: {
                            text: 'Phi / °',
                        },
                        fixedrange: true,
                        range: [360, 0],
                        showgrid: true,
                    },
                }}
            />
        </div>
    );
}
