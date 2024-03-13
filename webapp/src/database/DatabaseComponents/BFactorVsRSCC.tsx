import React, { useEffect, useState, lazy } from 'react';

const Plot = lazy(async () => await import('react-plotly.js'));

function calculatePoints(data): [number[], number[], string[], number[], number[], string[]] {
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
                    errorXAxis.push(sugars[j].bFactor as number);
                    errorYAxis.push(sugars[j].rscc as number);
                    errorText.push(sugars[j].sugarId as string);
                } else {
                    xAxis.push(sugars[j].bFactor as number);
                    yAxis.push(sugars[j].rscc as number);
                    text.push(sugars[j].sugarId as string);
                }
            }
        }
    }

    return [xAxis, yAxis, text, errorXAxis, errorYAxis, errorText];
}

export default function BFactorVsRSCC(props) {
    const [trace, setTrace] = useState({});
    const [corrTrace, setCorrTrace] = useState({});
    const [badTrace, setBadTrace] = useState({});

    useEffect(() => {
        const [xAxis, yAxis, text, errorXAxis, errorYAxis, errorText] = calculatePoints(props);

        const maxX = Math.max(Math.max(...xAxis), Math.max(...errorXAxis)) + 5;
        console.log(errorXAxis, xAxis)
        setCorrTrace({
            x: [0, maxX],
            y: [0.7, 0.7],
            text,
            hoverinfo: 'text',
            fill: 'tozeroy',
            fillcolor: 'rgba(173,181,189,0.3)',
            fillopacity: 0.1,
            marker: {
                size: 1,
                color: 'rgba(173,181,189,0.5)',
                symbol: 'x',
            },
            showlegend: false,
        });

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
                symbol: 'o',
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
                symbol: 'x',
            },
            name: 'Issues',
        });
    }, [props]);

    return (
        <div className="flex flex-col mx-auto">
            <span className="text-xl">BFactor vs RSCC</span>
            <Plot
                data={[trace, badTrace, corrTrace]}
                layout={{
                    showlegend: true,
                    width: 500,
                    height: 400,
                    title: '',
                    legend: {
                        x: 0,
                        xanchor: 'left',
                        y: 0,
                        bgcolor: 'rgba(0,0,0,0)',
                        borderwidth: 0.2,
                        bordercolor: 'gray',
                    },
                    plot_bgcolor: '#FFFFF',
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
                            text: 'RSCC',
                        },
                        fixedrange: true,
                        range: [0, 1],
                        showgrid: true,
                    },
                    xaxis: {
                        title: {
                            text: 'B Factor',
                        },
                        fixedrange: true,
                        // range: [0, 100],
                        showgrid: true,
                    },
                }}
            />
        </div>
    );
}
