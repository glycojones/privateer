import React, { useEffect, useState } from 'react';
import Loading from '../../shared/Loading/Loading.tsx';
import Plot from 'react-plotly.js';

export default function ConformationErrorsVsYear(props: { database: string }) {
    const [totalTrace, setTotalTrace] = useState();
    const [errorTrace, setErrorTrace] = useState();
    const [relativeTrace, setRelativeTrace] = useState();

    // const [depositedTrace, setDepositedTrace] = useState();
    const d = new Date();

    const [data, setData] = useState<Record<
        string,
        Record<string, Record<string, number>>
    > | null>(null);

    const [width, setWidth] = useState(1000);
    const [height, setHeight] = useState(1000);
    const [legendDown, setLegendDown] = useState(false);

    useEffect(() => {
        const url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/stats/conformational_errors_per_year.json';

        fetch(url)
            .then(async (response) => await response.json())
            .then((data) => {
                setData(
                    data as Record<
                        string,
                        Record<string, Record<string, number>>
                    >
                );
            })
            .catch((error) => {
                console.error('Error:', error);
            });

        function handleResize() {
            const currentWidth = window.innerWidth;

            if (currentWidth < 1000) {
                setWidth(0.75 * currentWidth);
                setHeight(600);
                setLegendDown(true);
            } else {
                setWidth(Math.max(0.75 * currentWidth, 1000));
                setHeight(0.5 * Math.max(0.75 * currentWidth, 1000));

                setLegendDown(false);
            }
        }

        window.addEventListener('resize', handleResize);
        handleResize();
    }, []);

    useEffect(() => {
        if (data === null) return;
        const newTotalTrace = {
            x: Object.keys(data[props.database]),
            y: Object.values(data[props.database]).map((e) => {
                return e.totalChair + e.totalNonChair;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'red'},
            name: 'Total Sugars',
        };

        setTotalTrace(newTotalTrace);

        const newErrorTrace = {
            x: Object.keys(data[props.database]),
            y: Object.values(data[props.database]).map((e) => {
                return e.totalNonChair;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'green'},
            name: 'Non Chair Sugars',
        };

        setErrorTrace(newErrorTrace);

        const newRelativeTrace = {
            x: Object.keys(data[props.database]),
            y: Object.values(data[props.database]).map((e) => {
                return (100 * e.totalErrors) / e.totalGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            yaxis: 'y2',
            // marker: {color: 'green'},
            name: 'Relative Validation Errors ',
        };

        setRelativeTrace(newRelativeTrace);
    }, [data, props.database]);

    return (
        <>
            {data === null ? (
                <Loading loadingText={'Crunching latest data...'} />
            ) : (
                <Plot
                    data={[
                        totalTrace,
                        errorTrace,
                        // relativeTrace
                    ]}
                    layout={{
                        autosize: true,
                        width,
                        height,
                        title: {
                            text: `Conformational Anomalies in ${
                                props.database === 'pdbredo'
                                    ? 'PDB-REDO'
                                    : 'the PDB'
                            } over time`,
                            x: 0.5,
                            // y: 1.1,
                            xanchor: 'auto', // or 'auto', which matches 'left' in this case
                            yanchor: 'bottom',
                            xref: 'paper',
                            yref: 'paper',
                        },
                        plot_bgcolor: '#FFFFFF',
                        paper_bgcolor: 'rgba(0,0,0,0)',
                        modebar: {
                            bgcolor: 'rgba(0,0,0,0)',
                            color: 'gray',
                            activecolor: 'black',
                        },
                        yaxis: {
                            title: {
                                text: 'Number',
                            },
                            tickformat: ',.0f',
                            linewidth: 2,
                            mirror: true,
                            automargin: true,
                            ticksuffix: ' ',
                            tickprefix: '    ',
                            range: [-50, 35000],
                        },
                        yaxis2: {
                            overlaying: 'y',
                            tickformat: ',.0f',
                            tickmode: 'auto',
                            // anchor: 'free',
                            side: 'right',
                            automargin: true,
                            tickprefix: '  ',
                            range: [0, 100],
                        },
                        xaxis: {
                            title: {
                                text: 'Year',
                            },
                            linecolor: 'black',
                            linewidth: 2,
                            mirror: true,
                            tickmode: 'auto',
                            range: [1980, d.getFullYear()],
                        },

                        legend: {
                            x: legendDown ? 0 : 1.15,
                            y: legendDown ? -0.6 : 0.5,
                            // bgcolor: '#FFFFFF',
                        },
                    }}
                    config={{
                        toImageButtonOptions: {
                            format: 'svg',
                            filename: 'validationErrorsOverTime',
                            height: 1000,
                            width: 1500,
                            scale: 1,
                        },
                    }}
                />
            )}
        </>
    );
}
