import React, { useEffect, useState } from 'react';
import Loading from '../../shared/Loading/Loading.tsx';
import Plot from 'react-plotly.js';

export default function GlycansVsYear() {
    const [totalTrace, setTotalTrace] = useState();
    const [nglycanTrace, setNGlycanTrace] = useState();
    const [oglycanTrace, setOGlycanTrace] = useState();
    const [cglycanTrace, setCGlycanTrace] = useState();
    const [sglycanTrace, setSGlycanTrace] = useState();
    const [depositedTrace, setDepositedTrace] = useState();

    const [data, setData] = useState<Record<
        string,
        Record<string, number>
    > | null>(null);

    const [width, setWidth] = useState(1000);
    const [height, setHeight] = useState(1000);
    const [legendDown, setLegendDown] = useState(false);

    useEffect(() => {
        const url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/stats/glycosylation_per_year.json';

        fetch(url)
            .then(async (response) => await response.json())
            .then((data) => {
                setData(data as Record<string, Record<string, number>>);
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
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.totalGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            marker: { color: 'red' },
            name: 'Total Carbohydrates',
            line: {
                width: 3,
            },
        };

        setTotalTrace(newTotalTrace);

        const newNGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.nGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            marker: { color: 'blue' },
            name: 'N-Glycans',
            line: {
                width: 3,
            },
        };

        setNGlycanTrace(newNGlycanTrace);

        const newOGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.oGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            marker: { color: 'gold' },
            name: 'O-Glycans',
            line: {
                width: 3,
            },
        };

        setOGlycanTrace(newOGlycanTrace);

        const newSGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.sGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'gold'},
            name: 'S-Glycans',
            line: {
                width: 3,
            },
        };

        setSGlycanTrace(newSGlycanTrace);

        const newCGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.cGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'blue'},
            name: 'C-Glycans',
            line: {
                width: 3,
            },
        };

        setCGlycanTrace(newCGlycanTrace);

        const newDepositedTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.totalDepositions;
            }),
            type: 'bar',
            name: 'Total Deposited',
            yaxis: 'y2',
            marker: {
                color: 'rgb(100,100,100)',
                opacity: 0.2,
            },
        };

        setDepositedTrace(newDepositedTrace);
    }, [data]);

    return (
        <>
            {data === null ? (
                <Loading loadingText={'Crunching latest data...'} />
            ) : (
                <Plot
                    data={[
                        depositedTrace,
                        totalTrace,
                        nglycanTrace,
                        oglycanTrace,
                        sglycanTrace,
                        cglycanTrace,
                        // ligandTrace,
                    ]}
                    layout={{
                        autosize: true,
                        width,
                        height,
                        title: {
                            text: '<b>Glycosylation in the PDB over time</b>',
                            x: 0.5,
                            // y: 1.1,
                            xanchor: 'auto', // or 'auto', which matches 'left' in this case
                            yanchor: 'bottom',
                            xref: 'paper',
                            font: {
                                size: width < 800 ? 18 : 24,
                                family: 'sans-serif',
                            },
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
                                text: 'Number of Depositions',
                                font: {
                                    size: 18,
                                    family: 'sans-serif',
                                },
                            },
                            tickformat: ',.0f',

                            linewidth: 2,
                            mirror: true,
                            automargin: true,
                            ticksuffix: ' ',
                            tickprefix: '    ',
                            range: [-50, 1800],
                        },
                        yaxis2: {
                            title: {
                                text: 'Total Number of Depositions',
                                font: {
                                    size: 18,
                                    family: 'sans-serif',
                                },
                            },
                            overlaying: 'y',
                            tickformat: ',.0f',
                            tickmode: 'auto',
                            // anchor: 'free',
                            side: 'right',
                            automargin: true,
                            tickprefix: '  ',
                            ticksuffix: '   ',
                            range: [-500, 18000],
                        },
                        xaxis: {
                            title: {
                                text: 'Year',
                                font: {
                                    size: 18,
                                    family: 'sans-serif',
                                },
                            },
                            linecolor: 'black',
                            linewidth: 2,
                            mirror: true,
                            tickmode: 'auto',
                        },

                        legend: {
                            x: legendDown ? 0 : 1.15,
                            y: legendDown ? -0.8 : 0.5,
                            // bgcolor: '#FFFFFF',
                            font: {
                                size: 14,
                            },
                        },
                    }}
                    config={{
                        toImageButtonOptions: {
                            format: 'png',
                            filename: 'glycosylationOverTime',
                            height: 800,
                            width: 1500,
                            scale: 1,
                        },
                    }}
                />
            )}
        </>
    );
}
