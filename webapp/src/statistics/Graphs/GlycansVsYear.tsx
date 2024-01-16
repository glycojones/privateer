import React, { lazy, useEffect, useState } from 'react';
const Plot = lazy(async () => await import('react-plotly.js'));

export default function GlycansVsYear() {
    const [totalTrace, setTotalTrace] = useState();
    const [nglycanTrace, setNGlycanTrace] = useState();
    const [oglycanTrace, setOGlycanTrace] = useState();
    const [cglycanTrace, setCGlycanTrace] = useState();
    const [sglycanTrace, setSGlycanTrace] = useState();
    const [ligandTrace, setLigandTrace] = useState();
    const [depositedTrace, setDepositedTrace] = useState();

    const [data, setData] = useState<Record<
        string,
        Record<string, number>
    > | null>(null);

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
            // marker: {color: 'red'},
            name: 'Total Glycans',
        };

        setTotalTrace(newTotalTrace);

        const newNGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.nGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'green'},
            name: 'N-Glycans',
        };

        setNGlycanTrace(newNGlycanTrace);

        const newOGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.oGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'blue'},
            name: 'O-Glycans',
        };

        setOGlycanTrace(newOGlycanTrace);

        const newSGlycanTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.sGlycans;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'blue'},
            name: 'S-Glycans',
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
        };

        setCGlycanTrace(newCGlycanTrace);

        const newLigandTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.ligands;
            }),
            type: 'scatter',
            mode: 'lines',
            // marker: {color: 'blue'},
            name: 'Ligands',
        };

        setLigandTrace(newLigandTrace);

        const newDepositedTrace = {
            x: Object.keys(data),
            y: Object.values(data).map((e) => {
                return e.totalDepositions;
            }),
            type: 'bar',
            name: 'Total Deposited',
            yaxis: 'y2',
            marker: {
                color: 'rgb(158,202,225)',
                opacity: 0.4,
            },
        };

        setDepositedTrace(newDepositedTrace);
    }, [data]);

    return (
        <Plot
            data={[
                depositedTrace,
                totalTrace,
                nglycanTrace,
                oglycanTrace,
                sglycanTrace,
                cglycanTrace,
                ligandTrace,
            ]}
            layout={{
                title: {
                    text: 'Glycosylation in the PDB over time',
                    x: 0.5,
                    y: 1.1,
                    xanchor: 'auto', // or 'auto', which matches 'left' in this case
                    yanchor: 'bottom',
                    xref: 'paper',
                    yref: 'paper',
                },
                plot_bgcolor: '#FFFFFF',
                paper_bgcolor: '#D6D9E5',
                yaxis: {
                    title: {
                        text: 'Number of Depositions',
                    },
                    tickformat: ',.0f',
                    linewidth: 2,
                    mirror: true,
                    automargin: true,
                    ticksuffix: " ",
                    tickprefix: "    "

                },
                yaxis2: {
                    overlaying: 'y',
                    tickformat: ',.0f',
                    tickmode: 'auto',
                    // anchor: 'free',
                    side: 'right',
                    automargin: true,
                    tickprefix: "  ",
                    range: [0, 16000]
                },
                xaxis: {
                    title: {
                        text: 'Year',
                    },
                    linecolor: 'black',
                    linewidth: 2,
                    mirror: true,
                    tickmode: 'auto',

                },

                legend: {
                    x: 1.15,
                    y: 0.5,
                    // bgcolor: '#FFFFFF',
                },
            }}
        />
    );
}
