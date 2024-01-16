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
                return e.totalglyco;
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
                return e['n-glycosylation'];
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
                return e['o-glycosylation'];
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
                return e['s-glycosylation'];
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
                return e['c-glycosylation'];
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
                return e.depositions;
            }),
            type: 'bar',
            name: 'Deposited',
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
                title: 'Glycosylation in the PDB over time',
                plot_bgcolor: '#FFFFFF',
                paper_bgcolor: '#D6D9E5',
                yaxis: {
                    title: {
                        text: 'Number of Depositions',
                    },
                    tickformat: ',.0f',
                },
                yaxis2: {
                    side: 'right',
                    overlaying: 'y',
                    tickformat: ',.0f',
                },
                xaxis: {
                    title: {
                        text: 'Year',
                    },
                },
                legend: {
                    x: 1.1,
                    y: 0.5,
                    // bgcolor: '#FFFFFF',
                },
            }}
        />
    );
}
