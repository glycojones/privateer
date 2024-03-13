import React, { useEffect, useState, lazy, type ReactElement } from 'react';

import { linkageDB, binDB } from '../../data/Constants.tsx';
const Plot = lazy(async () => await import('react-plotly.js'));

export default function TorsionPlot({
    linkageType,
    sortedTorsionList,
    size,
    background,
}: {
    linkageType: string;
    sortedTorsionList: any;
    size: any;
    background: string;
}): ReactElement {
    const [trace, setTrace] = useState({});
    const [overlay, setOverlay] = useState({});
    const [linkageFound, setLinkageFound] = useState(false);

    useEffect(() => {
        if (sortedTorsionList === undefined) {
            return;
        }

        fetch(linkageDB[linkageType])
            .then(async (response) => await response.text())
            .then((text) => {
                return JSON.parse(text.replace(/\bNaN\b/g, 'null'));
            })
            .then((responseJson) => {
                const xData: number[] = [];
                const yData: number[] = [];

                for (let i = 0; i < responseJson.length; i++) {
                    xData.push(parseFloat(responseJson[i].phi as string));
                    yData.push(parseFloat(responseJson[i].psi as string));
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
                    colorbar: {
                        title: 'Frequency',
                        side: 'bottom',
                    },
                    xbins: {
                        start: -180,
                        end: 180,
                        size: 4,
                    },
                    autobiny: false,
                    ybins: binDB[linkageType],
                });
                setLinkageFound(true);
            })
            .catch((error) => {
                console.error(error);
                setLinkageFound(false);
            });

        const overlayPhi: number[] = [];
        const overlayPsi: number[] = [];

        for (let i = 0; i < sortedTorsionList[linkageType].length; i++) {
            overlayPhi.push(sortedTorsionList[linkageType][i].phi as number);
            overlayPsi.push(sortedTorsionList[linkageType][i].psi as number);
        }
        setOverlay({
            x: overlayPhi,
            y: overlayPsi,
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 8,
                color: 'blue',
                symbol: ['x'],
            },
        });
    }, [linkageType, sortedTorsionList]);

    return !linkageFound ? (
        <h3>
            {linkageType} does not have enough datapoints to generate a torsion
            plot.
        </h3>
    ) : (
        <Plot
            data={[trace, overlay]}
            layout={{
                width: size,
                height: size,
                title: linkageType,
                plot_bgcolor: '#FFFFFF',
                // paper_bgcolor: ,
                paper_bgcolor: background,
                yaxis: {
                    title: {
                        text: 'ψ / °',
                    },
                    fixedrange: true,
                    range:
                        linkageType in binDB
                            ? [binDB[linkageType].start, binDB[linkageType].end]
                            : [-180, 180],
                    showgrid: true,
                    showline: true,
                    zeroline: true,
                    zerolinewidth: 4,
                    zerolinecolor: '#969696',
                },
                xaxis: {
                    showline: true,
                    zeroline: true,

                    title: {
                        text: 'φ / °',
                    },
                    fixedrange: true,
                    range: [-180, 180],
                    showgrid: true,
                    zerolinewidth: 4,
                    zerolinecolor: '#969696',
                },
            }}
        />
    );
}
