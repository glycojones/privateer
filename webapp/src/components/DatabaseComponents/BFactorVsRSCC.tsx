import React, { useEffect, useState, lazy } from 'react';

const Plot = lazy(async () => await import('react-plotly.js'));

function calculatePoints (data): [number[], number[], string[]] {
  const glycans = data.data.glycans;

  const xAxis: number[] = [];
  const yAxis: number[] = [];
  const text: string[] = [];

  for (const key in glycans) {
    const glycanType = glycans[key];
    for (let i = 0; i < glycanType.length; i++) {
      const sugars = glycanType[i].Sugars;
      for (let j = 0; j < sugars.length; j++) {
        xAxis.push(sugars[j].BFactor as number);
        yAxis.push(sugars[j].RSCC as number);
        text.push(sugars[j]['Sugar ID'] as string);
      }
    }
  }

  return [xAxis, yAxis, text];
}

export default function BFactorVsRSCC (props) {
  const [trace, setTrace] = useState({});
  const [corrTrace, setCorrTrace] = useState({});

  useEffect(() => {
    const [xAxis, yAxis, text] = calculatePoints(props);

    const maxX = Math.max(...xAxis) + 5;

    setCorrTrace({
      x: [0, maxX],
      y: [0.7, 0.7],
      fill: 'tozeroy',
      fillcolor: 'rgba(173,181,189,0.3)',
      fillopacity: 0.1,
      marker: {
        size: 1,
        color: 'rgba(173,181,189,0.5)',
        symbol: ['o']
      }
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
        color: 'green',
        symbol: ['o']
      }
    });
  }, []);

  return (
    <div className="flex flex-col mx-auto">
      <span className="text-xl">BFactor vs RSCC</span>
      <Plot
        data={[trace, corrTrace]}
        layout={{
          showlegend: false,
          width: 500,
          height: 400,
          title: '',
          plot_bgcolor: '#FFFFF',
          paper_bgcolor: '#D6D9E5',
          margin: {
            l: 50,
            r: 50,
            b: 50,
            t: 10,
            pad: 4
          },
          yaxis: {
            title: {
              text: 'RSCC'
            },
            fixedrange: true,
            range: [0, 1],
            showgrid: true
          },
          xaxis: {
            title: {
              text: 'B Factor'
            },
            fixedrange: true,
            // range: [0, 100],
            showgrid: true
          }
        }}
      />
    </div>
  );
}
