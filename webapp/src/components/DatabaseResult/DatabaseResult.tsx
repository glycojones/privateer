import React, { useEffect, useState } from 'react';
import { type DatabaseResultProps } from '../../interfaces/types';
import CremerPopleGraph from '../DatabaseComponents/CremerPopleGraph';
import BFactorVsRSCC from '../DatabaseComponents/BFactorVsRSCC';
import SNFGList from '../DatabaseComponents/SNFGList';
import SugarList from '../DatabaseComponents/SugarList';

export default function DatabaseResult (props: DatabaseResultProps) {
  const [data, setData] = useState();
  useEffect(() => {
    if (props.results === null) return;

    setData(props.results as undefined);
  }, []);

  return (
    <>
      {data !== undefined
        ? (
        <div className="flex flex-col space-y-6">
          <h2 className="text-center">Validation Report - {props.PDBCode}</h2>
          <div className="flex flex-wrap text-center">
            <CremerPopleGraph data={data} />
            <BFactorVsRSCC data={data} />
          </div>
          <SNFGList data={data} />
          <SugarList data={data} />
        </div>
          )
        : (
        <></>
          )}
    </>
  );
}
