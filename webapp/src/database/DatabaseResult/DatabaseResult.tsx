import React, { useEffect, useState } from 'react';
import { type DatabaseResultProps } from '../../interfaces/types.ts';
import CremerPopleGraph from '../DatabaseComponents/CremerPopleGraph.tsx';
import BFactorVsRSCC from '../DatabaseComponents/BFactorVsRSCC.tsx';
import SNFGList from '../DatabaseComponents/SNFGList.tsx';
import SugarList from '../DatabaseComponents/SugarList.tsx';

export default function DatabaseResult(props: DatabaseResultProps) {
    const [data, setData] = useState();
    useEffect(() => {
        if (props.results === null) return;

        setData(props.results as undefined);
    }, []);

    return (
        <>
            {data !== undefined ? (
                <div className="flex flex-col space-y-6">
                    <h2 className="text-center">
                        Validation Report - {props.PDBCode}
                    </h2>
                    <div className="flex flex-wrap text-center">
                        <CremerPopleGraph data={data} />
                        <BFactorVsRSCC data={data} />
                    </div>
                    <SNFGList data={data} />
                    <SugarList data={data} />
                </div>
            ) : (
                <></>
            )}
        </>
    );
}
