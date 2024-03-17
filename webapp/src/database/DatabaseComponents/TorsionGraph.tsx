import React, { useEffect, useState } from 'react';
import TorsionMultiPlot from '../../shared/TorsionPlot/TorsionMultiPlot.tsx';

function TorsionGraphTabs(props): React.JSX.Element[] {
    return props.keys.map((item, index) => {
        return (
            <li className="mr-2 text-lg" key={item}>
                <button
                    className={
                        props.tab === index
                            ? 'inline-block p-4 border-b-2 border-primary rounded-t-lg hover:scale-105 focus:outline-none'
                            : 'inline-block p-4 border-primary rounded-t-lg hover:scale-105 focus:outline-none'
                    }
                    onClick={() => {
                        props.setTab(index);
                    }}
                    onMouseDown={(e) => {
                        e.stopPropagation();
                    }}
                    onTouchStart={(e) => {
                        e.stopPropagation();
                    }}
                >
                    {item}
                </button>
            </li>
        );
    });
}

export default function TorsionGraph(data) {
    const [torsions, setTorsions] = useState<Record<any, any> | undefined>();
    const [torsionTab, setTorsionTab] = useState(0);
    const [glycanTab, setGlycanTab] = useState<number>(0);

    useEffect(() => {
        const glycans = data.data.glycans;
        const torsionList: Record<string, Array<Record<string, string>>> = {};

        for (const key in glycans) {
            const glycanType = glycans[key];
            for (let i = 0; i < glycanType.length; i++) {
                const chainID = glycanType[i].rootSugarChainId;
                const linkages = glycanType[i].linkages;

                for (const linkage in linkages) {
                    for (let j = 0; j < linkages[linkage].length; j++) {
                        const data = {
                            sugar_1: linkages[linkage][j].firstResidue,
                            sugar_2: linkages[linkage][j].secondResidue,
                            atom_number_1:
                                linkages[linkage][j].donorAtom.slice(-1),
                            atom_number_2:
                                linkages[linkage][j].acceptorAtom.slice(-1),
                            phi: linkages[linkage][j].phi,
                            psi: linkages[linkage][j].psi,
                        };
                        if (chainID in torsionList) {
                            torsionList[chainID].push(data);
                        } else {
                            torsionList[chainID] = [data];
                        }
                    }
                }
            }
        }
        // torsionList.sort;

        const sortedTorsionList = {};

        Object.keys(torsionList)
            .sort()
            .forEach(function (k, _) {
                sortedTorsionList[k] = torsionList[k];
            });

        if (Object.keys(sortedTorsionList).length === 0) {
            setTorsions(undefined);
        } else {
            setTorsions(sortedTorsionList);
        }
    }, [data]);

    useEffect(() => {
        setTorsionTab(0);
    }, [glycanTab]);

    return (
        <div className="flex flex-col w-full align-content-center">
            {torsions !== undefined ? (
                <>
                    <span className="text-xl">Linkage Torsion Analysis</span>
                    <div className="flex justify-center items-center ">
                        <span className="text-lg mt-2 mr-8">
                            <b>Chain:</b>
                        </span>
                        <ul className="flex flex-wrap -mb-px mt-2 justify-center  ">
                            <TorsionGraphTabs
                                keys={Object.keys(torsions)}
                                setTab={setGlycanTab}
                                tab={glycanTab}
                            />
                        </ul>
                    </div>

                    <TorsionMultiPlot
                        torsions={torsions[Object.keys(torsions)[glycanTab]]}
                        tab={torsionTab}
                        setTab={setTorsionTab}
                        size={500}
                        background={'#F4F9FF'}
                    />
                </>
            ) : (
                <>
                    <span className="text-xl">Linkage Torsion Analysis</span>
                    <div className="flex justify-center items-center text-center ">
                        <span className="text-lg mt-8 mb-6 text-center">
                            <b>There are no linkages detected in this model.</b>
                        </span>
                    </div>
                </>
            )}
        </div>
    );
}
