import React, { useEffect, lazy, type ReactElement } from 'react';

const TorsionPlot = lazy(async () => await import('./TorsionPlot.tsx'));

function sortTorsions(torsions) {
    const linkageSet = new Set<string>();

    torsions.forEach(
        (torsion: {
            sugar_1: string;
            atom_number_2: string;
            atom_number_1: string;
            sugar_2: string;
        }): void => {
            let linkageString = '';
            if (torsion.sugar_1 === 'ASN') {
                linkageString =
                    torsion.sugar_1 +
                    '-' +
                    torsion.atom_number_2 +
                    ',' +
                    torsion.atom_number_1 +
                    '-' +
                    torsion.sugar_2;
            } else {
                linkageString =
                    torsion.sugar_2 +
                    '-' +
                    torsion.atom_number_2 +
                    ',' +
                    torsion.atom_number_1 +
                    '-' +
                    torsion.sugar_1;
            }
            linkageSet.add(linkageString);
        }
    );

    const linkgeArray: string[] = Array.from(linkageSet);

    const sortedLinkageMap = {};

    linkgeArray.forEach((item): void => {
        sortedLinkageMap[item] = [];
    });

    torsions.forEach(
        (torsion: {
            sugar_1: string;
            atom_number_2: string;
            atom_number_1: string;
            sugar_2: string;
            phi: any;
            psi: any;
        }): void => {
            let linkageString = '';
            if (torsion.sugar_1 === 'ASN') {
                linkageString =
                    torsion.sugar_1 +
                    '-' +
                    torsion.atom_number_2 +
                    ',' +
                    torsion.atom_number_1 +
                    '-' +
                    torsion.sugar_2;
            } else {
                linkageString =
                    torsion.sugar_2 +
                    '-' +
                    torsion.atom_number_2 +
                    ',' +
                    torsion.atom_number_1 +
                    '-' +
                    torsion.sugar_1;
            }

            sortedLinkageMap[linkageString].push({
                phi: torsion.phi,
                psi: torsion.psi,
            });
        }
    );

    return [linkgeArray, sortedLinkageMap];
}

function TorsionMultiPlotTabs({ torsions, setTab }): ReactElement {
    const [linkageArray, _] = sortTorsions(torsions);

    return linkageArray.map((item, index) => {
        return (
            <li className="mr-2" key={item}>
                <button
                    className="inline-block p-4 border-b-2 border-primary rounded-t-lg hover:scale-105"
                    onClick={() => {
                        setTab(index);
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

export default function TorsionMultiPlot({
    torsions,
    tab,
    setTab,
    size,
}: {
    torsions: any;
    tab: string;
    setTab: any;
    size: any;
}): ReactElement {
    const [linkageArray, sortedLinkageArray] = sortTorsions(torsions);

    useEffect(() => {
        setTab(0);
    }, []);

    return (
        <div className="flex flex-col align-middle justify-center items-center space-y-6 ">
            <div className="text-sm font-medium text-center text-gray-500 border-gray-200 text-gray-400 border-gray-700">
                <ul className="flex flex-wrap -mb-px mt-2 justify-center">
                    <TorsionMultiPlotTabs torsions={torsions} setTab={setTab} />
                </ul>
            </div>

            <TorsionPlot
                linkageType={linkageArray[tab]}
                sortedTorsionList={sortedLinkageArray}
                size={size}
            />
        </div>
    );
}
