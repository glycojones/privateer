import GlycansVsYear from './GlycansVsYear.tsx';
import React, { Suspense, useEffect, useState } from 'react';
import Loading from '../../shared/Loading/Loading.tsx';
import ErrorsVsYear from './old/ErrorsVsYear.tsx';
import ConformationErrorsVsYear from './ConformationErrorsVsYear.tsx';

import ErrorsVsResolution from './old/ErrorsVsResolution.tsx';
import PDBToggle from '../../shared/PDBToggle/PDBToggle.tsx';
import BorderElement from '../../layouts/BorderElement.tsx';
import ConformationalErrorsVsResolution from './ConformationalErrorsVsResolution.tsx';

export default function Graphs() {
    const [lastUpdated, setLastUpdated] = useState<string>();

    useEffect(() => {
        const url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/stats/last_updated.json';

        fetch(url)
            .then(async (response) => await response.json())
            .then((data) => {
                setLastUpdated(data.date as string);
            })
            .catch((error) => {
                console.log('Error reading last updated date', error);
            });
    }, []);

    const [confVsResDBSwitch, setConfVsResDBSwitch] = useState(false);
    const [confVsYearDBSwitch, setConfVsYearDBSwitch] = useState(false);
    const [errorVsYearDBSwitch, setErrorVsYearDBSwitch] = useState(false);
    const [errorVsResDBSwitch, setErrorVsResDBSwitch] = useState(false);

    return (
        <>
            <Suspense
                fallback={<Loading loadingText={'Getting latest data'} />}
            >
                <div className="flex flex-col items-center justify-center">
                    <h2 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                        Protein Data Bank Statistics
                    </h2>
                    <h4 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                        Last Updated - {lastUpdated}
                    </h4>

                    <GlycansVsYear />

                    <BorderElement
                        bottomColor={'#F4F9FF'}
                        topColor={'#D6D9E5'}
                        reverse={false}
                    ></BorderElement>

                    <div className="flex flex-col sm:flex-row flex-wrap sm:justify-between w-full items-center align-center bg-tertiary">
                        <div>
                            <h2 className="w-full text-center sm:text-left pl-2 mt-8 sm:pl-12">
                                Validation Statistics
                            </h2>
                            <h4 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                                Last Updated - {lastUpdated}
                            </h4>
                        </div>
                    </div>

                    {/*<div className="w-full bg-tertiary text-center">*/}
                    {/*    <div className="w-full sm:px-12 text-center sm:text-left bg-tertiary">*/}
                    {/*        <PDBToggle*/}
                    {/*            checkState={errorVsYearDBSwitch}*/}
                    {/*            setCheckState={setErrorVsYearDBSwitch}*/}
                    {/*        />*/}
                    {/*    </div>*/}
                    {/*    <ErrorsVsYear*/}
                    {/*        database={errorVsYearDBSwitch ? 'pdbredo' : 'pdb'}*/}
                    {/*    />*/}
                    {/*</div>*/}

                    {/*<BorderElement*/}
                    {/*    topColor={'#F4F9FF'}*/}
                    {/*    bottomColor={'#D6D9E5'}*/}
                    {/*    reverse={true}*/}
                    {/*></BorderElement>*/}
                    <div className="w-full text-center bg-tertiary">
                        <div className="w-full sm:px-12 text-center sm:text-left ">
                            <PDBToggle
                                checkState={confVsYearDBSwitch}
                                setCheckState={setConfVsYearDBSwitch}
                            />
                        </div>

                        <ConformationErrorsVsYear
                            database={confVsYearDBSwitch ? 'pdbredo' : 'pdb'}
                        />
                    </div>

                    {/*<BorderElement*/}
                    {/*    bottomColor={'#F4F9FF'}*/}
                    {/*    topColor={'#D6D9E5'}*/}
                    {/*    reverse={false}*/}
                    {/*></BorderElement>*/}

                    {/*<div className="w-full bg-tertiary text-center">*/}
                    {/*    <div className="w-full sm:px-12 text-center sm:text-left">*/}
                    {/*        <PDBToggle*/}
                    {/*            checkState={errorVsResDBSwitch}*/}
                    {/*            setCheckState={setErrorVsResDBSwitch}*/}
                    {/*        />*/}
                    {/*    </div>*/}
                    {/*    <ErrorsVsResolution*/}
                    {/*        database={errorVsResDBSwitch ? 'pdbredo' : 'pdb'}*/}
                    {/*    />*/}
                    {/*</div>*/}

                    <BorderElement
                        topColor={'#F4F9FF'}
                        bottomColor={'#D6D9E5'}
                        reverse={true}
                    ></BorderElement>

                    <div className="w-full  text-center">
                        <div className="w-full sm:px-12 text-center sm:text-left">
                            <PDBToggle
                                checkState={confVsResDBSwitch}
                                setCheckState={setConfVsResDBSwitch}
                            />
                        </div>
                        <ConformationalErrorsVsResolution
                            database={confVsResDBSwitch ? 'pdbredo' : 'pdb'}
                        />
                    </div>
                </div>
            </Suspense>
        </>
    );
}
