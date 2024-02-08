import React, { lazy, type ReactElement, Suspense } from 'react';
import PrivateerResults from '../PrivateerResults/PrivateerResults.tsx';

import { type HeaderProps } from '../../interfaces/types.ts';
import { MoorhenReduxProvider } from 'moorhen';

const Upload = lazy(async () => await import('../../shared/Upload/Upload.tsx'));
const Loading = lazy(
    async () => await import('../../shared/Loading/Loading.tsx')
);
const NavBar = lazy(async () => await import('../../layouts/NavBar.tsx'));
const NoGlycans = lazy(
    async () => await import('../../shared/NoGlycans/NoGlycans.tsx')
);

export function MainHeader(props: HeaderProps): ReactElement {
    let filename = '';
    if (props.PDBCode !== '') {
        filename = props.PDBCode;
    } else if (props.coordinateFile !== null) {
        filename = props.coordinateFile.name;
    }

    return (
        <div className="bg-gray text-primary text-center">
            <NavBar />

            {!props.fallback ? (
                <Suspense
                    fallback={<Loading loadingText={'Loading Content...'} />}
                >
                    {!props.submit ? (
                        <Upload {...props} />
                    ) : props.tableData === null ? (
                        <Loading loadingText={props.loadingText} />
                    ) : (
                        <MoorhenReduxProvider>
                            <PrivateerResults
                                {...props}
                                filename={filename}
                            ></PrivateerResults>
                        </MoorhenReduxProvider>
                    )}
                </Suspense>
            ) : (
                <NoGlycans
                    setResetApp={props.setResetApp}
                    text={props.failureText}
                />
            )}
        </div>
    );
}
