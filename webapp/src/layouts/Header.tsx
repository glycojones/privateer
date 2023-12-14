import React, { lazy, type ReactElement, Suspense } from 'react';
import SNFG from '../components/PrivateerDisplay/SNFG.tsx';

import { type HeaderProps } from '../interfaces/types';
import { MoorhenReduxProvider } from 'moorhen';

const Upload = lazy(
    async () => await import('../components/Upload/Upload.tsx')
);
const Loading = lazy(
    async () => await import('../components/Loading/Loading.tsx')
);
const NavBar = lazy(async () => await import('./NavBar.tsx'));
const NoGlycans = lazy(
    async () => await import('../components/NoGlycans/NoGlycans.tsx')
);

export function Header(props: HeaderProps): ReactElement {
    let filename = '';
    if (props.PDBCode !== '') {
        filename = props.PDBCode;
    } else if (props.coordinateFile !== null) {
        filename = props.coordinateFile.name;
    }

    return (
        <div className="bg-gray text-primary text-center">
            <NavBar setResetApp={props.setResetApp} />

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
                            <SNFG {...props} filename={filename}></SNFG>
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
