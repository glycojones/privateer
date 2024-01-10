import React, { lazy, type ReactElement, Suspense } from 'react';
import { type DatabaseHeaderProps } from '../interfaces/types';

const Loading = lazy(async () => await import('../shared/Loading/Loading.tsx'));
const NavBar = lazy(async () => await import('./NavBar.tsx'));
const NoGlycans = lazy(
    async () => await import('../shared/NoGlycans/NoGlycans.tsx')
);
const DatabaseFetch = lazy(
    async () => await import('../database/DatabaseFetch/DatabaseFetch.jsx')
);
const DatabaseResult = lazy(
    async () => await import('../database/DatabaseResult/DatabaseResult.jsx')
);

export function DatabaseHeader(props: DatabaseHeaderProps): ReactElement {
    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={props.setResetApp} />
            <div className="flex justify-center mb-6">
                {!props.fallback ? (
                    <Suspense
                        fallback={
                            <Loading loadingText={'Loading Content...'} />
                        }
                    >
                        {props.pdbResults === '' ? (
                            <DatabaseFetch
                                PDBCode={props.PDBCode}
                                setPDBCode={props.setPDBCode}
                                submitPressed={props.setSubmit}
                            />
                        ) : (
                            <DatabaseResult
                                pdbResults={props.pdbResults}
                                pdbRedoResults={props.pdbRedoResults}
                                PDBCode={props.PDBCode}
                            />
                        )}
                    </Suspense>
                ) : (
                    <NoGlycans
                        setResetApp={props.setResetApp}
                        text={props.failureText}
                    />
                )}{' '}
            </div>
        </div>
    );
}
