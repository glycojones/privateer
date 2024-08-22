import React, { lazy, type ReactElement, Suspense } from 'react';
import { type DatabaseHeaderProps } from '../../interfaces/types.ts';
const Loading = lazy(
    async () => await import('../../shared/Loading/Loading.tsx')
);
const NavBar = lazy(async () => await import('../../layouts/NavBar.tsx'));
const NoGlycans = lazy(
    async () => await import('../../shared/NoGlycans/NoGlycans.tsx')
);
const DatabaseInput = lazy(
    async () => await import('../DatabaseInput/DatabaseInput.jsx')
);
const DatabaseResult = lazy(
    async () => await import('../DatabaseResult/DatabaseResult.tsx')
);

export function DatabaseHeader(props: DatabaseHeaderProps): ReactElement {
    return (
        <div className="bg-gray text-primary">
            <NavBar />
            <div className="flex w-full justify-center mb-6">
                {!props.fallback ? (
                    <Suspense
                        fallback={
                            <Loading loadingText={'Loading Content...'} />
                        }
                    >
                        {props.pdbResults === '' ? (
                            <DatabaseInput
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