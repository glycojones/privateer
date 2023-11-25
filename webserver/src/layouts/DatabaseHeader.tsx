import { lazy, Suspense } from 'react'

const Upload = lazy(() => import('../components/Upload/Upload.tsx'));
const Loading = lazy(() => import('../components/Loading/Loading.tsx'));
const NavBar = lazy(() => import('./NavBar.tsx'));
const NoGlycans = lazy(() => import("../components/NoGlycans/NoGlycans.tsx"))
const DatabaseFetch = lazy(() => import("../components/DatabaseFetch/DatabaseFetch.jsx"))
const DatabaseResult = lazy(() => import("../components/DatabaseResult/DatabaseResult.jsx"))

import {DatabaseHeaderProps} from "../interfaces/types"

export function DatabaseHeader(props: DatabaseHeaderProps) {

    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={props.setResetApp}/>
            <div className="flex justify-center mb-6">
                {props.fallback !== true ?
                    <Suspense fallback={<Loading loadingText={"Loading Content..."} />}>
                        {!props.results ?
                            <DatabaseFetch PDBCode={props.PDBCode} setPDBCode={props.setPDBCode} submitPressed={props.setSubmit}/>
                            : <DatabaseResult results={props.results} PDBCode={props.PDBCode}/>
                            }
                    </Suspense>
                    : <NoGlycans setResetApp={props.setResetApp} text={props.failureText} />
                } </div>
        </div>);
}
