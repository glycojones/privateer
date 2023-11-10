import { lazy, Suspense } from 'react'

import SNFG from '../components/PrivateerDisplay/SNFG.tsx'

const Upload = lazy(() => import('../components/Upload/Upload.tsx'));
const Loading = lazy(() => import('../components/Loading/Loading.tsx'));
const NavBar = lazy(() => import('./NavBar.tsx'));
const NoGlycans = lazy(() => import("../components/NoGlycans/NoGlycans.tsx"))

import {HeaderProps} from "../interfaces/types"

export function Header(props: HeaderProps) {
    
    let filename = ""
    if (props.PDBCode != "") { 
        filename = props.PDBCode
    }
    else if (props.coordinateFile) { 
        filename = props.coordinateFile.name
    }

    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={props.setResetApp}/>
            <div className="flex justify-center mb-6">
                {props.fallback !== true ?
                    <Suspense fallback={<Loading loadingText={"Loading Content..."} />}>
                        {props.submit === false ?
                            <Upload {...props} />
                                : props.tableData === null ?
                                    <Loading loadingText={props.loadingText} /> :
                                    <SNFG {...props} filename={filename}
                                    ></SNFG>}
                    </Suspense>
                    : <NoGlycans setResetApp={props.setResetApp} text={props.failureText} />
                } </div>
        </div>);
}
