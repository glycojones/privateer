import { lazy, Suspense } from 'react'

import SNFG from '../components/PrivateerDisplay/SNFG'

const Upload = lazy(() => import('../components/Upload/Upload'));
const Loading = lazy(() => import('../components/Loading/Loading'));
const NavBar = lazy(() => import('../layouts/NavBar'));

const NoGlycans = lazy(() => import("../components/NoGlycans/NoGlycans"))

export function Header(props) {

    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={props.setResetApp}/>
            <div className="flex justify-center mb-6">
                {props.fallback != true ?
                    <Suspense fallback={<Loading loadingText={"Loading content..."} />}>
                        {props.submit == null ?
                            <Upload {...props} />
                                : props.tableData == null ?
                                    <Loading loadingText={props.loadingText} /> :
                                    <SNFG {...props} filename={props.PDBCode != "" ? props.PDBCode : props.coordinateFile.name}
                                    ></SNFG>}
                    </Suspense>
                    : <NoGlycans setResetApp={props.setResetApp} text={props.failureText} />
                } </div>
        </div>);
}
