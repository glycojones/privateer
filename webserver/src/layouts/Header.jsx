import { lazy, Suspense } from 'react'

import SNFG from '../components/PrivateerDisplay/SNFG'
// import Upload from "../common/Upload";
// import Submit from "../common/Submit";
// import Loading from "../common/Loading";
import { GITHUB_REPO } from "../data/Constants"

const Upload = lazy(() => import('../common/Upload'));
const Submit = lazy(() => import('../common/Submit'));
const Loading = lazy(() => import('../common/Loading'));
const NavBar = lazy(() => import('../layouts/NavBar'));

const NoGlycans = lazy(() => import("../components/NoGlycans/NoGlycans"))

export function Header({
    setResetApp,
    coordinateFile,
    setCoordinateFile,
    reflectionFile,
    setReflectionFile,
    submit,
    setSubmit,
    tableData,
    loadingText,
    fileContent,
    fallback
}) {
    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={setResetApp}/>
            <div className="flex justify-center mb-6">
                {fallback != true ?
                    <Suspense fallback={<Loading loadingText={"Loading"} />}>
                        {submit == null ?
                            <Upload coordinateFile={coordinateFile} setCoordinateFile={setCoordinateFile} reflectionFile={reflectionFile} setReflectionFile={setReflectionFile} submitPressed={setSubmit} setResetApp={setResetApp} />
                            // : submit == null ?
                            //     <Submit coordinateFile={coordinateFile} submitPressed={setSubmit} setResetApp={setResetApp} />
                                : tableData == null ?
                                    <Loading loadingText={loadingText} /> :
                                    <SNFG tableData={tableData} fileName={coordinateFile.name} pdbString={fileContent}></SNFG>}
                    </Suspense>
                    : <NoGlycans setResetApp={setResetApp} />
                } </div>
        </div>);
}
