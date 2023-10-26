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
    resetApp,
    setResetApp,
    PDBCode, 
    setPDBCode,
    coordinateFile,
    setCoordinateFile,
    reflectionFile,
    setReflectionFile,
    submit,
    setSubmit,
    tableData,
    loadingText,
    fileContent,
    fallback, 
    mtzData,
    failureText
}) {
    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={setResetApp}/>
            <div className="flex justify-center mb-6">
                {fallback != true ?
                    <Suspense fallback={<Loading loadingText={"Loading"} />}>
                        {submit == null ?
                            <Upload PDBCode={PDBCode} setPDBCode={setPDBCode} coordinateFile={coordinateFile} setCoordinateFile={setCoordinateFile} reflectionFile={reflectionFile} setReflectionFile={setReflectionFile} submitPressed={setSubmit} resetApp={resetApp} setResetApp={setResetApp} />
                                : tableData == null ?
                                    <Loading loadingText={loadingText} /> :
                                    <SNFG tableData={tableData} fileName={PDBCode != "" ? PDBCode : coordinateFile.name} PDBCode={PDBCode} pdbString={fileContent} mtzData={mtzData}></SNFG>}
                    </Suspense>
                    : <NoGlycans setResetApp={setResetApp} text={failureText} />
                } </div>
        </div>);
}
