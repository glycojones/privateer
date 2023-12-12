import { lazy, useEffect, useState } from "react";
import { Information } from '../../components/Information/Information.tsx';

const Footer = lazy(() => import('../../layouts/Footer.tsx'));
const BorderElement = lazy(() => import('../../layouts/BorderElement.tsx'));

import { DatabaseHeaderProps } from "../../interfaces/types"
import { DatabaseHeader } from "../../layouts/DatabaseHeader.tsx";

import pako from "pako"

export default function DatabaseSection(props) {

    const [PDBCode, setPDBCode] = useState<string>("")
    const [submit, setSubmit] = useState<boolean>(false);
    const [loadingText, setLoadingText] = useState<string>("Validating Glycans...");
    const [resetApp, setResetApp] = useState<boolean>(false)
    const [fallback, setFallBack] = useState<boolean>(false)
    const [failureText, setFailureText] = useState<string>("")
    const [results, setResults] = useState<string>("")
    const [failure, setFailure] = useState<boolean>(false)

    async function handle_database_lookup(PDBCode) {
        let pdb_code = PDBCode.toLowerCase()
        let middlefix = pdb_code.substring(1, 3)

        let url = `https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/${middlefix}/${pdb_code}.json`
        
        try {
            const response = await fetch(url)
            const data = await response.json()
            setResults(data)
        } catch {
            console.log("not in db")
            setFallBack(true)
            setFailureText("This PDB is not in the database")
        }
    }

    useEffect(() => {
        if (PDBCode != "") {
            setLoadingText(`Fetching ${PDBCode.toUpperCase()} from the database`)
            handle_database_lookup(PDBCode)
            props.query.set("pdb", PDBCode)
            props.setSearchParams({["pdb"]:PDBCode})
        }
    }, [submit])

    useEffect(() => {
        setSubmit(false)
        setFallBack(false)
        setResetApp(false)
        setPDBCode("")
        setResults("")
    }, [resetApp])


    useEffect(() => {
        if (props.query.get("pdb") != null) {
            handle_database_lookup(props.query.get("pdb"))
            setPDBCode(props.query.get("pdb"))
        }
    }, [])


    

    const main_props: DatabaseHeaderProps = {
        resetApp: resetApp,
        setResetApp: setResetApp,
        PDBCode: PDBCode,
        setPDBCode: setPDBCode,
        submit: submit,
        setSubmit: setSubmit,
        loadingText: loadingText,
        fallback: fallback,
        failureText: failureText,
        results: results,
    }

    return (
        <>
            <DatabaseHeader {...main_props} />
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"} reverse={false}></BorderElement>
            <Information />
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}