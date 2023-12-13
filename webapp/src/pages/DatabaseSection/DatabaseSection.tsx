import React, { lazy, type ReactElement, useEffect, useState } from "react";
import { Information } from "../../components/Information/Information.tsx";
import { type DatabaseHeaderProps } from "../../interfaces/types";
import { DatabaseHeader } from "../../layouts/DatabaseHeader.tsx";

const Footer = lazy(async () => await import("../../layouts/Footer.tsx"));
const BorderElement = lazy(
  async () => await import("../../layouts/BorderElement.tsx"),
);

export default function DatabaseSection(props: {
  query: any;
  setSearchParams: any;
}): ReactElement {
  const [PDBCode, setPDBCode] = useState<string>("");
  const [submit, setSubmit] = useState<boolean>(false);
  const [loadingText, setLoadingText] = useState<string>(
    "Validating Glycans...",
  );
  const [resetApp, setResetApp] = useState<boolean>(false);
  const [fallback, setFallBack] = useState<boolean>(false);
  const [failureText, setFailureText] = useState<string>("");
  const [results, setResults] = useState<string>("");
  // const [failure, setFailure] = useState<boolean>(false);

  async function handleDatabaseLookup(PDBCode: string): Promise<void> {
    const pdbCode = PDBCode.toLowerCase();
    const middlefix = pdbCode.substring(1, 3);

    const url = `https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/${middlefix}/${pdbCode}.json`;

    try {
      const response = await fetch(url);
      const data: string = await response.json();
      setResults(data);
    } catch {
      console.log("not in db");
      setFallBack(true);
      setFailureText("This PDB is not in the database");
    }
  }

  useEffect(() => {
    if (PDBCode !== "") {
      setLoadingText(`Fetching ${PDBCode.toUpperCase()} from the database`);
      handleDatabaseLookup(PDBCode).then(
        () => {},
        () => {},
      );
      props.query.set("pdb", PDBCode);
      props.setSearchParams({ pdb: PDBCode });
    }
  }, [submit]);

  useEffect(() => {
    setSubmit(false);
    setFallBack(false);
    setResetApp(false);
    setPDBCode("");
    setResults("");
    props.setSearchParams({});
  }, [resetApp]);

  useEffect(() => {
    if (props.query.get("pdb") != null) {
      const pdb: string = props.query.get("pdb");
      handleDatabaseLookup(pdb).then(
        () => {},
        () => {},
      );
      setPDBCode(pdb);
    }
  }, []);

  const mainProps: DatabaseHeaderProps = {
    resetApp,
    setResetApp,
    PDBCode,
    setPDBCode,
    submit,
    setSubmit,
    loadingText,
    fallback,
    failureText,
    results,
  };

  return (
    <>
      <DatabaseHeader {...mainProps} />
      <BorderElement
        topColor={"#D6D9E5"}
        bottomColor={"#F4F9FF"}
        reverse={false}
      ></BorderElement>
      <Information />
      <BorderElement
        topColor={"#F4F9FF"}
        bottomColor={"#D6D9E5"}
        reverse={true}
      ></BorderElement>
      <Footer></Footer>
    </>
  );
}
