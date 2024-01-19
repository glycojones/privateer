import { type Dispatch, type SetStateAction } from 'react';

export interface TorsionEntry {
    sugar_1: string;
    sugar_2: string;
    atom_number_1: string;
    atom_number_2: string;
    phi: number;
    psi: number;
}

export interface TableDataEntry {
    svg: string;
    wurcs: string;
    chain: string;
    glyconnect_id: string;
    glytoucan_id: string;
    id: string;
    torsion_err: number;
    conformation_err: number;
    anomer_err: number;
    puckering_err: number;
    chirality_err: number;
    torsions: TorsionEntry[];
}

export interface HeaderProps {
    resetApp: boolean;
    setResetApp: Dispatch<SetStateAction<boolean>>;
    PDBCode: string;
    setPDBCode: Dispatch<SetStateAction<string>>;
    coordinateFile: File | null;
    setCoordinateFile: Dispatch<SetStateAction<File | null>>;
    reflectionFile: File | null;
    setReflectionFile: Dispatch<SetStateAction<File | null>>;
    submit: boolean;
    setSubmit: Dispatch<SetStateAction<boolean>>;
    tableData: TableDataEntry[] | null;
    loadingText: string;
    fileContent: string | ArrayBuffer;
    fallback: boolean;
    mtzData: Uint8Array | null;
    failureText: string;
}

export interface DatabaseHeaderProps {
    resetApp: boolean;
    setResetApp: Dispatch<SetStateAction<boolean>>;
    PDBCode: string;
    setPDBCode: Dispatch<SetStateAction<string>>;
    submit: boolean;
    setSubmit: Dispatch<SetStateAction<boolean>>;
    loadingText: string;
    fallback: boolean;
    failureText: string;
    pdbResults: any;
    pdbRedoResults: any;
}

export interface StatisticsHeaderProps {
    resetApp: boolean;
    setResetApp: Dispatch<SetStateAction<boolean>>;
}

export interface DatabaseResultProps {
    pdbResults: any;
    pdbRedoResults: any;
    PDBCode: string;
}

export interface SNFGProps extends HeaderProps {
    filename: string;
}

export interface UploadButtonProps {
    setCoordinateFile: Dispatch<SetStateAction<File | null>>;
    setReflectionFile: Dispatch<SetStateAction<File | null>>;
}

export interface GlycanDetailProps {
    tableData: TableDataEntry[];
    hideMoorhen: boolean;
    setHideMoorhen: Dispatch<SetStateAction<boolean>>;
    rowID: number;
    controls: any;
    map: any;
    moorhenProps: any;
}

export interface NoGlycansProps {
    setResetApp: Dispatch<SetStateAction<boolean>>;
    text: string;
}

export interface SVGTableProps {
    tableData: TableDataEntry[];
    allowRowClick: boolean;
    rowClick: boolean;
    setRowClicked: Dispatch<SetStateAction<boolean>>;
    setRowID: Dispatch<SetStateAction<number>>;
}
