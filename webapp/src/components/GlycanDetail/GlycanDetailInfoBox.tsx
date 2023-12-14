import React from 'react';

interface GlycanDetailInfoBoxProps {
    key: string;
    row: any;
}
export default function GlycanDetailInfoBox({
    key,
    row,
}: GlycanDetailInfoBoxProps) {
    return (
        <div key={key} className="flex flex-col flex-wrap px-8">
            <h3 className="text-left text-xl font-bold mt-2">
                Validation Report
            </h3>
            <div className="flex  justify-between mt-3">
                <h4>
                    Glycan ID: <b>{row.id}</b>
                </h4>
                <h4>
                    GlyTouCan ID:{' '}
                    <a
                        href={`https://glytoucan.org/Structures/Glycans/${row.glytoucan_id}`}
                    >
                        <b>{row.glytoucan_id}</b>
                    </a>
                </h4>
            </div>
            <div className="flex flex-col mt-3">
                <h4>
                    Number of conformation issues: <b>{row.conformation_err}</b>
                </h4>
                <h4>
                    Number of anomer issues: <b>{row.anomer_err}</b>
                </h4>
                <h4>
                    Number of torsion issues: <b>{row.torsion_err}</b>
                </h4>
                <h4>
                    Number of pucker issues: <b>{row.puckering_err}</b>
                </h4>
                <h4>
                    Number of chirality issues: <b>{row.chirality_err}</b>
                </h4>
            </div>
        </div>
    );
}
