import GlycansVsYear from './GlycansVsYear.tsx';
import React, { useEffect, useState } from 'react';
export default function Graphs() {
    const [lastUpdated, setLastUpdated] = useState<string>();

    useEffect(() => {
        const url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/stats/last_updated.json';

        fetch(url)
            .then(async (response) => await response.json())
            .then((data) => {
                setLastUpdated(data.date as string);
            })
            .catch((error) => {
                console.log('Error reading last updated date', error);
            });
    }, []);

    return (
        <>
            <div className="flex flex-col items-center justify-center">
                <h2 className="w-full text-left pl-12">
                    Protein Data Bank Statistics
                </h2>
                <h4 className="w-full text-left pl-12">
                    Last Updated - {lastUpdated}
                </h4>

                <GlycansVsYear />
            </div>
        </>
    );
}
