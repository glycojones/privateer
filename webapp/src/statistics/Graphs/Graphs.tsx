import GlycansVsYear from './GlycansVsYear.tsx';
import React from 'react';
export default function Graphs() {
    return (
        <>
            <div className="flex flex-col items-center justify-center">
                <h2 className="mx-auto">Statistics</h2>

                <GlycansVsYear />
            </div>
        </>
    );
}
