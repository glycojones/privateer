import React, { lazy, type ReactElement } from 'react';
import Graphs from '../statistics/Graphs/Graphs.tsx';

const NavBar = lazy(async () => await import('./NavBar.tsx'));

export default function StatisticsHeader(): ReactElement {
    return (
        <div className="bg-gray text-primary">
            <NavBar />
            <Graphs />
        </div>
    );
}
