import React, { lazy, type ReactElement, Suspense } from 'react';
import { type StatisticsHeaderProps } from '../interfaces/types';
import Graphs from '../statistics/Graphs/Graphs.tsx';

const Loading = lazy(async () => await import('../shared/Loading/Loading.tsx'));
const NavBar = lazy(async () => await import('./NavBar.tsx'));

export default function StatisticsHeader(
    props: StatisticsHeaderProps
): ReactElement {
    return (
        <div className="bg-gray text-primary">
            <NavBar setResetApp={props.setResetApp} />
            <Graphs />
        </div>
    );
}
