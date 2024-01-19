import React, { useMemo, Suspense } from 'react';
import './App.css';
import HomeSection from './pages/Home/HomeSection';
import DatabaseSection from './pages/DatabaseSection/DatabaseSection';
import Statistics from './pages/Statistics/Statistics.tsx';

import PageLoad from './shared/Loading/PageLoad';
import { Routes, Route, useSearchParams, useLocation } from 'react-router-dom';
// import APIForwarding from "./components/APIComponent/APIForwarding";

function useQuery() {
    const { search } = useLocation();
    return useMemo(() => new URLSearchParams(search), [search]);
}

function App() {
    const query = useQuery();
    const [_, setSearchParams] = useSearchParams();

    return (
        <Suspense fallback={<PageLoad />}>
            <div id="main" className="flex flex-col">
                <Routes>
                    {/* @ts-expect-error */}
                    <Route index path="/" element={<HomeSection />} />
                    <Route
                        path="/database"
                        element={
                            <DatabaseSection
                                query={query}
                                setSearchParams={setSearchParams}
                            />
                        }
                    />
                    <Route path="/statistics" element={<Statistics />} />
                    {/* <Route path="/api" element={<APIForwarding query={query} />} /> */}
                </Routes>
            </div>
        </Suspense>
    );
}

export default App;
