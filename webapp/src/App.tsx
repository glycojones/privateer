import React, { useMemo, Suspense } from 'react';
import './App.css';
import Home from './routes/Home/Home.tsx';
import Database from './routes/Database/Database.tsx';
import Statistics from './routes/Statistics/Statistics.tsx';
import PageLoad from './shared/Loading/PageLoad';
import { Routes, Route, useSearchParams, useLocation } from 'react-router-dom';

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
                    <Route index path="/" element={<Home />} />
                    <Route
                        path="/database"
                        element={
                            <Database
                                query={query}
                                setSearchParams={setSearchParams}
                            />
                        }
                    />
                    <Route path="/statistics" element={<Statistics />} />
                </Routes>
            </div>
        </Suspense>
    );
}

export default App;
