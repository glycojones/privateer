import { useMemo, Suspense } from 'react'
import './App.css'
import HomeSection from './pages/Home/HomeSection'
import DatabaseSection from './pages/DatabaseSection/DatabaseSection'

import PageLoad from './components/Loading/PageLoad'
import { Routes, Route, useSearchParams } from "react-router-dom";
import {
  BrowserRouter as Router,
  Link,
  useLocation
} from "react-router-dom";
import APIForwarding from './components/APIComponent/APIForwarding'

function useQuery() {
  const { search } = useLocation();
  return useMemo(() => new URLSearchParams(search), [search]);
}


function App() {
  let query = useQuery();
  const [searchParams, setSearchParams] = useSearchParams();

  return (
    <Suspense fallback={<PageLoad />}>
      <div className='flex flex-col'>

        <Routes>
          <Route index path="/" element={<HomeSection />} />
          <Route path="/database" element={<DatabaseSection query={query} setSearchParams={setSearchParams} />} />
          <Route path="/api" element={<APIForwarding query={query}/>} />

        </Routes>
      </div>
    </Suspense>
  )
}

export default App