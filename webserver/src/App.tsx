import { useEffect, useState, Suspense } from 'react'
import './App.css'
import HomeSection from './pages/Home/HomeSection'
import DatabaseSection from './pages/DatabaseSection/DatabaseSection'

import PageLoad from './components/Loading/PageLoad'
import {  Routes, Route } from "react-router-dom";

function App() {
  return (
    <Suspense fallback={<PageLoad />}>
      <div className='flex flex-col'>

          <Routes>
            <Route index path="/" element={<HomeSection />}/>
            <Route path="/database" element={<DatabaseSection />} />
          </Routes>
      </div>
    </Suspense>
  )
}

export default App