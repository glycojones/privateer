import { useEffect, useState, Suspense } from 'react'
import './App.css'
import HomeSection from './pages/Home/HomeSection'

function App() {
  const [resetApp, setResetApp] = useState(false)

  return (
    <Suspense fallback={<div>Loading... </div>}>

    <div className='flex flex-col'>
      <HomeSection/>
     </div>
     </Suspense>
  )
}

export default App