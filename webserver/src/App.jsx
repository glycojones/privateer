import { useEffect, useState, Suspense } from 'react'
import './App.css'
import HomeSection from './pages/Home/HomeSection'
import PageLoad from './components/Loading/PageLoad'

function App() {
  const [resetApp, setResetApp] = useState(false)

  return (
    <Suspense fallback={<PageLoad/>}>

    <div className='flex flex-col'>
      <HomeSection/>
     </div>
     </Suspense>
  )
}

export default App