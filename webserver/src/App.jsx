import { useEffect, useState } from 'react'
import './App.css'
import HomeSection from './pages/Home/HomeSection'

function App() {
  const [resetApp, setResetApp] = useState(false)

  return (
    <div className='flex flex-col'>
      <HomeSection/>
     </div>
  )
}

export default App
