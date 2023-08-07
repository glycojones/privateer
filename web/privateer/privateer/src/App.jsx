import { useEffect, useState } from 'react'
import './App.css'
import HomeSection from './layouts/HomeSection'

function App() {
  const [resetApp, setResetApp] = useState(false)

  return (
    <div className='flex flex-col'>
      <HomeSection/>
    {/* //   <NavBar setResetApp={setResetApp}></NavBar>
    //   <div className='mb-auto'>
    //     <Main resetApp={resetApp}/>
    //   </div>
    //   <Footer></Footer> */}
     </div>
  )
}

export default App
