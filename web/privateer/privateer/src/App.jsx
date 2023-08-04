import { useEffect, useState } from 'react'
import './App.css'
import NavBar from './layouts/Navbar'
import Upload from './components/Upload'
import Main from './pages/Main/Main'
import Footer from './layouts/Footer'

function App() {
  const [resetApp, setResetApp] = useState(false)

  return (
    <div className='h-screen flex flex-col justify-between'>
      <NavBar setResetApp={setResetApp}></NavBar>
      <div className='mb-auto'>
        <Main resetApp={resetApp}/>
      </div>
      <Footer></Footer>
    </div>
  )
}

export default App
