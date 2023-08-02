import { useEffect, useState } from 'react'
import './App.css'
import NavBar from './components/navbar.jsx'
import Upload from './components/Upload'
import Main from './components/Main'
import Footer from './components/Footer'

function App() {
  const [count, setCount] = useState(0)

  // useEffect(() => {
  //   privateer_module().then((Module) => { 
  //     console.log("Privateer Module loaded!")
  //   });
  // })

  return (
    <div className='h-screen flex flex-col justify-between'>
      <NavBar></NavBar>
      <div className='mb-auto'>
        <Main/>
      </div>
      <Footer></Footer>
    </div>
  )
}

export default App
