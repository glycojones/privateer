import { useEffect, useState } from 'react'
import './App.css'
import NavBar from './components/navbar.jsx'
import Upload from './components/Upload'
import Main from './components/Main'


function App() {
  const [count, setCount] = useState(0)

  // useEffect(() => {
  //   privateer_module().then((Module) => { 
  //     console.log("Privateer Module loaded!")
  //   });
  // })

  return (
    <div className='h-screen'>
      <NavBar></NavBar>
      <Main/>
    </div>
  )
}

export default App
