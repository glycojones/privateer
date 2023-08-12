import './App.css';
import {useEffect} from "react"
import privateer_module from '../wasm/privateer'; 

function App() {
  useEffect(() => {

    privateer_module().then((module_) => {
      console.log("WASM LOADED");
      
    });
  })

  return (
    <div className="App">
      <h1>Hello World!</h1>
    </div>
  );
}

export default App;
