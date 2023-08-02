import {MoorhenContextProvider, MoorhenMolecule, MoorhenContainer, itemReducer} from 'moorhen'
import SNFG from './SNFG'
import { useState, useReducer, useRef } from "react";


export default function Main() { 

    return (
            <div className="flex flex-grow bg-slate-700 justify-center items-center mb-5">
                <SNFG></SNFG>

                {/* <input onChange={handleFileSelected} type="file" /> */}

                {/* <MoorhenContextProvider>


                    <MoorhenContainer setMoorhenDimensions={(() => {return [200,400]})}/>
                    
                </MoorhenContextProvider> */}

            </div>
    )

}