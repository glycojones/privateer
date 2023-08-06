import { useState } from "react"
import { GITHUB_REPO, GENERAL_CITATION } from "../data/Constants"

export default function NavBar({setResetApp}) {

  return (
    <nav className="">
      <div className="h-16 bg-gray-800 shadow border-bottom max-w-full py-2 mb-2">
        <div className="flex justify-between px-4 py-2 first:mr-2 last:ml-2">
          <div className="flex flex-1 justify-center">
            <h1 className="text-3xl font-title"><button id="title" title="Home" onClick={() => {setResetApp(true)}}>Privateer</button></h1>
          </div>
          
          <div className="flex flex-1 justify-center"> 
            <ul className="flex items-center justify-between">
            <li className="pr-4 hover:scale-110">
              <a href={GENERAL_CITATION} className=" transition-all">
                Cite
                </a>
            </li>
            <li className="pl-4">
              <a href={GITHUB_REPO}>
                <img className="h-6 hover:scale-125 transition-all hidden dark:block" src="github-mark-white.png"/>
                <img className="h-6 hover:scale-125 transition-all block dark:hidden" src="github-mark.png"/>
              </a>
            </li>
            </ul>
            </div>
          </div>
        </div>
    </nav>
  )
}
    
  