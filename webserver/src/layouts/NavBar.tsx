import { Dispatch, SetStateAction } from "react"
import {GENERAL_CITATION, GITHUB_REPO} from "../data/Constants.tsx"

export default function NavBar({setResetApp}: {setResetApp: Dispatch<SetStateAction<boolean>>}) {

    return (
        <div className="flex flex-col sm:flex-row sm:justify-between">
        <div className="text-center sm:text-left px-12 pt-12 sm:p-12 flex flex-col">
            <span
                className="font-primary text-xl text-secondary sm:text-3xl">Validate your carbohydrates online with</span>
            <span className="font-body text-4xl sm:text-6xl my-2 sm:my-1"><button id="title" title="Home"
                                                                                  className="bg-gray"
                                                                                  onClick={() => {setResetApp(true)}}>Privateer</button></span>
            <span className="font-primary text-l text-secondary sm:text-xl sm:my-4 my-2">The Swiss Army knife for carbohydrate structure validation, refinement and analysis</span>
        </div>
        <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-12 flex items-center ">
            <a href={GITHUB_REPO}>
                <img className="w-full hover:scale-125 transition-all hidden dark:block" src="/github-mark.png"/>
                <img className="w-full hover:scale-125 transition-all block dark:hidden"
                     src="/github-mark-white.png"/>
            </a>
        </div>
    </div>
    )
}
    
  