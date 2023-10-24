import {GENERAL_CITATION, GITHUB_REPO} from "./Constants"

export default function NavBar() {

    return (
        <nav className="">
            <div className="h-16 bg-dark-900 border-bottom max-w-full mx-2 py-2 mb-2">
                <div className="flex justify-between px-4 py-2 first:mr-2 last:ml-2">
                    <div className="flex flex-1 justify-center">
                        <h1 className="text-3xl">Privateer</h1>
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
                                    <img className="h-6 hover:scale-125 transition-all hidden dark:block"
                                         src="github-mark-white.png"/>
                                    <img className="h-6 hover:scale-125 transition-all block dark:hidden"
                                         src="github-mark.png"/>
                                </a>
                            </li>
                        </ul>
                    </div>
                </div>
            </div>
        </nav>
    )
}
    
  