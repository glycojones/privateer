import { useEffect, useState } from "react";

export default function PDBFetch({PDBCode, setPDBCode, submitPressed}) {

    const [pdb, setPDB] = useState("")

    useEffect( () => { 
        if (pdb.length != 4) { return }
        if (!pdb.match(/^[a-z0-9]+$/i)) { return }

        setPDBCode(pdb)
        submitPressed(true)
    }, [pdb])

    return (
        <>
        {PDBCode != true ? 
            <div className="flex items-center justify-center m-12 w-64 ">
            <label className="flex flex-col items-center justify-center w-full p-12 h-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer border-gray-600">
                <div className="flex flex-col items-center justify-center pt-5 pb-6 text-center" >
                    <svg className="w-6 h-6 mb-4 text-gray-500 dark:text-gray-400" xmlns="http://www.w3.org/2000/svg"
                        xmlnsXlink="http://www.w3.org/1999/xlink" version="1.1" id="Capa_1" x="0px" y="0px"
                        viewBox="0 0 56.966 56.966" xmlSpace="preserve"
                        width="512px" height="512px">
                        <path
                        d="M55.146,51.887L41.588,37.786c3.486-4.144,5.396-9.358,5.396-14.786c0-12.682-10.318-23-23-23s-23,10.318-23,23  s10.318,23,23,23c4.761,0,9.298-1.436,13.177-4.162l13.661,14.208c0.571,0.593,1.339,0.92,2.162,0.92  c0.779,0,1.518-0.297,2.079-0.837C56.255,54.982,56.293,53.08,55.146,51.887z M23.984,6c9.374,0,17,7.626,17,17s-7.626,17-17,17  s-17-7.626-17-17S14.61,6,23.984,6z" />
                    </svg>
                    <p className="mb-2 text-md text-gray-500 dark:text-gray-400"><span className="font-semibold">Fetch from PDB</span>
                    </p>
                    <input type="text" id="code" className="bg-gray-50 border border-gray-300 text-center text-gray-900 text-sm rounded-lg focus:border-3 block w-full p-2.5 my-2 " placeholder="5FJI" required/>
                    <button type="button" id="fetch" className="bg-gray-50 border font-bold  border-gray-300 text-gray-900 text-sm rounded-lg  focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                    onClick={
                        () => {
                            let element = document.getElementById("code")
                            setPDB(element.value)
                        }
                    }
                    >Fetch</button>

                </div>
            </label>
        </div> :
         <></>
        }
        </>
    )
}