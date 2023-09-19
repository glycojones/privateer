import { useEffect, useState } from "react";
import TorsionPlot from "./TorsionPlot";


function sortTorsions(torsions) { 
    const linkage_set = new Set();

    torsions.map((torsion) => {
        let linkage_string = torsion.sugar_1 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_2
        linkage_set.add(linkage_string)
     })

    const linkage_array = Array.from(linkage_set)
    
    const sorted_linkage_array = {}

    linkage_array.map((item) => {
        sorted_linkage_array[item] = []
     })

    torsions.map((torsion) => {
        let linkage_string = torsion.sugar_1 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_2
        sorted_linkage_array[linkage_string].push({"phi": torsion.phi, "psi": torsion.psi})
    })

    return [linkage_array, sorted_linkage_array]
}

function TorsionMultiPlotTabs({torsions, setTab}) { 

    const [linkage_array, sorted_linkage_array] = sortTorsions(torsions)

    console.log(linkage_array, sorted_linkage_array)

    return (
       linkage_array.map((item, index) => { 
        return (
            <div className="px-5"> 
                <button className="bg-tertiary text-sec" onClick={
                    () => {
                        setTab(index)}
                        }> <p className="p-1">{item}</p> </button>
            </div>
        )
       })
    )
}

export default function TorsionMultiPlot({torsions}) { 

    const [tab, setTab] = useState(0)

    const [linkage_array, sorted_linkage_array] = sortTorsions(torsions)
    console.log(linkage_array, sorted_linkage_array)

    useEffect(() => {console.log("TAB UODATED", tab)}, [tab])
    
    return (
        <>  
            <div className="flex">
                <TorsionMultiPlotTabs torsions={torsions} setTab={setTab}/>   
            </div>

            <TorsionPlot linkage_type={linkage_array[tab]} sorted_torsion_list={sorted_linkage_array}/>
        </>
    )
}