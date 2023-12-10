import { useEffect, lazy} from "react";

const TorsionPlot = lazy(() => import('./TorsionPlot'));


function sortTorsions(torsions) { 
    const linkage_set = new Set();

    torsions.map((torsion) => {
        let linkage_string = ""
        if (torsion.sugar_1 == "ASN") {
            linkage_string = torsion.sugar_1 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_2
        }
        else { 
            linkage_string = torsion.sugar_2 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_1
        }
        linkage_set.add(linkage_string)
     })

    const linkage_array = Array.from(linkage_set)
    
    const sorted_linkage_array = {}

    linkage_array.map((item) => {
        sorted_linkage_array[item] = []
     })

    torsions.map((torsion) => {
        let linkage_string = ""
        if (torsion.sugar_1 == "ASN") {
            linkage_string = torsion.sugar_1 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_2
        }
        else { 
            linkage_string = torsion.sugar_2 + "-" + torsion.atom_number_2 + "," + torsion.atom_number_1 + "-" + torsion.sugar_1
        }        
        
        sorted_linkage_array[linkage_string].push({"phi": torsion.phi, "psi": torsion.psi})
    })

    return [linkage_array, sorted_linkage_array]
}

function TorsionMultiPlotTabs({torsions, setTab}) { 

    const [linkage_array, sorted_linkage_array] = sortTorsions(torsions)

    return (
       linkage_array.map((item, index) => { 
        return (
            <li className="mr-2" key={item}>
                    <button className="inline-block p-4 border-b-2 border-transparent border-secondary rounded-t-lg hover:scale-105" onClick={() => {setTab(index)}}  onMouseDown={(e) => {e.stopPropagation()}}
            onTouchStart={(e) => {e.stopPropagation()}}>{item}</button>
            </li>

        )
       })
    )
}

export default function TorsionMultiPlot({torsions, tab, setTab}) { 


    const [linkage_array, sorted_linkage_array] = sortTorsions(torsions)

    useEffect(() => {
        setTab(0)
    }, [])
    
    return (
        <div className="flex flex-col align-middle justify-center items-center space-y-6 ">  
            <div className="text-sm font-medium text-center text-gray-500 border-gray-200 text-gray-400 border-gray-700">
                <ul className="flex flex-wrap -mb-px mt-2">                
                    <TorsionMultiPlotTabs torsions={torsions} setTab={setTab}/>   
                </ul>
            </div>

            <TorsionPlot linkage_type={linkage_array[tab]} sorted_torsion_list={sorted_linkage_array}/>
        </div>
    )
}