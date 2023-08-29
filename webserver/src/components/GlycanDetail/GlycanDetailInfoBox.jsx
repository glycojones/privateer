

export default function GlycanDetailInfoBox({row}) { 
    
    
        
    return (<div className="">
    <h3 className="text-left text-xl">Validation Report</h3>
    <div className="flex justify-between mt-3"> 
        <h4>Glycan ID: <b>{row.id}</b></h4>
        <h4>Glytoucan ID: <b>{row.glytoucan_id}</b></h4>
    </div>
    <div className="flex flex-col mt-3"> 
        <h4>Number of conformation issues: <b>{row.conformation_err}</b></h4>
        <h4>Number of anomer issues:  <b>{row.anomer_err}</b></h4>
        <h4>Number of torsion issues:  <b>{row.torsion_err}</b></h4>
        <h4>Number of pucker issues:  <b>{row.puckering_err}</b></h4>
        <h4>Number of chirality issues:  <b>{row.chirality_err}</b></h4>

    </div>

    </div>)
}