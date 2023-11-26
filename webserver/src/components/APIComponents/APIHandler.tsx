import { useEffect, useState } from "react"
import { Link } from "react-router-dom"

export default function APIHandler(props) {

    const [redirectString, setRedirectString] = useState("/database")
    const [doRedirect, setDoRedirect] = useState(false)

    useEffect(() => {
        console.log(props)
        if (props.query.get("pdb") != null) {
            let pdb = props.query.get("pdb")
            let middle_fix = pdb.substring(1, 3)

            console.log(pdb)

            let redirect_string = `https://dialpuri.github.io/PrivateerDatabase/${middle_fix}/${pdb}.json.gz`
            window.location.replace(redirect_string);

            setRedirectString(redirect_string)
            console.log(redirect_string)
            setDoRedirect(true)
        }
    }, [])

    return (<>
        {doRedirect ?
            <Link to={redirectString} /> :
            <>Nothing to see</>
        }
    </>
    )

}