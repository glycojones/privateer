
import "react-responsive-carousel/lib/styles/carousel.min.css"; 
// import { Carousel } from 'react-responsive-carousel';
import { useCallback, useState } from "react";

export default function SVGCarousel({svgs}) {

    let [index, setIndex] = useState(0);

    const next = (() => { 
        let newIndex = index + 1 
        if (newIndex >= svgs.length) { 
            newIndex = 0;
        }
        setIndex(newIndex);

    })

    const prev = (() => { 
        let newIndex = index - 1 
        if (newIndex < 0) { 
            newIndex = svgs.length-1
        }
        setIndex(newIndex);
    })

    const handleMoorhenClick = ((xyz) => {
        // parse_id()
        // sendToMoorhen()
    })

    const ref = useCallback((node) => {
        let useList = document.querySelectorAll('use')

        for (let i = 0; i < useList.length; i++) {
            useList[i].addEventListener("click", (e) => {console.log(useList[i].id)})
        }

        // document.getElementById('cercle1').setAttribute("height", "10px");
        // document.getElementsByClassName("svg")[0].setAttribute("width", )

    })

    return (
        <div className="flex flex-col items-center text-left">

            <div className="w-full">
                <h2 className="">Glycans</h2>
            </div>

            <div className="flex w-96 h-96 justify-center">
                {/* <SVG src ={svgs[index]}/> */}
                <div className="my-auto" dangerouslySetInnerHTML={{ __html: svgs[index]}} ref={ref} />

            </div>

            <div className="my-4">
                    <button onClick={prev} className="mx-4">Prev</button>
                    <span>{index+1}/{svgs.length}</span>
                    <button onClick={next} className="mx-4">Next</button>
                </div>
        </div>

);
}