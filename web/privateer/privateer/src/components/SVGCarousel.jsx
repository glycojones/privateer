
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
        let use_list = document.querySelectorAll('use')

        for (let i = 0; i < x.length; i++) { 
            use_list[i].addEventListener("click", (e) => {console.log(use_list[i].id)})
        }

    })

    return (
        <div className="flex flex-col items-center ">
            <div className="flex w-128 h-96 justify-center">
                {/* <SVG src ={svgs[index]}/> */}
                <div className="w-full my-auto" dangerouslySetInnerHTML={{ __html: svgs[index]}} ref={ref} />

            </div>

            <div className="my-4">
                    <button onClick={prev} className="mx-4">Prev</button>
                    <button onClick={next} className="mx-4">Next</button>
                </div>
        </div>

);
}