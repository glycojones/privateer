
import "react-responsive-carousel/lib/styles/carousel.min.css"; 
// import { Carousel } from 'react-responsive-carousel';
import SVG from 'react-inlinesvg';
import { useState } from "react";

export default function SVGCarousel({svgs}) {

    let [index, setIndex] = useState(0);

    const next = (() => { 
        setIndex(index++);
        if (index >= svgs.length) { 
            setIndex(0);
        }
    })

    const prev = (() => { 
        setIndex(index--);
        if (index < 0) { 
            setIndex(svgs.length-1);
        }
    })

    return (
        <div className="flex flex-row">
            <button onClick={prev}>Prev</button>
            <SVG src ={svgs[index]}/>
            <button onClick={next}>Next</button>
        </div>

);
}