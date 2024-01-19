import { type ReactElement } from 'react';
import React from 'react';
export default function BorderElement({
    topColor,
    bottomColor,
    reverse,
}: {
    topColor: string;
    bottomColor: string;
    reverse: boolean;
}): ReactElement {
    const divStyle = {
        height: '90px',
        backgroundImage: `linear-gradient(to bottom right, ${topColor}, ${topColor} 50%, ${bottomColor} 50%, ${bottomColor})`,
    };

    if (reverse) {
        divStyle.backgroundImage = `linear-gradient(to top right, ${bottomColor}, ${bottomColor} 50%, ${topColor} 50%, ${topColor})`;
    }

    return <div style={divStyle} className="w-full"></div>;
}
