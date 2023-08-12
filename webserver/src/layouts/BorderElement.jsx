import styled from 'styled-components'
const Styles = styled.div`
div {
    height: 90px;
    background-image: linear-gradient(to bottom right, #eef2ff, #eef2ff 50%, #F4F9FF 50%, #F4F9FF);
  }
    `

export default function BorderElement({topColor, bottomColor, reverse}) {

    const divStyle = {
        height: "90px",
        backgroundImage: `linear-gradient(to bottom right, ${topColor}, ${topColor} 50%, ${bottomColor} 50%, ${bottomColor})`
      };

    if (reverse) { 
      divStyle['backgroundImage'] = `linear-gradient(to top right, ${bottomColor}, ${bottomColor} 50%, ${topColor} 50%, ${topColor})`
    }

    return (
        <div style={divStyle}></div>   
    )
}