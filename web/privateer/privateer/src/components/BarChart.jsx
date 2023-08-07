import { useState, useEffect, useRef } from 'react';
import * as d3 from 'd3'

export default function BarChart() {

    const svgRef = useRef();

    useEffect(() => {
        
        var margin = { top: 20, right: 30, bottom: 30, left: 30 },
            width = 460 - margin.left - margin.right,
            height = 400 - margin.top - margin.bottom;

        // append the svg object to the body of the page
        var svg = d3.select("#data")
            .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .append("g")
            .attr("transform",
                "translate(" + margin.left + "," + margin.top + ")")
        let url = "https://raw.githubusercontent.com/glycojones/privateer/master/data/linkage_torsions/unprocessed_files/ASN-NAG_reduced.json"
        fetch(url).then(response => response.json()).then((jsonData) => {

            var xLim = [-180, 180]
            var yLim = [-180, 180]

            var x = d3.scaleLinear()
                .nice()
                .domain(xLim)
                .range([0, width]);
            svg.append("g")
                .attr("transform", "translate(0," + height + ")")
                .call(d3.axisBottom(x));

            var y = d3.scaleLinear()
                .nice()
                .domain(yLim)
                .range([height, 0]);
            svg.append("g")
                .call(d3.axisLeft(y));

            var inputForRectBinning = []


            for (let i = 0; i < jsonData.length; i++) { 
                inputForRectBinning.push([+jsonData[i].Phi, +jsonData[i].Psi])  // Note that we had the transform value of X and Y !
            }
            
            console.log(inputForRectBinning)

            


            (function() {
                d3.rectbin = function() {
                  var dx = 0.1,
                      dy = 0.1, 
                      x = rectbinX,
                      y = rectbinY;
              
                  function rectbin(points) {
                    var binsById = {};
                    var xExtent = d3.extent(points, function(d, i){ return x.call(rectbin, d, i) ;});
                    var yExtent = d3.extent(points, function(d, i){ return y.call(rectbin, d, i) ;});
              
                    d3.range(yExtent[0], yExtent[1] + dx, dy).forEach(function(Y){
                      d3.range(xExtent[0], xExtent[1] + dx, dx).forEach(function(X){
                        var py = Y / dy; 
                        var pj = trunc(py);
                        var px = X / dx;
                        var pi = trunc(px);
                        var id = pi + '-' + pj;
                        var bin = binsById[id] = [];
                        bin.i = pi;
                        bin.j = pj;
                        bin.x = pi * dx;
                        bin.y = pj * dy;
                      });
                    });
                    points.forEach(function(point, i) {
                      var py = y.call(rectbin, point, i) / dy;
                      var pj = trunc(py);
                      var px = x.call(rectbin, point, i) / dx;
                      var pi = trunc(px);
              
                      var id = pi + '-' + pj;
                      var bin = binsById[id];
                      bin.push(point);
                    });
                    return d3.values(binsById);
                  }
              
                  rectbin.x = function(_) {
                    if (!arguments.length) return x;
                    x = _;
                    return rectbin;
                  };
              
                  rectbin.y = function(_) {
                    if (!arguments.length) return y;
                    y = _;
                    return rectbin;
                  };
              
                  rectbin.dx = function(_) {
                    if (!arguments.length) return dx;
                    dx = _;
                    return rectbin;
                  };
              
                  rectbin.dy = function(_) {
                    if (!arguments.length) return dy;
                    dy = _;
                    return rectbin;
                  };
              
              
                  return rectbin;
                };
              
                var rectbinX = function(d) { return d[0]; },
                    rectbinY = function(d) { return d[1]; };
              
              })();
              
              function trunc(x) {
                return x < 0 ? Math.ceil(x) : Math.floor(x);
              }

            var size = 0.5
            var rectbinData = d3.rectbin()
                .dx(size)
                .dy(size)
                (inputForRectBinning)

            // Prepare a color palette
            var color = d3.scaleLinear()
                .domain([0, 350]) // Number of points in the bin?
                .range(["transparent", "#69a3b2"])

            // What is the height of a square in px?
            heightInPx = y(yLim[1] - size)

            // What is the width of a square in px?
            var widthInPx = x(xLim[0] + size)

            // Now we can add the squares
            svg.append("clipPath")
                .attr("id", "clip")
                .append("rect")
                .attr("width", width)
                .attr("height", height)
            svg.append("g")
                .attr("clip-path", "url(#clip)")
                .selectAll("myRect")
                .data(rectbinData)
                .enter().append("rect")
                .attr("x", function (d) { return x(d.x) })
                .attr("y", function (d) { return y(d.y) - heightInPx })
                .attr("width", widthInPx)
                .attr("height", heightInPx)
                .attr("fill", function (d) { return color(d.length); })
                .attr("stroke", "black")
                .attr("stroke-width", "0.4")
        });
    });

    return (<div id='data'></div>)
}

