import React, { useEffect, useRef, useState } from 'react';
import { MoorhenContainer } from 'moorhen';
// import LeftMouseClick from '../../assets/LeftClick';

export function GlycanDetailMoorhenView (props: {
  key: string
  onChange: (e) => Promise<void>
  onSymmetryChange: (e) => any
  moorhenProps: any
  moorhenDimensions: () => [number, number]
  mapContour: number
}) {
  const size = Math.min(800, 0.9 * window.innerWidth);
  const [dimensions, setDimensions] = useState<Record<string, number>>(
    size >= 800
      ? { width: size, height: size * (3 / 4) }
      : { width: size, height: size }
  );

  const dimensionRef = useRef();

  // @ts-expect-error
  dimensionRef.current = dimensions;
  useEffect(() => {
    function handleResize () {
      const size = Math.min(800, 0.9 * window.innerWidth);
      if (size >= 800) {
        setDimensions({ width: size, height: size * (3 / 4) });
      } else {
        setDimensions({ width: size, height: size });
      }
    }

    window.addEventListener('resize', handleResize);
  }, []);

  const moorhenDimensionCallback = (): [number, number] => {
    // @ts-expect-error
    return [dimensionRef.current.width, dimensionRef.current.height];
  };

  return (
    <div key={props.key} className="px-8 flex flex-col items-center">
      <h3 className="text-left text-xl w-full font-bold mt-2">
        Visualise with <i>Moorhen</i>{' '}
        <a href="https://moorhen.org" title="Go to Moorhen.org">
          <img className="inline h-8" src="./moorhen_logo.png"></img>
        </a>
      </h3>

      <div className="mx-auto mt-4 items-center" style={{ transform: 'translateY(-60px)' }}>
        <MoorhenContainer
          {...props.moorhenProps}
          setMoorhenDimensions={moorhenDimensionCallback}
          viewOnly={true}
        />
      </div>

      <div className="flex gap-8" style={{ transform: 'translateY(-40px)' }}>
        <div className="block">
          <label
            htmlFor="contour-range-text"
            className="block mt-2 text-sm font-medium text-gray-909"
          >
            Map Contour - {props.mapContour}{' '}
          </label>
          <input
            id="contour-range"
            type="range"
            min="0"
            max="2"
            step="0.05"
            defaultValue="0.2"
            className="w-36 mt-2 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
            onChange={props.onChange}
            onMouseDown={(e) => {
              e.stopPropagation();
            }}
            onTouchStart={(e) => {
              e.stopPropagation();
            }}
          />
        </div>

        <div className="mt-4">
          <label className="relative inline-flex items-center cursor-pointer">
            <input
              type="checkbox"
              value=""
              className="sr-only peer"
              onMouseDown={(e) => {
                e.stopPropagation();
              }}
              onTouchStart={(e) => {
                e.stopPropagation();
              }}
            />
            <div
              onMouseDown={(e) => {
                e.stopPropagation();
              }}
              onTouchStart={(e) => {
                e.stopPropagation();
              }}
              onClick={props.onSymmetryChange}
              className="w-11 h-6 bg-tertiary peer-focus:outline-none peer-focus:ring-4 peer-focus:ring-blue-300 dark:peer-focus:ring-blue-800 rounded-full peer dark:bg-gray-700 peer-checked:after:translate-x-full rtl:peer-checked:after:-translate-x-full peer-checked:after:border-white after:content-[''] after:absolute after:top-[2px] after:start-[2px] after:bg-white after:border-gray-300 after:border after:rounded-full after:h-5 after:w-5 after:transition-all dark:border-gray-600 peer-checked:bg-blue-600"
            ></div>
            <span className="ms-3 text-sm font-medium text-gray-900 dark:text-gray-300">
              Toggle Symmetry
            </span>
          </label>
        </div>
      </div>
      {/* <div className="flex mt-2 items-center justify-center"><LeftMouseClick /> Rotate Molecule</div> */}
    </div>
  );
}
