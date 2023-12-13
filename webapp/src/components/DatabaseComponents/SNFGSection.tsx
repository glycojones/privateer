import React, { useEffect, useState } from "react";

export default function SNFGSection(props) {
  const [showSection, setShowSection] = useState<boolean>(false);

  useEffect(() => {
    console.log(props);
  }, []);

  return (
    <div className="flex flex-col">
      <span className="text-xl text-left sm:text-center">
        Chain: {props.id}
        <button
          className="ml-2"
          onClick={() => {
            setShowSection((x) => !x);
          }}
        >
          {showSection ? (
            <svg
              xmlns="http://www.w3.org/2000/svg"
              height="1em"
              viewBox="0 0 320 512"
            >
              <path d="M182.6 137.4c-12.5-12.5-32.8-12.5-45.3 0l-128 128c-9.2 9.2-11.9 22.9-6.9 34.9s16.6 19.8 29.6 19.8H288c12.9 0 24.6-7.8 29.6-19.8s2.2-25.7-6.9-34.9l-128-128z" />
            </svg>
          ) : (
            <svg
              xmlns="http://www.w3.org/2000/svg"
              height="1em"
              viewBox="0 0 320 512"
            >
              <path d="M137.4 374.6c12.5 12.5 32.8 12.5 45.3 0l128-128c9.2-9.2 11.9-22.9 6.9-34.9s-16.6-19.8-29.6-19.8L32 192c-12.9 0-24.6 7.8-29.6 19.8s-2.2 25.7 6.9 34.9l128 128z" />
            </svg>
          )}
        </button>
      </span>

      {showSection ? (
        props.item.map((array_item, array_index) => {
          return <SNFGItem item={array_item} index={array_index} />;
        })
      ) : (
        <></>
      )}
    </div>
  );
}

export function SNFGItem(props) {
  const [showWURCS, setShowWURCS] = useState<boolean>(false);

  return (
    <div className="flex flex-col items-center justify-center collapsible mt-2">
      <div className="" key={`wurcs-${props.index}`}>
        WURCS:
        <button
          className="w-4 mr-4"
          onClick={() => {
            setShowWURCS((x) => !x);
          }}
        >
          {!showWURCS ? (
            <svg
              className="w-4 h-4 text-gray-500 "
              xmlns="http://www.w3.org/2000/svg"
              height="1em"
              viewBox="0 0 256 512"
            >
              <path d="M246.6 278.6c12.5-12.5 12.5-32.8 0-45.3l-128-128c-9.2-9.2-22.9-11.9-34.9-6.9s-19.8 16.6-19.8 29.6l0 256c0 12.9 7.8 24.6 19.8 29.6s25.7 2.2 34.9-6.9l128-128z" />
            </svg>
          ) : (
            <svg
              className="w-4 h-4 text-gray-500 "
              xmlns="http://www.w3.org/2000/svg"
              height="1em"
              viewBox="0 0 256 512"
            >
              <path d="M9.4 278.6c-12.5-12.5-12.5-32.8 0-45.3l128-128c9.2-9.2 22.9-11.9 34.9-6.9s19.8 16.6 19.8 29.6l0 256c0 12.9-7.8 24.6-19.8 29.6s-25.7 2.2-34.9-6.9l-128-128z" />
            </svg>
          )}
        </button>
        {showWURCS ? props.item.WURCS : <></>}
      </div>
      <div
        className="mt-4 py-4"
        id="svgContainer"
        key={`snfg-${props.index}`}
        dangerouslySetInnerHTML={{
          __html: props.item.SNFG,
        }}
      />
    </div>
  );
}
