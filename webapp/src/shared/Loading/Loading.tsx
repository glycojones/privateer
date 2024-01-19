import ClipLoader from 'react-spinners/ClipLoader';
import React from 'react';
export default function Loading(props: { loadingText: string }) {
    return (
        <div className="flex flex-col items-center text-center justify-center mx-auto px-12 m-12 h-64 w-64 border-2 border-gray-300 rounded-lg cursor-pointer bg-gray-50 dark:hover:bg-bray-800 dark:bg-gray-700 dark:border-gray-600">
            <h3 className="my-6">{props.loadingText}</h3>
            <ClipLoader
                loading={true}
                size={100}
                aria-label="Loading Spinner"
                data-testid="loader"
            />
        </div>
    );
}
