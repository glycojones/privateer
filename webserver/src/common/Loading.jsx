import ClipLoader from "react-spinners/ClipLoader";

export default function Loading({loadingText}) {

    return (
        <div
            className="flex flex-col items-center justify-center px-12 m-12 h-64 min-w-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer bg-gray-50 dark:hover:bg-bray-800 dark:bg-gray-700 dark:border-gray-600">
            <h3 className="my-6">{loadingText}</h3>
            <ClipLoader
                // color={}
                loading={true}
                size={100}
                aria-label="Loading Spinner"
                data-testid="loader"
            />

        </div>
    )
}