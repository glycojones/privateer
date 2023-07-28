import ClipLoader from "react-spinners/ClipLoader";

export default function Loading({file, submitPressed}) { 
   
    return (
        <div class="flex flex-col items-center justify-center w-full m-12 h-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer bg-gray-50 dark:hover:bg-bray-800 dark:bg-gray-700 dark:border-gray-600">
            <h3 className="my-6">Validating glycans...</h3>
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