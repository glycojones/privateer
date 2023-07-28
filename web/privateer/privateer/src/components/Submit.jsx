

export default function Submit({file, submitPressed}) { 
   
    return (
        <div class="flex flex-col items-center justify-center w-full m-12 h-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer bg-gray-50 dark:hover:bg-bray-800 dark:bg-gray-700 dark:border-gray-600">
            {file.name} uploaded successfully!
            <div class="py-6">
                <button class="" onClick={submitPressed}>Submit for analysis</button>
            </div>
        </div> 
    )
}