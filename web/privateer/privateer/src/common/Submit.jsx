

export default function Submit({file, submitPressed}) { 
   
    return (
        <div className="flex flex-col m-12 p-12 items-center text-center justify-between h-64 border-2 border-dashed transition-all border-gray-300 rounded-lg cursor-pointer bg-gray-50 dark:hover:bg-bray-800 dark:bg-gray-700 dark:border-gray-600">
            {file.name} uploaded successfully!
            <div className="py-6">
                <button className="bg-gray hover:bg-hover border-gray-300 border-2 text-white font-bold py-2 px-4 rounded" onClick={submitPressed}>Submit for analysis</button>
            </div>
        </div> 
    )
}