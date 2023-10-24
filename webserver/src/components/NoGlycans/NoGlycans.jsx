

export default function NoGlycans({
    setResetApp, text
}) {
    return (
        <div
            className="flex flex-col m-12 px-12 pt-8 w-64 
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg 
            bg-gray-50 flex-grow-0">

            <p>{text}</p>

            <div className="flex space-x-4 py-6">
                <button
                    className="bg-gray hover:bg-hover border-gray-800 border-2 text-primary opacity-60 font-bold py-2 px-4 rounded"
                    onClick={() => { 

                        setResetApp(true)
                         }}>Retry</button>
            </div>

        </div>
    )
}