import { useEffect, useState } from 'react';

export default function DatabaseFetch({ PDBCode, setPDBCode, submitPressed }) {
    const [pdb, setPDB] = useState('');

    useEffect(() => {
        if (pdb.length != 4) {
            return;
        }
        if (!pdb.match(/^[a-z0-9]+$/i)) {
            return;
        }

        setPDBCode(pdb);
        submitPressed(true);
    }, [pdb]);

    return (
        <>
            {PDBCode != true ? (
                <div className="flex items-center justify-center m-12 w-64 ">
                    <label className="flex flex-col items-center justify-center w-full p-12 h-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer border-gray-600">
                        <div className="flex flex-col items-center justify-center pt-5 pb-6 text-center">
                            <svg
                                className="w-6 h-6 mb-4 text-gray-500 dark:text-gray-400"
                                xmlns="http://www.w3.org/2000/svg"
                                height="1em"
                                viewBox="0 0 448 512"
                            >
                                <path d="M448 80v48c0 44.2-100.3 80-224 80S0 172.2 0 128V80C0 35.8 100.3 0 224 0S448 35.8 448 80zM393.2 214.7c20.8-7.4 39.9-16.9 54.8-28.6V288c0 44.2-100.3 80-224 80S0 332.2 0 288V186.1c14.9 11.8 34 21.2 54.8 28.6C99.7 230.7 159.5 240 224 240s124.3-9.3 169.2-25.3zM0 346.1c14.9 11.8 34 21.2 54.8 28.6C99.7 390.7 159.5 400 224 400s124.3-9.3 169.2-25.3c20.8-7.4 39.9-16.9 54.8-28.6V432c0 44.2-100.3 80-224 80S0 476.2 0 432V346.1z" />
                            </svg>
                            <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                                <span className="font-semibold">
                                    Fetch from the <i>Privateer</i> Database
                                </span>
                            </p>
                            <input
                                type="text"
                                id="code"
                                className="bg-gray-50 border border-gray-300 text-center text-gray-900 text-sm rounded-lg focus:border-3 block w-full p-2.5 my-2 "
                                placeholder="5FJI"
                                required
                                onKeyDown={(e) => {
                                    if (e.key == 'Enter') {
                                        let element =
                                            document.getElementById('code');
                                        setPDB(element.value);
                                    }
                                }}
                            />
                            <button
                                type="button"
                                id="fetch"
                                className="bg-gray-50 border font-bold  border-gray-300 text-gray-900 text-sm rounded-lg  focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                                onClick={() => {
                                    let element =
                                        document.getElementById('code');
                                    setPDB(element.value);
                                }}
                            >
                                Fetch
                            </button>
                        </div>
                    </label>
                </div>
            ) : (
                <></>
            )}
        </>
    );
}
