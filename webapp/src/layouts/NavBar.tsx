import { type Dispatch, type ReactElement, type SetStateAction } from 'react';
import { GENERAL_CITATION, GITHUB_REPO } from '../data/Constants.tsx';
import React from 'react';
export default function NavBar({
    setResetApp,
}: {
    setResetApp: Dispatch<SetStateAction<boolean>>;
}): ReactElement {
    return (
        <div className="flex flex-col w-full sm:flex-row sm:justify-between">
            <div className="text-center sm:text-left px-12 pt-12 sm:p-12 flex flex-col">
                <span className="font-primary text-xl text-secondary sm:text-3xl">
                    Validate your carbohydrates online with
                </span>
                <span className="font-body text-4xl sm:text-6xl my-2 sm:my-1">
                    <button
                        id="title"
                        title="Home"
                        className="bg-gray hover:scale-105 transition-all"
                        onClick={() => {
                            setResetApp(true);
                        }}
                    >
                        Privateer
                    </button>
                </span>
                <span className="font-primary text-l text-secondary sm:text-xl sm:my-4 my-2">
                    The Swiss Army knife for carbohydrate structure validation,
                    refinement and analysis
                </span>
            </div>

            <div className="flex">
                <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-6 flex items-center ">
                    <a href="/database">
                        <svg
                            className="w-6 h-6 text-gray-500 hover:scale-125 transition-all"
                            xmlns="http://www.w3.org/2000/svg"
                            height="1em"
                            viewBox="0 0 448 512"
                        >
                            <path d="M448 80v48c0 44.2-100.3 80-224 80S0 172.2 0 128V80C0 35.8 100.3 0 224 0S448 35.8 448 80zM393.2 214.7c20.8-7.4 39.9-16.9 54.8-28.6V288c0 44.2-100.3 80-224 80S0 332.2 0 288V186.1c14.9 11.8 34 21.2 54.8 28.6C99.7 230.7 159.5 240 224 240s124.3-9.3 169.2-25.3zM0 346.1c14.9 11.8 34 21.2 54.8 28.6C99.7 390.7 159.5 400 224 400s124.3-9.3 169.2-25.3c20.8-7.4 39.9-16.9 54.8-28.6V432c0 44.2-100.3 80-224 80S0 476.2 0 432V346.1z" />
                        </svg>
                    </a>
                </div>
                <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-6 flex items-center ">
                    <a href="/">
                        <svg
                            className="w-8 h-8 text-gray-500 hover:scale-125 transition-all"
                            aria-hidden="true"
                            xmlns="http://www.w3.org/2000/svg"
                            fill="none"
                            viewBox="0 0 20 16"
                        >
                            <path
                                stroke="currentColor"
                                strokeLinecap="round"
                                strokeLinejoin="round"
                                strokeWidth="2"
                                d="M13 13h3a3 3 0 0 0 0-6h-.025A5.56 5.56 0 0 0 16 6.5 5.5 5.5 0 0 0 5.207 5.021C5.137 5.017 5.071 5 5 5a4 4 0 0 0 0 8h2.167M10 15V6m0 0L8 8m2-2 2 2"
                            />
                        </svg>
                    </a>
                </div>
                <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-6 flex items-center ">
                    <a href={GENERAL_CITATION}>
                        <svg
                            className="w-6 h-6 text-gray-500 hover:scale-125 transition-all"
                            xmlns="http://www.w3.org/2000/svg"
                            height="1em"
                            viewBox="0 0 384 512"
                        >
                            <path d="M0 48V487.7C0 501.1 10.9 512 24.3 512c5 0 9.9-1.5 14-4.4L192 400 345.7 507.6c4.1 2.9 9 4.4 14 4.4c13.4 0 24.3-10.9 24.3-24.3V48c0-26.5-21.5-48-48-48H48C21.5 0 0 21.5 0 48z" />
                        </svg>
                    </a>
                </div>
                <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-12 flex items-center ">
                    <a href={GITHUB_REPO}>
                        <img
                            className="w-full hover:scale-125 transition-all hidden dark:block"
                            src="/github-mark.png"
                        />
                        <img
                            className="w-full hover:scale-125 transition-all block dark:hidden"
                            src="/github-mark-white.png"
                        />
                    </a>
                </div>
            </div>
        </div>
    );
}
