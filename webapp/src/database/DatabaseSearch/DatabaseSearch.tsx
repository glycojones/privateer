import React, {
    type Dispatch,
    type SetStateAction,
    useEffect,
    useState,
} from 'react';
import DatabaseSearchTable from '../DatabaseSearchTable/DatabaseSearchTable.tsx';
import Loading from '../../shared/Loading/Loading.tsx';
import { sugarLinkageMap } from '../../data/Constants.tsx';
function ViewAllEntriesButton(props: {
    text: string;
    label: any;
    onClickMethod: any;
    secondLabel: string | undefined;
    onSecondClick: any;
}) {
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-12 pt-8 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                className="w-6 h-6 mb-4 text-gray-500 dark:text-gray-400"
                viewBox="0 0 576 512"
            >
                <path d="M64 32C64 14.3 49.7 0 32 0S0 14.3 0 32v96V384c0 35.3 28.7 64 64 64H256V384H64V160H256V96H64V32zM288 192c0 17.7 14.3 32 32 32H544c17.7 0 32-14.3 32-32V64c0-17.7-14.3-32-32-32H445.3c-8.5 0-16.6-3.4-22.6-9.4L409.4 9.4c-6-6-14.1-9.4-22.6-9.4H320c-17.7 0-32 14.3-32 32V192zm0 288c0 17.7 14.3 32 32 32H544c17.7 0 32-14.3 32-32V352c0-17.7-14.3-32-32-32H445.3c-8.5 0-16.6-3.4-22.6-9.4l-13.3-13.3c-6-6-14.1-9.4-22.6-9.4H320c-17.7 0-32 14.3-32 32V480z" />
            </svg>
            <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">{props.text}</span>
            </p>
            <div className="flex w-full">
                <button
                    type="button"
                    id="fetch"
                    className="bg-gray-50 border font-bold mb-8 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                    onClick={() => {
                        props.onClickMethod((value: any) => !value);
                    }}
                >
                    {props.label}
                </button>
                {props.secondLabel !== undefined ? (
                    <button
                        type="button"
                        id="fetch"
                        className="bg-gray-50 border ml-2 font-bold mb-8 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                        onClick={() => {
                            props.onSecondClick((value: any) => !value);
                        }}
                    >
                        {props.secondLabel}
                    </button>
                ) : (
                    <></>
                )}
            </div>
        </div>
    );
}

function TypeFilterBox(props: { selected: string; onClickMethod: any }) {
    const sugars = ['N-glycans', 'O-glycans', 'S-glycans', 'C-glycans'];
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-2 pt-6 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                viewBox="0 0 448 512"
                className="w-6 h-6 mb-2 text-gray-500 dark:text-gray-400"
            >
                <path d="M192 64v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V64c0-17.7-14.3-32-32-32H224c-17.7 0-32 14.3-32 32zM82.7 207c-15.3 8.8-20.5 28.4-11.7 43.7l32 55.4c8.8 15.3 28.4 20.5 43.7 11.7l55.4-32c15.3-8.8 20.5-28.4 11.7-43.7l-32-55.4c-8.8-15.3-28.4-20.5-43.7-11.7L82.7 207zM288 192c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V224c0-17.7-14.3-32-32-32H288zm64 160c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H352zM160 384v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H192c-17.7 0-32 14.3-32 32zM32 352c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32H96c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H32z" />
            </svg>

            <p className="mb-1 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">Filter Sugars</span>
            </p>

            <div className="w-full flex flex-wrap my-auto overflow-x-hidden justify-center align-items-center ">
                {sugars.map((item, index) => {
                    return (
                        <button
                            type="button"
                            key={index}
                            id="fetch"
                            className={
                                props.selected === item
                                    ? 'bg-gray-50 mx-2 my-2 border w-24 font-bold border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                                    : 'bg-gray-50 mx-2 my-2 border w-24 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                            }
                            onClick={() => {
                                props.onClickMethod(item);
                            }}
                        >
                            {item}
                        </button>
                    );
                })}
            </div>
        </div>
    );
}

function LinkageFilterBox(props: {
    selected: string;
    sugars: string[];
    onClickMethod: Dispatch<SetStateAction<unknown>>;
}) {
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-2 pt-6 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                viewBox="0 0 640 512"
                className="w-6 h-6 mb-2 text-gray-500 dark:text-gray-400"
            >
                <path d="M579.8 267.7c56.5-56.5 56.5-148 0-204.5c-50-50-128.8-56.5-186.3-15.4l-1.6 1.1c-14.4 10.3-17.7 30.3-7.4 44.6s30.3 17.7 44.6 7.4l1.6-1.1c32.1-22.9 76-19.3 103.8 8.6c31.5 31.5 31.5 82.5 0 114L422.3 334.8c-31.5 31.5-82.5 31.5-114 0c-27.9-27.9-31.5-71.8-8.6-103.8l1.1-1.6c10.3-14.4 6.9-34.4-7.4-44.6s-34.4-6.9-44.6 7.4l-1.1 1.6C206.5 251.2 213 330 263 380c56.5 56.5 148 56.5 204.5 0L579.8 267.7zM60.2 244.3c-56.5 56.5-56.5 148 0 204.5c50 50 128.8 56.5 186.3 15.4l1.6-1.1c14.4-10.3 17.7-30.3 7.4-44.6s-30.3-17.7-44.6-7.4l-1.6 1.1c-32.1 22.9-76 19.3-103.8-8.6C74 372 74 321 105.5 289.5L217.7 177.2c31.5-31.5 82.5-31.5 114 0c27.9 27.9 31.5 71.8 8.6 103.9l-1.1 1.6c-10.3 14.4-6.9 34.4 7.4 44.6s34.4 6.9 44.6-7.4l1.1-1.6C433.5 260.8 427 182 377 132c-56.5-56.5-148-56.5-204.5 0L60.2 244.3z" />
            </svg>

            <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">Filter Linkages</span>
            </p>

            <div className="w-full flex flex-wrap overflow-scroll my-auto overflow-x-hidden justify-center align-items-center ">
                {props.sugars.map((item, index) => {
                    return (
                        <button
                            key={index}
                            type="button"
                            id="fetch"
                            className={
                                props.selected === item
                                    ? 'bg-gray-50 mx-2 my-2 border w-32 font-bold border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                                    : 'bg-gray-50 mx-2 my-2 border w-32 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                            }
                            onClick={() => {
                                props.onClickMethod(item);
                            }}
                        >
                            {item}
                        </button>
                    );
                })}
            </div>
        </div>
    );
}

function FilterZone(props: { setSearchBegin: any }) {
    const [search, setSearch] = useState<boolean>(false);
    const [linkage, setLinkage] = useState<string>('Any');
    const [type, setType] = useState<string>('N-glycans');
    const [text, setText] = useState<string>('');

    useEffect(() => {
        setLinkage(sugarLinkageMap[type][0]);
    }, [type]);

    useEffect(() => {
        let textString = 'Find';

        if (type === 'Any') {
            textString += ' any carbohydrate';
        } else {
            textString += ' ' + type;
        }
        textString += ' with ';
        if (linkage === 'Any') {
            textString += ' any linkage type';
        } else {
            textString += ' ' + linkage + ' linkages';
        }

        setText(textString);
    }, [linkage, type]);

    useEffect(() => {
        if (!search) {
            return;
        }
        let formattedType = '';
        let url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/linkages/';
        if (type === 'N-glycans') {
            url += 'n-glycan/';
            formattedType = 'n-glycan';
        }
        if (type === 'O-glycans') {
            url += 'o-glycan/';
            formattedType = 'o-glycan';
        }
        if (type === 'S-glycans') {
            url += 's-glycan/';
            formattedType = 's-glycan';
        }
        if (type === 'C-glycans') {
            url += 'c-glycan/';
            formattedType = 'c-glycan';
        }
        if (type === 'Ligands') {
            url += 'ligand/';
            formattedType = 'ligand';
        }
        if (linkage === 'Any') {
            url += 'any';
        } else {
            url += linkage.replace(',', '%2C');
        }
        url += '.json';

        void fetch(url)
            .then(async (response) => await response.json())
            .then((json: Record<any, any>) => {
                let formattedData: Record<string, any>;
                if (linkage === 'Any') {
                    formattedData = Object.keys(json).flatMap((item) => {
                        return json[item].map((element) => {
                            return {
                                pdb: element.pdb,
                                count: element.count,
                                resolution: element.resolution,
                                linkage: item,
                                type: formattedType,
                                link:
                                    'https://privateer.york.ac.uk/database?pdb=' +
                                    element.pdb,
                            };
                        });
                    });
                } else {
                    formattedData = json.map((item) => {
                        return {
                            pdb: item.pdb,
                            count: item.count,
                            resolution: item.resolution,
                            linkage,
                            type: formattedType,
                            link:
                                'https://privateer.york.ac.uk/database?pdb=' +
                                item.pdb,
                        };
                    });
                }
                setData(formattedData);
            })
            .catch(() => {});
    }, [search]);

    const [data, setData] = useState<Record<string, string | number> | null>(
        null
    );

    return !search ? (
        <>
            {ViewAllEntriesButton({
                label: 'Back',
                text,
                onClickMethod: props.setSearchBegin,
                secondLabel: 'Search',
                onSecondClick: setSearch,
            })}
            {TypeFilterBox({ onClickMethod: setType, selected: type })}
            {LinkageFilterBox({
                onClickMethod: setLinkage,
                sugars: sugarLinkageMap[type],
                selected: linkage,
            })}
        </>
    ) : data === null ? (
        <Loading loadingText={'Crunching Numbers'} />
    ) : (
        <div className="flex flex-col w-full text-center">
            <span className="font-semibold">Query: {text}</span>
            <button
                className="align-self-start w-32 mb-2"
                onClick={() => {
                    setSearch(false);
                }}
            >
                &#8592; Back to Search
            </button>
            <DatabaseSearchTable data={data} />
        </div>
    );
}

export default function DatabaseSearch(props: {
    setSearchBegin: any;
    searchBegin: boolean;
}) {
    return !props.searchBegin
        ? ViewAllEntriesButton({
              label: 'View',
              text: (
                  <>
                      View all entries in the <i>Privateer</i> Database
                  </>
              ),
              onClickMethod: props.setSearchBegin,
          })
        : FilterZone({ setSearchBegin: props.setSearchBegin });
}
