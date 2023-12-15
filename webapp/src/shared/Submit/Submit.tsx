import Molecule from '../../assets/Molecule';
import Reflections from '../../assets/Reflections.tsx';
import React, {
    type Dispatch,
    type ReactElement,
    type SetStateAction,
} from 'react';
interface FileLineProps {
    icon: any;
    name: string;
}
function FileLine(props: FileLineProps): ReactElement {
    return (
        <div className="flex flex-row w-64 justify-center px-12 py-2 rounded-lg min-w-0">
            <div className="pr-2">{props.icon}</div>
            <p className="text-ellipsis overflow-hidden">{props.name}</p>
        </div>
    );
}

interface SubmitProps {
    coordinateFile: File;
    reflectionFile: File;
    submitPressed: boolean;
    setResetApp: Dispatch<SetStateAction<boolean>>;
    allowSubmit: boolean;
}
export default function Submit({
    coordinateFile,
    reflectionFile,
    submitPressed,
    setResetApp,
    allowSubmit,
}: SubmitProps): ReactElement {
    function getFileList(): ReactElement[] {
        const array: ReactElement[] = [];
        if (coordinateFile !== null) {
            array.push(
                <FileLine
                    icon={<Molecule />}
                    name={coordinateFile.name}
                    key={'c'}
                />
            );
        }
        if (reflectionFile !== null) {
            array.push(
                <FileLine
                    icon={<Reflections />}
                    name={reflectionFile.name}
                    key={'r'}
                />
            );
        }
        return array;
    }

    return (
        <div
            className="flex flex-col m-12 px-12 pt-8 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0"
        >
            <p>File uploaded:</p>
            {getFileList()}

            <div className="flex space-x-4 py-6">
                <button
                    className="bg-gray hover:bg-hover border-gray-800 border-2 text-primary opacity-60 font-bold py-2 px-4 rounded"
                    onClick={() => {
                        setResetApp(true);
                    }}
                >
                    Cancel
                </button>

                {allowSubmit ? (
                    <button
                        className="bg-gray hover:bg-hover border-gray-300 border-2 text-primary font-bold py-2 px-4 rounded"
                        onClick={submitPressed}
                    >
                        Submit
                    </button>
                ) : (
                    <button
                        className="bg-gray border-gray-800 border-2 text-primary opacity-40 font-bold py-2 px-4 rounded"
                        onClick={submitPressed}
                        disabled
                    >
                        Submit
                    </button>
                )}
            </div>
        </div>
    );
}
