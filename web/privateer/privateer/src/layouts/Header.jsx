import SNFG from '../components/SNFG'
import Upload from "../components/Upload";
import Submit from "../common/Submit";
import Loading from "../common/Loading";
import { GITHUB_REPO } from "../data/Constants"

export function Header({
  setResetApp,
  file,
  setFile,
  submit,
  setSubmit,
  tableData,
  loadingText,
  fileContent
}) {
  return <div className="bg-gray text-primary">
                <div className="flex flex-col sm:flex-row sm:justify-between">
                    <div className="text-center sm:text-left px-12 pt-12 sm:p-12 flex flex-col">
                        <span className="font-primary text-xl text-secondary sm:text-3xl">Validate your N-glycans online with</span>
                        <span className="font-body text-4xl sm:text-6xl my-2 sm:my-1"><button id="title" title="Home" onClick={() => {
            setResetApp(true);
          }}>Privateer</button></span>
                        <span className="font-primary text-l text-secondary sm:text-xl sm:my-4 my-2">The Swiss Army knife for carbohydrate structure validation, refinement and analysis</span>
                    </div>
                    <div className="h-12 w-12 mx-auto my-4 sm:w-12 sm:mt-12 sm:mr-12 flex items-center ">
                        <a href={GITHUB_REPO}>
                            <img className="w-full hover:scale-125 transition-all hidden dark:block" src="../../public/github-mark.png" />
                            <img className="w-full hover:scale-125 transition-all block dark:hidden" src="../../public/github-mark-white.png" />
                        </a>
                    </div>
                </div>

                <div className="flex justify-center mb-6">
                {file == null ? <Upload setFile={setFile} /> : submit == null ? <Submit file={file} submitPressed={setSubmit} /> : tableData == null ? <Loading loadingText={loadingText} /> : <SNFG tableData={tableData} fileName={file.name} pdbString={fileContent}></SNFG>}
                </div>

            </div>;
}
  