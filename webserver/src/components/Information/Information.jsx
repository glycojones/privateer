import {TORSIONS_CITATION, GLYCOMICS_CITATION} from "../../data/Constants"

export function Information({}) {
    return <div className="text-secondary bg-tertiary pb-12">
        <div className="flex flex-col">
            <div className="text-center sm:text-left sm:text-center mb-4 px-12 pt-12 sm:p-12 flex flex-col">
                <span className="font-body text-xl sm:text-3xl"><b>How does it work?</b></span>
                <span className="font-body text-m mt-4">
                Carbohydrates, including O- and N-glycans attached to protein and lipid structures, are increasingly important in cellular biology. Crystallographic refinement of sugars is, however, very poorly performed, with thousands of incorrect structures polluting the Protein Data Bank. <i>Privateer</i> is a software that aims to detect and prevent conformational, regiochemical and stereochemical anomalies in cyclic carbohydrate structures. Multiple features can be used to assess the validity of a glycan chain. </span>
            </div>

            <div className="flex flex-col mx-24 gap-6">
                <div className="flex mr-auto gap-6">
                    <img className="w-74 h-44 hidden sm:block" src="/Cremer-Pople.png"></img>
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Stereochemical Validation</b></span>
                        <span className="text-m font-body">
                        Cyclic carbohydrates usually have clear conformational preferences in terms of energy. Whilst enzymes can force sugars into higher-energy conformations in order to achieve catalysis, sugars tend to stay in their original lowest-energy conformation, which for pyranoses is either a <sup>4</sup>C<sub>1</sub> or a <sup>1</sup>C<sub>4</sub> chair. Privateer checks each glycan for the correct stereochemistry and highlights any outliers.
                                </span>
                    </div>
                </div>

                <div className="flex ml-auto gap-6">
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Real Space Correlation Coefficient Validation</b></span>
                        <span className="text-m font-body">The real space correlation coefficient metric quantitatively describes how well a modelled residue fits its associated experimental density. It is defined by calculating the difference between observed (experimental) structural factors and calculated structural factors associated with the fitted model. For monosaccharides that are modelled as components of glycans, RSCC values of above 0.80 are considered to signify a good model fit to its associated experimental density. A real space correlation coefficient against omit mFo-DFc electron density is reported as a quality of fit indicator.</span>
                    </div>
                    <img className="w-76 h-48 hidden sm:block" src="/refmac_figure.png"></img>

                </div>

                <div className="flex mr-auto gap-6">
                    <img className="w-74 h-64 hidden sm:block" src="/glycomics.png"></img>
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Glycomics-Powered Validation</b>
                        <a href={GLYCOMICS_CITATION}>
                                <svg className="inline pl-2 h-4 transition-all hover:h-5 fill-secondary" xmlns="http://www.w3.org/2000/svg" height="1em" viewBox="0 0 640 512"><path d="M579.8 267.7c56.5-56.5 56.5-148 0-204.5c-50-50-128.8-56.5-186.3-15.4l-1.6 1.1c-14.4 10.3-17.7 30.3-7.4 44.6s30.3 17.7 44.6 7.4l1.6-1.1c32.1-22.9 76-19.3 103.8 8.6c31.5 31.5 31.5 82.5 0 114L422.3 334.8c-31.5 31.5-82.5 31.5-114 0c-27.9-27.9-31.5-71.8-8.6-103.8l1.1-1.6c10.3-14.4 6.9-34.4-7.4-44.6s-34.4-6.9-44.6 7.4l-1.1 1.6C206.5 251.2 213 330 263 380c56.5 56.5 148 56.5 204.5 0L579.8 267.7zM60.2 244.3c-56.5 56.5-56.5 148 0 204.5c50 50 128.8 56.5 186.3 15.4l1.6-1.1c14.4-10.3 17.7-30.3 7.4-44.6s-30.3-17.7-44.6-7.4l-1.6 1.1c-32.1 22.9-76 19.3-103.8-8.6C74 372 74 321 105.5 289.5L217.7 177.2c31.5-31.5 82.5-31.5 114 0c27.9 27.9 31.5 71.8 8.6 103.9l-1.1 1.6c-10.3 14.4-6.9 34.4 7.4 44.6s34.4 6.9 44.6-7.4l1.1-1.6C433.5 260.8 427 182 377 132c-56.5-56.5-148-56.5-204.5 0L60.2 244.3z"/></svg>
                            </a>
                        </span>
                        <span className="text-m font-body">When modelling glycans,  it is possible to produce incorrect glycan compositions that do nor conform with general glycan biosynthesis knowledge. In Privateer, users can check their glycan structure and composition against glycomics databases. This allows the identification of inconsistent linkages, and following this, alternative compositions will be suggested to provide a more accurate, complete structure.</span>
                    </div>
                </div>

                <div className="flex ml-auto gap-6">
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Torsion Angles</b> 
                            <a href={TORSIONS_CITATION}>
                                <svg className="inline pl-2 h-4 transition-all hover:h-5 fill-secondary" xmlns="http://www.w3.org/2000/svg" height="1em" viewBox="0 0 640 512"><path d="M579.8 267.7c56.5-56.5 56.5-148 0-204.5c-50-50-128.8-56.5-186.3-15.4l-1.6 1.1c-14.4 10.3-17.7 30.3-7.4 44.6s30.3 17.7 44.6 7.4l1.6-1.1c32.1-22.9 76-19.3 103.8 8.6c31.5 31.5 31.5 82.5 0 114L422.3 334.8c-31.5 31.5-82.5 31.5-114 0c-27.9-27.9-31.5-71.8-8.6-103.8l1.1-1.6c10.3-14.4 6.9-34.4-7.4-44.6s-34.4-6.9-44.6 7.4l-1.1 1.6C206.5 251.2 213 330 263 380c56.5 56.5 148 56.5 204.5 0L579.8 267.7zM60.2 244.3c-56.5 56.5-56.5 148 0 204.5c50 50 128.8 56.5 186.3 15.4l1.6-1.1c14.4-10.3 17.7-30.3 7.4-44.6s-30.3-17.7-44.6-7.4l-1.6 1.1c-32.1 22.9-76 19.3-103.8-8.6C74 372 74 321 105.5 289.5L217.7 177.2c31.5-31.5 82.5-31.5 114 0c27.9 27.9 31.5 71.8 8.6 103.9l-1.1 1.6c-10.3 14.4-6.9 34.4 7.4 44.6s34.4 6.9 44.6-7.4l1.1-1.6C433.5 260.8 427 182 377 132c-56.5-56.5-148-56.5-204.5 0L60.2 244.3z"/></svg>
                            </a>
                        </span>

                        <span className="text-m font-body">                            
                            Torsion angle analysis allows the validation of the overall conformation of N-glycans. Using a compiled set of glycosidic linkage torsional preferences harvested from a curated set of glycoprotein models, a Z score for each proposed glycan linkage is calculated and displayed within the Symbol Nomenclature For Glycan (SNFG) image.</span>
                    </div>
                    <img className="w-78 h-48 hidden sm:block" src="/torsion_figure.png"></img>
                </div>
            </div>
        </div>
    </div>;
}
  