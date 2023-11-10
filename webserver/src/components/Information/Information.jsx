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
                    <img className="w-54 h-48 hidden sm:block" src="/glycomics.png"></img>
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Glycomics-Powered Validation</b></span>
                        <span className="text-m font-body">When modelling glycans,  it is possible to produce incorrect glycan compositions that do nor conform with general glycan biosynthesis knowledge. In Privateer, users can check their glycan structure and composition against glycomics databases. This allows the identification of inconsistent linkages, and following this, alternative compositions will be suggested to provide a more accurate, complete structure.</span>
                    </div>
                </div>

                <div className="flex ml-auto gap-6">
                    <div className="flex flex-col gap-6 justify-center">
                        <span className="text-xl font-body"><b>Torsion Angles</b></span>
                        <span className="text-m font-body">Torsion angle analysis allows the validation of the overall conformation of N-glycans, using a newly compiled set of glycosidic linkage torsional preferences harvested from a curated set of glycoprotein models.</span>
                    </div>
                    <img className="w-78 h-48 hidden sm:block" src="/torsion_figure.png"></img>
                </div>
            </div>
        </div>
    </div>;
}
  