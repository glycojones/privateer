export function Information({}) {
  return <div className="text-secondary bg-tertiary pb-12">
                <div className="flex flex-col">
                    <div className="text-center sm:text-left sm:text-center mb-4 px-12 pt-12 sm:p-12 flex flex-col">
                        <span className="font-body text-xl sm:text-3xl">How does it work?</span>
                        <span className="font-primary text-l mt-6 sm:text-xl mx-12">Privateer mainly uses <b>three</b> properties to assess the validity of a glycan chain. </span>
                    </div>

                    <div className="flex flex-col mx-24 gap-6">
                        <div className="flex mr-auto gap-6">
                            <img className="w-48 h-48 hidden sm:block" src="/Cremer-Pople.png"></img>
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Stereochemical Validation</span>
                                <span className="text-m font-body">
                                Cyclic carbohydrates usually have clear conformational preferences in terms of energy. Enzymes can force sugars into higher-energy conformations in order to achieve catalysis, and these conformations should indeed be kept in mind when modelling a sugar in the -1 site of an enzyme. However, most of the time ( i.e. in the rest of the cases) sugars stay in their original lowest-energy conformation, which for pyranoses is either a <sup>4</sup>C<sub>1</sub> or a <sup>1</sup>C<sub>4</sub> chair. 
                                Privateer checks each glycan for the correct stereochemistry and highlights when and where there are issues.
                                </span>
                            </div>
                        </div>

                        <div className="flex ml-auto gap-6">
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Real Space Correlation Coefficient Validation </span>
                                <span className="text-m font-body">The real space correlation coefficient metric in quantitative terms describes how well a modelled residue fits its associated experimental density, by calculating the difference between observed(experimental) structural factors and calculated structural factors associated with the fitted model. For monosaccharides that are modelled as components of glycans, RSCC values of above 0.80 are considered to signify a good model fit to its associated experimental density.</span>
                            </div>
                            <img className="w-48 h-48 hidden sm:block" src="/Cremer-Pople.png"></img>

                        </div>

                        <div className="flex mr-auto gap-6">
                            <img className="w-48 h-48 hidden sm:block" src="/Cremer-Pople.png"></img>
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Anomer Validation</span>
                                <span className="text-m font-body">Occasionally, the glycans modelled in a structure can exhibit anomeric configurations which have not been seen before, and are most likely the result of slightly incorrect glycan positioning. Privateer can highlight this to allow a user to correct the mistakes.</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>;
}
  