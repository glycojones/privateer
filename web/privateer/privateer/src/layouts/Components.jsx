import React from "react";
export function Components({}) {
  return <div className="text-secondary bg-tertiary pb-12">
                <div className="flex flex-col">
                    <div className="text-center sm:text-left sm:text-center mb-4 px-12 pt-12 sm:p-12 flex flex-col">
                        <span className="font-body text-xl sm:text-3xl">How does it work?</span>
                        <span className="font-primary text-l mt-6 sm:text-xl mx-12">Privateer mainly uses <b>three</b> properties to assess the validity of a glycan chain. </span>
                    </div>

                    <div className="flex flex-col mx-24 gap-6">
                        <div className="flex mr-auto gap-6">
                            <img className="w-48 h-48 hidden sm:block" src="../../public/Cremer-Pople.png"></img>
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Stereochemistry</span>
                                <span className="text-m font-body">Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.</span>
                            </div>
                        </div>

                        <div className="flex ml-auto gap-6">
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Real Space Correlation Coefficient</span>
                                <span className="text-m font-body">Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.</span>
                            </div>
                            <img className="w-48 h-48 hidden sm:block" src="../../public/Cremer-Pople.png"></img>

                        </div>

                        <div className="flex mr-auto gap-6">
                            <img className="w-48 h-48 hidden sm:block" src="../../public/Cremer-Pople.png"></img>
                            <div className="flex flex-col gap-6 justify-center">
                                <span className="text-xl font-body">Anomeric</span>
                                <span className="text-m font-body">Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industry's standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>;
}
  