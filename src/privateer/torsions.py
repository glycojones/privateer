import os
from privateer import privateer_core as pvtcore
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
from datetime import datetime
import argparse


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)


class PrivateerTorsionResultsOutputParser:
    def __init__(self, privateer_output):
        self.__type = "uninitialized"
        self.__privateer_output = privateer_output
        self.__type = self.__validateInputJSON()
        self.__parsedTorsionsOfStructureOutputForPlotting = None
        if self.__type == "structure":
            self.__parsedTorsionsOfStructureOutputForPlotting = (
                self.__parseOutputIfTorsionsStructureProvided()
            )

    def getTypeOfOutput(self):
        return self.__type

    def getSortedOutputOfTorsionsForEntireStructure(self):
        return self.__parsedTorsionsOfStructureOutputForPlotting

    def __validateInputJSON(self):
        expected_keys_structure = ["glycanIndex", "WURCS", "all_torsions_in_structure"]
        expected_keys_glycan = [
            "first_residue",
            "second_residue",
            "root_descr",
            "detected_torsions",
            "database_phi",
            "database_psi",
        ]
        if len(self.__privateer_output):
            item = self.__privateer_output[0]
            rootKeys = item.keys()
            if sorted(rootKeys) == sorted(expected_keys_glycan):
                self.__type = "glycan"
            elif sorted(rootKeys) == sorted(expected_keys_structure):
                glycan = item["all_torsions_in_structure"]
                glycanItem = glycan[0]
                glycanKeys = glycanItem.keys()
                if sorted(glycanKeys) == sorted(expected_keys_glycan):
                    self.__type = "structure"
                else:
                    raise AttributeError(
                        f"Root structure of imported output resembles that of a glycan, but keys of individual glycans do not match.\nReceived keys {sorted(glycanKeys)}\nExpected keys {expected_keys_glycan}"
                    )
            else:
                raise AttributeError(
                    f"Imported Privateer output does not match expected output from Privateer."
                )

        else:
            raise AttributeError("Imported torsions input from Privateer is empty.")

        return self.__type

    def __parseOutputIfTorsionsStructureProvided(self):
        def collectUniqueResiduePairs(input):
            pairs = []
            for rootItem in input:
                all_torsions = rootItem["all_torsions_in_structure"]
                for item in all_torsions:
                    pair = [
                        item["first_residue"],
                        item["second_residue"],
                    ]
                    if pair not in pairs:
                        pairs.append(pair)
                    else:
                        continue

            return pairs

        unique_pairs = collectUniqueResiduePairs(self.__privateer_output)
        output = []
        for unique_pair in unique_pairs:
            first_unique_residue = unique_pair[0]
            second_unique_residue = unique_pair[1]
            alreadyAdded = False
            unique_pair_dict = {
                "first": None,
                "second": None,
                "database_phi": None,
                "database_psi": None,
                "detected_torsions_in_structures": [],
            }
            for rootItem in self.__privateer_output:
                all_torsions = rootItem["all_torsions_in_structure"]
                for item in all_torsions:
                    if (
                        item["first_residue"] == first_unique_residue
                        and item["second_residue"] == second_unique_residue
                    ):
                        detected_torsions = {
                            "glycanIndex": rootItem["glycanIndex"],
                            "WURCS": rootItem["WURCS"],
                            "glycan_torsions": item["detected_torsions"],
                        }
                        if alreadyAdded:
                            unique_pair_dict["detected_torsions_in_structures"].append(
                                detected_torsions
                            )
                        else:
                            unique_pair_dict["first"] = item["first_residue"]
                            unique_pair_dict["second"] = item["second_residue"]
                            unique_pair_dict["database_phi"] = item["database_phi"]
                            unique_pair_dict["database_psi"] = item["database_psi"]
                            unique_pair_dict["detected_torsions_in_structures"].append(
                                detected_torsions
                            )
                            alreadyAdded = True
            output.append(unique_pair_dict)

        return output


class TorsionVisualiser:
    def __init__(self, outputMasterFolder):
        self.outputMasterFolder = outputMasterFolder
        self.glycan_focused_view_output_folder = os.path.join(
            outputMasterFolder, "glycan_perspective"
        )
        self.structure_focused_view_output_folder = os.path.join(
            outputMasterFolder, "structure_perspective"
        )
        CreateFolder(self.glycan_focused_view_output_folder)
        CreateFolder(self.structure_focused_view_output_folder)

    def plot_single_pair_torsions(self, glycan_description, glycanIndex=None):
        for torsion_pair in glycan_description:
            root_description = torsion_pair["root_descr"]
            outputFolderDescription = (
                f"[{glycanIndex}]__{root_description}"
                if glycanIndex is not None
                else f"{root_description}"
            )
            outputFolder = os.path.join(
                self.glycan_focused_view_output_folder, outputFolderDescription
            )
            CreateFolder(outputFolder)
            glycan_torsions = torsion_pair["detected_torsions"]
            bonds = [d["glycan_bond"] for d in glycan_torsions]
            color_labels = list(set(bonds))
            col_values = sns.color_palette("hls", len(color_labels))
            color_map = dict(zip(color_labels, col_values))
            colors = [color_map[label] for label in bonds]

            fig = plt.figure(figsize=(6.4, 3.6), dpi=300, facecolor="gainsboro")
            ax = fig.add_subplot(111)

            plt.hist2d(
                torsion_pair["database_phi"],
                torsion_pair["database_psi"],
                bins=(180, 180),
                cmap=plt.get_cmap("viridis"),
                range=[[-185, 185], [-185, 185]],
                norm=mcolors.PowerNorm(0.7),
            )
            first_residue_name = torsion_pair["first_residue"]
            second_residue_name = torsion_pair["second_residue"]
            plt.title(
                f"{first_residue_name}-{second_residue_name} in {root_description}"
            )
            plt.axhline(linewidth=0.8, color="white")
            plt.axvline(linewidth=0.8, color="white")
            plt.xlim((-181, 181))
            plt.ylim((-181, 181))
            plt.xlabel("φ(Phi) / °", size=14)
            plt.ylabel("ψ(Psi) / °", size=14)
            plt.rc("xtick", labelsize=10)
            plt.rc("ytick", labelsize=10)
            plt.xticks(range(-180, 181, 60))
            plt.yticks(range(-180, 181, 60))
            ax.set_aspect("equal")
            cbar = plt.colorbar()
            cbar.set_label("Frequency", size=14)

            for index, torsion in enumerate(glycan_torsions):
                currentPhi = torsion["Phi"]
                currentPsi = torsion["Psi"]
                plt.plot(
                    currentPhi,
                    currentPsi,
                    marker="x",
                    markersize=4,
                    path_effects=[
                        pe.Stroke(linewidth=1.3, foreground="w"),
                        pe.Normal(),
                    ],
                    linestyle="None",
                    color=colors[index],
                )
            plt.tight_layout()
            plt.legend(
                prop={"size": 6},
                loc=2,
                ncol=1,
                bbox_to_anchor=(-0.8, 0.8),
                title="Linkage:",
                title_fontsize="xx-small",
                handles=[
                    mlines.Line2D(
                        [],
                        [],
                        color=v,
                        marker="x",
                        markersize=4,
                        path_effects=[
                            pe.Stroke(linewidth=1.3, foreground="w"),
                            pe.Normal(),
                        ],
                        linestyle="None",
                        label=k,
                    )
                    for k, v in color_map.items()
                ],
            )
            imageFileName = (
                f"{first_residue_name}-{second_residue_name}_{root_description}.png"
            )

            plt.savefig(
                os.path.join(outputFolder, imageFileName),
                facecolor=fig.get_facecolor(),
                transparent=True,
            )
            plt.cla()
            plt.clf()
            plt.close(fig)

    def plot_whole_structure_torsions(self, parsed_structure_output):
        for torsion_pair in parsed_structure_output:
            first_residue = torsion_pair["first"]
            second_residue = torsion_pair["second"]
            outputFolderDescription = f"{first_residue}-{second_residue}"
            outputFolder = os.path.join(
                self.structure_focused_view_output_folder, outputFolderDescription
            )
            CreateFolder(outputFolder)
            glycan_torsions = torsion_pair["detected_torsions_in_structures"]
            glycanIndices = [d["glycanIndex"] for d in glycan_torsions]
            color_labels = list(set(glycanIndices))
            col_values = sns.color_palette("hls", len(color_labels))
            color_map = dict(zip(color_labels, col_values))
            colors = [color_map[label] for label in glycanIndices]

            fig = plt.figure(figsize=(6.4, 3.6), dpi=300, facecolor="gainsboro")
            ax = fig.add_subplot(111)

            plt.hist2d(
                torsion_pair["database_phi"],
                torsion_pair["database_psi"],
                bins=(180, 180),
                cmap=plt.get_cmap("viridis"),
                range=[[-185, 185], [-185, 185]],
                norm=mcolors.PowerNorm(0.7),
            )
            plt.title(f"{first_residue}-{second_residue} linkage torsions")
            plt.axhline(linewidth=0.8, color="white")
            plt.axvline(linewidth=0.8, color="white")
            plt.xlim((-181, 181))
            plt.ylim((-181, 181))
            plt.xlabel("φ(Phi) / °", size=14)
            plt.ylabel("ψ(Psi) / °", size=14)
            plt.rc("xtick", labelsize=10)
            plt.rc("ytick", labelsize=10)
            plt.xticks(range(-180, 181, 60))
            plt.yticks(range(-180, 181, 60))
            ax.set_aspect("equal")
            cbar = plt.colorbar()
            cbar.set_label("Frequency", size=14)

            imageFileName = f"{first_residue}-{second_residue}.png"
            plt.tight_layout()
            plt.savefig(
                os.path.join(
                    outputFolder, f"{first_residue}-{second_residue}_reference.png"
                ),
                facecolor=fig.get_facecolor(),
                transparent=True,
            )

            for index, glycan in enumerate(glycan_torsions):
                torsions = glycan["glycan_torsions"]
                for torsion in torsions:
                    currentPhi = torsion["Phi"]
                    currentPsi = torsion["Psi"]
                    plt.plot(
                        currentPhi,
                        currentPsi,
                        marker="x",
                        markersize=4,
                        path_effects=[
                            pe.Stroke(linewidth=1.3, foreground="w"),
                            pe.Normal(),
                        ],
                        linestyle="None",
                        color=colors[index],
                    )

                    plt.tight_layout()

            plt.legend(
                prop={"size": 4},
                loc=2,
                ncol=3,
                bbox_to_anchor=(-0.8, 1.05),
                title="Glycan ID:",
                title_fontsize="xx-small",
                handles=[
                    mlines.Line2D(
                        [],
                        [],
                        color=v,
                        marker="x",
                        markersize=4,
                        path_effects=[
                            pe.Stroke(linewidth=1.3, foreground="w"),
                            pe.Normal(),
                        ],
                        linestyle="None",
                        label=k,
                    )
                    for k, v in color_map.items()
                ],
            )
            plt.savefig(
                os.path.join(outputFolder, imageFileName),
                facecolor=fig.get_facecolor(),
                transparent=True,
            )
            plt.cla()
            plt.clf()
            plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="torsions.py",
        usage="%(prog)s [options]. Python -pdbin 5fjj.pdb",
        description=f"For an input PDB file generates Ramachandran-like plots of torsion for every individual glycan or entire structure for linkages found in privateer_torsion_database.json",
    )
    parser.add_argument(
        "-pdbin",
        action="store",
        default=None,
        required=True,
        dest="user_inputPath",
        help="Path to the input PDB file",
    )

    args = parser.parse_args()
    inputPath = args.user_inputPath

    fileName = os.path.basename(inputPath)[:4]

    print(f"Analyzing '{fileName}' located in '{inputPath}'")

    PRIVATEERRESULTSPATH = os.getenv("PRIVATEERRESULTS")
    if PRIVATEERRESULTSPATH is not None:
        dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
        currentStructureResultsPath = os.path.join(
            PRIVATEERRESULTSPATH, fileName + "__" + dt_string
        )
        os.mkdir(currentStructureResultsPath)
    else:
        raise EnvironmentError(
            "Unable to retrieve 'PRIVATEERResults' environment variable. Please try sourcing the ccp4.envsetup-sh file again"
        )
    structure = pvtcore.GlycosylationComposition_memsafe(inputPath)
    torsionsdb = pvtcore.OfflineTorsionsDatabase()
    total_torsions_in_structure = structure.get_torsions_summary(torsionsdb)
    wrapper = PrivateerTorsionResultsOutputParser(total_torsions_in_structure)

    parsed_structure_output = wrapper.getSortedOutputOfTorsionsForEntireStructure()
    visualizer = TorsionVisualiser(currentStructureResultsPath)
    print(f"Outputting produced figures to {os.path.join(currentStructureResultsPath)}")
    visualizer.plot_whole_structure_torsions(parsed_structure_output)
    for item in total_torsions_in_structure:
        glycan = item["all_torsions_in_structure"]
        glycanIndex = item["glycanIndex"]
        visualizer.plot_single_pair_torsions(glycan, glycanIndex)

    print("Task Completed Successfully!")
