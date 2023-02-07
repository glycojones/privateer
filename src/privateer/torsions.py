import argparse
import json
import multiprocessing
import os
import time
from dataclasses import dataclass, field
from datetime import datetime
from logging import root
from typing import List, Tuple

import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from privateer import privateer_core as pvtcore


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)


@dataclass
class TorsionEntry:
    Phi: float
    Psi: float
    glycan_bond: str
    sugar_1: str
    sugar_2: str
    donor_position: str
    acceptor_position: str
    glycanIndex: int


@dataclass
class GlycanTorsion:
    glycanIndex: int
    WURCS: str
    torsions: List[TorsionEntry]
    root_description: str


@dataclass
class TorsionSet:
    sugar_1: str = ""
    sugar_2: str = ""
    donor_position: str = ""
    acceptor_position: str = ""
    torsions: List = field(default_factory=list)  # List[GlycanTorsion]
    database_phi: List = field(default_factory=list)  # List[Float]
    database_psi: List = field(default_factory=list)  # List[Float]


class PrivateerTorsionResultsOutputParser:

    def __init__(self, privateer_output):
        self.__type = "uninitialized"
        self.__privateer_output = privateer_output
        self.__type = self.__validateInputJSON()
        self.__parsedTorsionsOfStructureOutputForPlotting: List[
            TorsionSet] = None
        if self.__type == "structure":
            self.__parsedTorsionsOfStructureOutputForPlotting = (
                self.__parseOutputIfTorsionsStructureProvided())

    def getTypeOfOutput(self):
        return self.__type

    def getSortedOutputOfTorsionsForEntireStructure(self):
        return self.__parsedTorsionsOfStructureOutputForPlotting

    def __validateInputJSON(self):
        expected_keys_structure = [
            "glycanIndex", "WURCS", "all_torsions_in_structure"
        ]
        expected_keys_glycan = [
            "first_residue",
            "second_residue",
            "first_number",
            "second_number",
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
            raise AttributeError(
                "Imported torsions input from Privateer is empty.")

        return self.__type

    def __parseOutputIfTorsionsStructureProvided(self) -> List[TorsionSet]:

        def collectUniqueResiduePairs(input):
            pairs = []
            for rootItem in input:
                all_torsions = rootItem["all_torsions_in_structure"]
                for item in all_torsions:
                    pair = [
                        item["first_residue"], item["second_residue"],
                        item["first_number"], item["second_number"]
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
            donor_position = unique_pair[2]
            acceptor_position = unique_pair[3]
            alreadyAdded = False

            unique_pair = TorsionSet()

            for rootItem in self.__privateer_output:
                all_torsions = rootItem["all_torsions_in_structure"]
                for item in all_torsions:
                    if (item["first_residue"] == first_unique_residue and
                            item["second_residue"] == second_unique_residue):
                        detected_torsions = {
                            "glycanIndex": rootItem["glycanIndex"],
                            "WURCS": rootItem["WURCS"],
                            "glycan_torsions": item["detected_torsions"],
                        }

                        for index, torsion in enumerate(
                                item["detected_torsions"]):
                            item["detected_torsions"][index][
                                "sugar_1"] = first_unique_residue
                            item["detected_torsions"][index][
                                "sugar_2"] = second_unique_residue

                            item["detected_torsions"][index][
                                "donor_position"] = donor_position
                            item["detected_torsions"][index][
                                "acceptor_position"] = acceptor_position

                            item["detected_torsions"][index][
                                "glycanIndex"] = rootItem["glycanIndex"]

                        glycan_torsions = [
                            TorsionEntry(**data)
                            for data in item["detected_torsions"]
                        ]

                        detected_torsions = GlycanTorsion(
                            glycanIndex=rootItem["glycanIndex"],
                            WURCS=rootItem["WURCS"],
                            torsions=glycan_torsions,
                            root_description=item["root_descr"],
                        )

                        if alreadyAdded:
                            unique_pair.torsions.append(detected_torsions)
                        else:
                            unique_pair = TorsionSet(
                                sugar_1=item["first_residue"],
                                sugar_2=item["second_residue"],
                                donor_position=item['first_number'],
                                acceptor_position=item['second_number'],
                                database_phi=item["database_phi"],
                                database_psi=item["database_psi"],
                            )
                            unique_pair.torsions.append(detected_torsions)
                            alreadyAdded = True

            output.append(unique_pair)
        return output


class TorsionVisualiser:
    stats_cache = {"Sugar 1": "", "Sugar 2": "", "Stats": ()}

    def __init__(self,
                 master_folder: str,
                 combine_legends: bool = False) -> None:
        self.combined_legends = combine_legends
        self._create_dir_structure(master_folder)

    def plot_single_pair_torsions(self, torsion_set: TorsionSet):

        for glycan in torsion_set.torsions:
            outputFolderDescription = f"{torsion_set.sugar_1}-{torsion_set.acceptor_position},{torsion_set.donor_position}-{torsion_set.sugar_2}"
            outputFolder = os.path.join(self.glycan_focused_view_output_folder,
                                        outputFolderDescription)

            CreateFolder(outputFolder)

            for torsion in glycan.torsions:
                self._init_plot()

                title = f"{torsion_set.sugar_1}-{torsion_set.acceptor_position},{torsion_set.donor_position}-{torsion_set.sugar_2} in {glycan.root_description}"

                label = self._get_label(glycan, torsion_set)

                self._draw_base_plot(torsion_set, title)
                points = self._draw_annotations(torsions=[torsion],
                                                colours=None,
                                                labels=[label])
                self._draw_legends(points)

                imageFileName = (
                    f"{outputFolderDescription}__{glycan.root_description}.png"
                )

                self._save_figure(outputFolder, imageFileName)

    def plot_all(self, torsion_set: TorsionSet):

        torsion_list = []
        label_list = []

        for glycan in torsion_set.torsions:
            for torsions in glycan.torsions:
                torsion_list.append(torsions)
                label = self._get_label(glycan, torsion_set)
                label_list.append(label)

        outputFolderDescription = f"{torsion_set.sugar_1}-{torsion_set.acceptor_position},{torsion_set.donor_position}-{torsion_set.sugar_2}"
        outputFolder = os.path.join(self.structure_focused_view_output_folder,
                                    outputFolderDescription)

        CreateFolder(outputFolder)

        col_values = sns.color_palette("hls", len(torsion_list))

        title = f"{torsion_set.sugar_1}-{torsion_set.acceptor_position},{torsion_set.donor_position}-{torsion_set.sugar_2} linkage torsions"

        self._init_plot()
        self._draw_base_plot(torsion_set, title)
        points = self._draw_annotations(torsions=torsion_list,
                                        colours=col_values,
                                        labels=label_list)
        self._draw_legends(points)
        imageFileName = f"{torsion_set.sugar_1}-{torsion_set.acceptor_position},{torsion_set.donor_position}-{torsion_set.sugar_2}_{torsion_set.torsions[0].root_description}.png"
        self._save_figure(output_dir=outputFolder, image_name=imageFileName)

    def _get_label(self, glycan, torsion_set):

        sugar_1_tmp_descr = glycan.root_description.split("[")
        sugar_1_tmp_descr = sugar_1_tmp_descr[1].split("]")
        sugar_1_tmp_descr = sugar_1_tmp_descr[0]

        sugar_2_tmp_descr = glycan.root_description.split("[")
        sugar_2_tmp_descr = sugar_2_tmp_descr[2].split("]")
        sugar_2_tmp_descr = sugar_2_tmp_descr[0]

        label = (
            f"{glycan.glycanIndex} {torsion_set.sugar_1}/{sugar_2_tmp_descr}-{torsion_set.sugar_2}/{sugar_1_tmp_descr}"
            if (torsion_set.sugar_1 == "ASN"
                and torsion_set.sugar_2 == "NAG") else
            f"{glycan.glycanIndex} {torsion_set.sugar_1}/{sugar_1_tmp_descr}-{torsion_set.sugar_2}/{sugar_1_tmp_descr}"
        )

        return label

    def _create_dir_structure(self, master_folder):
        self.glycan_focused_view_output_folder = os.path.join(
            master_folder, "glycan_perspective")
        self.structure_focused_view_output_folder = os.path.join(
            master_folder, "structure_perspective")

        CreateFolder(self.glycan_focused_view_output_folder)
        CreateFolder(self.structure_focused_view_output_folder)

    def _init_plot(self):
        self.fig = plt.figure(figsize=(8, 4.5))
        self.ax = self.fig.add_subplot(111)

    def _draw_base_plot(self, torsion_set: TorsionSet, title):

        rng = (np.array([(-180, 180), (0, 360)]) if
               (torsion_set.sugar_1 == "ASN"
                and torsion_set.sugar_2 == "NAG") else np.array([(-180, 180),
                                                                 (-180, 180)]))
        plt.hist2d(
            torsion_set.database_phi,
            torsion_set.database_psi,
            bins=(180, 180),
            cmap=plt.get_cmap("gist_heat_r"),
            range=rng,
            norm=mcolors.PowerNorm(0.7),
        )

        plt.axhline(linewidth=0.8, color="black")
        plt.axvline(linewidth=0.8, color="black")
        plt.xlim((-180, 180))
        y_lim = ((0, 360) if (torsion_set.sugar_1 == "ASN"
                              and torsion_set.sugar_2 == "NAG") else
                 (-180, 180))
        plt.title(title)
        plt.ylim(y_lim)
        plt.xlabel("φ / °", size=14)
        plt.ylabel("ψ / °", size=14)
        plt.rc("xtick", labelsize=14)
        plt.rc("ytick", labelsize=14)
        self.ax.set_aspect("equal")
        cbar = plt.colorbar()
        cbar.set_label("Frequency", size=14)

    def _draw_annotations(self, torsions: List[TorsionEntry], colours,
                          labels) -> Tuple[List[mlines.Line2D]]:

        if colours == None:
            colours = "blue"

        inliers: List[TorsionEntry] = []
        outliers: List[TorsionEntry] = []
        inlier_points: List[mlines.Line2D] = []
        outlier_points: List[mlines.Line2D] = []

        for torsion_pair in torsions:
            if self._is_outlier(torsion_pair):
                inliers.append(torsion_pair)
            else:
                outliers.append(torsion_pair)

        point_index = 0

        for inlier in inliers:
            (annotation, ) = self.ax.plot(
                inlier.Phi,
                inlier.Psi,
                marker="x",
                path_effects=[
                    pe.Stroke(linewidth=1.3, foreground="w"),
                    pe.Normal(),
                ],
                linestyle="None",
                label=labels[point_index],
                color=colours[point_index],
            )
            point_index += 1
            inlier_points.append(annotation)

        outlier_colours = sns.color_palette("cubehelix", len(outliers))

        for index, outlier in enumerate(outliers):
            (annotation, ) = self.ax.plot(
                outlier.Phi,
                outlier.Psi,
                marker="*",
                path_effects=[
                    pe.Stroke(linewidth=1.3, foreground="w"),
                    pe.Normal(),
                ],
                linestyle="None",
                label=labels[point_index],
                color=outlier_colours[index],
            )
            point_index += 1
            outlier_points.append(annotation)

        return (inlier_points, outlier_points)

    def _draw_legends(self, points: Tuple[List[mlines.Line2D],
                                          List[mlines.Line2D]]):

        inlier_points, outlier_points = points

        if self.combined_legends:
            combinded_points = inlier_points + outlier_points

            combined_legends = self.ax.legend(
                title="Linkages:",
                title_fontsize="small",
                prop={"size": 10}
                if len(combinded_points) < 5 else {"size": 6},
                handles=combinded_points,
                bbox_to_anchor=(-0.3, 0.5),
                loc="right",
                labelcolor=[
                    "red" if index >= len(inlier_points) else "black"
                    for index, point in enumerate(combinded_points)
                ],
                ncol=3,
            )

            self.ax.add_artist(combined_legends)
            return

        if inlier_points:
            linkage_legend = self.ax.legend(
                title="Linkages:",
                title_fontsize="small",
                prop={"size": 10} if len(inlier_points) < 1 else {"size": 6},
                handles=inlier_points,
                bbox_to_anchor=(-0.3, 0.5),
                ncol=3,
                loc="upper right" if outlier_points else "right",  #
            )
            self.ax.add_artist(linkage_legend)

        if outlier_points:
            outlier_legend = self.ax.legend(
                title="Outliers:",
                title_fontsize="small",
                prop={"size": 10} if len(outlier_points) < 1 else {"size": 6},
                handles=outlier_points,
                bbox_to_anchor=(-0.3, 0.5),
                loc="lower right" if inlier_points else "right",
                ncol=3,
                labelcolor="red",
            )

            self.ax.add_artist(outlier_legend)

    def _show(self):
        plt.show()

    def _save_figure(self, output_dir, image_name):
        plt.savefig(
            os.path.join(output_dir, image_name),
            facecolor=self.fig.get_facecolor(),
            transparent=True,
            dpi=600,
            bbox_inches="tight",
        )
        plt.close(self.fig)

    def _is_outlier(self, torsion: TorsionEntry) -> bool:

        if (self.stats_cache["Sugar 1"] == torsion.sugar_1
                and self.stats_cache["Sugar 2"] == torsion.sugar_2):
            (phi_mean, psi_mean, phi_std, psi_std) = self.stats_cache["Stats"]
        else:
            if os.getenv("PRIVATEERDATA", None) is not None:
                ROOTPATH = os.getenv("PRIVATEERDATA", None)
            else:
                ROOTPATH = os.getenv("CLIBD", None)
                if ROOTPATH is None:
                    raise EnvironmentError(
                        "Unable to retrieve 'PRIVATEERDATA' nor 'CLIBD' environment variable. Please try sourcing the ccp4.envsetup-sh file again"
                    )
                ROOTPATH = os.path.join(ROOTPATH, "privateer_data")
            (phi_mean, psi_mean, phi_std,
             psi_std) = self._get_outliers_from_file(
                 file_source=os.path.join(ROOTPATH, "linkage_torsions",
                                          "privateer_torsion_statistics.json"),
                 sugar_1=torsion.sugar_1,
                 sugar_2=torsion.sugar_2,
             )
            self.stats_cache["Stats"] = (phi_mean, psi_mean, phi_std, psi_std)
            self.stats_cache["Sugar 1"] = torsion.sugar_1
            self.stats_cache["Sugar 2"] = torsion.sugar_2

        number_of_stdevs = 2
        lower_phi_range = phi_mean - (number_of_stdevs * phi_std)
        higher_phi_range = phi_mean + (number_of_stdevs * phi_std)
        lower_psi_range = psi_mean - (number_of_stdevs * psi_std)
        higher_psi_range = psi_mean + (number_of_stdevs * psi_std)

        if (lower_phi_range <= torsion.Phi <= higher_phi_range) and (
                lower_psi_range <= torsion.Psi <= higher_psi_range):
            return True
        else:
            return False

    # Could obtain this from C++ similar to how the database torsions are passed up.
    def _get_outliers_from_file(self, file_source, sugar_1, sugar_2):
        with open(file_source, "r") as stats_file:
            stats = json.load(stats_file)

            for stat in stats:
                if stat["Linkage"] == f"{sugar_1}-{sugar_2}":
                    return (
                        stat["Stats"]["Phi Mean"],
                        stat["Stats"]["Psi Mean"],
                        stat["Stats"]["Phi Std"],
                        stat["Stats"]["Psi Std"],
                    )


if __name__ == "__main__":

    t_0 = time.time()

    parser = argparse.ArgumentParser(
        prog="torsions.py",
        usage="%(prog)s [options]. Python -pdbin 5fjj.pdb",
        description=
        f"For an input PDB file generates Ramachandran-like plots of torsion for every individual glycan or entire structure for linkages found in privateer_torsion_database.json",
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

    PRIVATEERRESULTSPATH = os.getenv("PRIVATEERRESULTS", None)
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")

    if PRIVATEERRESULTSPATH is not None:
        currentStructureResultsPath = os.path.join(PRIVATEERRESULTSPATH,
                                                   fileName + "__" + dt_string)
    else:
        currentStructureResultsPath = fileName + "__" + dt_string

    os.mkdir(currentStructureResultsPath)
    structure = pvtcore.GlycosylationComposition_memsafe(inputPath)
    torsionsdb = pvtcore.OfflineTorsionsDatabase()
    total_torsions_in_structure = structure.get_torsions_summary(torsionsdb)

    wrapper = PrivateerTorsionResultsOutputParser(total_torsions_in_structure)

    parsed_structure_output = wrapper.getSortedOutputOfTorsionsForEntireStructure(
    )
    t_1 = time.time()

    print(f"Time taken to initialise - {t_1 - t_0} seconds ")

    visualizer = TorsionVisualiser(master_folder=currentStructureResultsPath,
                                   combine_legends=False)

    t_2 = time.time()

    print(f"Time taken to create visualisation -  {t_2 - t_1} seconds")

    multiprocessing_enabled = True

    processes = []

    t_3 = time.time()

    if multiprocessing_enabled:
        for torsion_set in parsed_structure_output:
            p1 = multiprocessing.Process(target=visualizer.plot_all,
                                         args=(torsion_set, ))
            p1.start()
            processes.append(p1)

            p2 = multiprocessing.Process(
                target=visualizer.plot_single_pair_torsions,
                args=(torsion_set, ))
            p2.start()
            processes.append(p2)

        for process in processes:
            process.join()
    else:
        for torsion_set in parsed_structure_output:
            visualizer.plot_all(torsion_set=torsion_set)
            visualizer.plot_single_pair_torsions(torsion_set=torsion_set)

    t_4 = time.time()

    print(f"Total time taken -  {t_4 - t_0} seconds")

    print(
        f"Outputting produced figures to {os.path.join(currentStructureResultsPath)}"
    )
    print("Task Completed Successfully!")

# Single thread time taken - 220 seconds / 3 minutes 40 seconds
# Multithread time taken -  82 seconds with 5fjj / 15 seconds with 2wah
