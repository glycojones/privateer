import { libcootApi } from "../../src/types/libcoot"
import { emscriptem } from "../../src/types/emscriptem"

// @ts-ignore
importScripts('./wasm/moorhen.js')
// @ts-ignore
importScripts('./wasm/web_example.js')

let cootModule: libcootApi.CootModule;
let molecules_container: libcootApi.MoleculesContainerJS;
let ccp4Module: any;

const guid = () => {
    var d = Date.now();
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
        var r = (d + Math.random() * 16) % 16 | 0;
        d = Math.floor(d / 16);
        return (c === 'x' ? r : (r & 0x3 | 0x8)).toString(16);
    });
    return uuid;
}

// @ts-ignore
let print = (stuff) => {
    console.log(stuff)
    postMessage({ consoleMessage: JSON.stringify(stuff) })
}

const instancedMeshToMeshData = (instanceMesh: libcootApi.InstancedMeshT, perm: boolean, toSpheres: boolean = false, maxZSize: number = 10000.): libcootApi.InstancedMeshJS => {
    //maxZSize is arguably a hack to deal with overlong bonds. It is set to 5 incall to this function.

    let totIdxs: number[][] = []
    let totPos: number[][] = []
    let totNorm: number[][] = []
    let totInstance_sizes: number[][] = []
    let totInstance_colours: number[][] = []
    let totInstance_origins: number[][] = []
    let totInstance_orientations: number[][] = []
    let totInstanceUseColours: boolean[] = []
    let totInstancePrimTypes: string[] = []

    const geom = instanceMesh.geom
    const markup = instanceMesh.markup
    const geomSize = geom.size()
    for (let i = 0; i < geomSize; i++) {
        let thisToSpheres = toSpheres;
        let thisIdxs: number[] = []
        let thisPos: number[] = []
        let thisNorm: number[] = []
        let thisInstance_sizes: number[] = []
        let thisInstance_colours: number[] = []
        let thisInstance_origins: number[] = []
        let thisInstance_orientations: number[] = []
        const inst = geom.get(i);
        if (inst.name === "spherical-atoms") thisToSpheres = true;
        const vertices = inst.vertices;
        const triangles = inst.triangles;
        const trianglesSize = triangles.size()
        for (let i = 0; i < trianglesSize; i++) {
            const triangle = triangles.get(i)
            const idxs = triangle.point_id
            if (perm) {
                thisIdxs.push(idxs[0])
                thisIdxs.push(idxs[2])
                thisIdxs.push(idxs[1])
            } else {
                thisIdxs.push(idxs[0])
                thisIdxs.push(idxs[1])
                thisIdxs.push(idxs[2])
            }
        }
        triangles.delete()

        const verticesSize = vertices.size()
        for (let i = 0; i < verticesSize; i++) {
            const vert = vertices.get(i);
            const vertPos = vert.pos
            thisPos.push(vertPos[0])
            thisPos.push(vertPos[1])
            thisPos.push(vertPos[2])
            const vertNormal = vert.normal
            thisNorm.push(vertNormal[0])
            thisNorm.push(vertNormal[1])
            thisNorm.push(vertNormal[2])
            vert.delete()
        }
        vertices.delete()

        const As = inst.instancing_data_A;
        const Asize = As.size();

        if (Asize > 0) {
            for (let j = 0; j < Asize; j++) {
                const inst_data = As.get(j)

                const instDataPosition = inst_data.position
                thisInstance_origins.push(instDataPosition[0])
                thisInstance_origins.push(instDataPosition[1])
                thisInstance_origins.push(instDataPosition[2])

                const instDataColour = inst_data.colour
                thisInstance_colours.push(instDataColour[0])
                thisInstance_colours.push(instDataColour[1])
                thisInstance_colours.push(instDataColour[2])
                thisInstance_colours.push(instDataColour[3])

                const instDataSize = inst_data.size
                thisInstance_sizes.push(instDataSize[0])
                thisInstance_sizes.push(instDataSize[1])
                thisInstance_sizes.push(instDataSize[2])

                thisInstance_orientations.push(...[
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    0.0, 0.0, 0.0, 1.0,
                ])

                inst_data.delete()
            }
        }
        As.delete()

        const Bs = inst.instancing_data_B;
        const Bsize = Bs.size();
        if (Bsize > 0) {
            for (let j = 0; j < Bsize; j++) {
                const inst_data = Bs.get(j)

                const instDataSize = inst_data.size
                if (instDataSize[2] > maxZSize) continue;
                thisInstance_sizes.push(instDataSize[0])
                thisInstance_sizes.push(instDataSize[1])
                thisInstance_sizes.push(instDataSize[2])

                const instDataPosition = inst_data.position
                thisInstance_origins.push(instDataPosition[0])
                thisInstance_origins.push(instDataPosition[1])
                thisInstance_origins.push(instDataPosition[2])

                const instDataColour = inst_data.colour
                thisInstance_colours.push(instDataColour[0])
                thisInstance_colours.push(instDataColour[1])
                thisInstance_colours.push(instDataColour[2])
                thisInstance_colours.push(instDataColour[3])

                const instDataOrientation = inst_data.orientation

                thisInstance_orientations.push(instDataOrientation[0][0])
                thisInstance_orientations.push(instDataOrientation[0][1])
                thisInstance_orientations.push(instDataOrientation[0][2])
                thisInstance_orientations.push(instDataOrientation[0][3])

                thisInstance_orientations.push(instDataOrientation[1][0])
                thisInstance_orientations.push(instDataOrientation[1][1])
                thisInstance_orientations.push(instDataOrientation[1][2])
                thisInstance_orientations.push(instDataOrientation[1][3])

                thisInstance_orientations.push(instDataOrientation[2][0])
                thisInstance_orientations.push(instDataOrientation[2][1])
                thisInstance_orientations.push(instDataOrientation[2][2])
                thisInstance_orientations.push(instDataOrientation[2][3])

                thisInstance_orientations.push(instDataOrientation[3][0])
                thisInstance_orientations.push(instDataOrientation[3][1])
                thisInstance_orientations.push(instDataOrientation[3][2])
                thisInstance_orientations.push(instDataOrientation[3][3])

                inst_data.delete()
            }
        }
        Bs.delete()
        inst.delete()

        totNorm.push(thisNorm)
        totPos.push(thisPos)
        totIdxs.push(thisIdxs)
        totInstance_sizes.push(thisInstance_sizes)
        totInstance_origins.push(thisInstance_origins)
        totInstance_orientations.push(thisInstance_orientations)
        totInstance_colours.push(thisInstance_colours)
        totInstanceUseColours.push(true)
        if (thisToSpheres)
            totInstancePrimTypes.push("PERFECT_SPHERES")
        else
            totInstancePrimTypes.push("TRIANGLES")

    }

    geom.delete()
    const simpleMeshData = simpleMeshToMeshData(markup) // simpleMeshToMeshData should do the "delete"
    instanceMesh.delete()

    if (simpleMeshData.idx_tri.length > 0 && simpleMeshData.idx_tri[0].length > 0 && simpleMeshData.idx_tri[0][0].length > 0) {
        if (toSpheres) {
            return {
                prim_types: [totInstancePrimTypes, simpleMeshData.prim_types[0]],
                idx_tri: [totIdxs, simpleMeshData.idx_tri[0]],
                vert_tri: [totInstance_origins, simpleMeshData.vert_tri[0]],
                norm_tri: [totNorm, simpleMeshData.norm_tri[0]],
                col_tri: [totInstance_colours, simpleMeshData.col_tri[0]],
                instance_use_colors: [totInstanceUseColours, null],
                instance_sizes: [totInstance_sizes, null],
                instance_origins: [totInstance_origins, null],
                instance_orientations: [totInstance_orientations, null]
            }
        } else {
            return {
                prim_types: [totInstancePrimTypes, simpleMeshData.prim_types[0]],
                idx_tri: [totIdxs, simpleMeshData.idx_tri[0]],
                vert_tri: [totPos, simpleMeshData.vert_tri[0]],
                norm_tri: [totNorm, simpleMeshData.norm_tri[0]],
                col_tri: [totInstance_colours, simpleMeshData.col_tri[0]],
                instance_use_colors: [totInstanceUseColours, null],
                instance_sizes: [totInstance_sizes, null],
                instance_origins: [totInstance_origins, null],
                instance_orientations: [totInstance_orientations, null]
            }
        }
    } else {
        return {
            prim_types: [totInstancePrimTypes],
            idx_tri: [totIdxs],
            vert_tri: [totPos],
            norm_tri: [totNorm],
            col_tri: [totInstance_colours],
            instance_use_colors: [totInstanceUseColours],
            instance_sizes: [totInstance_sizes],
            instance_origins: [totInstance_origins],
            instance_orientations: [totInstance_orientations]
        }
    }
}

const simpleMeshToMeshData = (simpleMesh: libcootApi.SimpleMeshT, perm: boolean = false): libcootApi.SimpleMeshJS => {
    const vertices = simpleMesh.vertices;
    const triangles = simpleMesh.triangles;
    let totIdxs: number[] = [];
    let totPos: number[] = [];
    let totNorm: number[] = [];
    let totCol: number[] = [];

    const trianglesSize = triangles.size()
    for (let i = 0; i < trianglesSize; i++) {
        const triangle = triangles.get(i)
        const idxs = triangle.point_id;
        if (perm)
            totIdxs.push(...[idxs[0], idxs[2], idxs[1]]);
        else
            totIdxs.push(...idxs);
    }
    triangles.delete()

    const verticesSize = vertices.size()
    for (let i = 0; i < verticesSize; i++) {
        const vert = vertices.get(i);
        const vertPos = vert.pos
        const vertNormal = vert.normal
        const vertColor = vert.color
        totPos.push(...vertPos);
        if (perm)
            totNorm.push(...[-vertNormal[0], -vertNormal[1], -vertNormal[2]]);
        else
            totNorm.push(...vertNormal);
        totCol.push(...vertColor);
        vert.delete()
    }
    vertices.delete()

    simpleMesh.delete()

    return {
        prim_types: [["TRIANGLES"]],
        idx_tri: [[totIdxs]],
        vert_tri: [[totPos]],
        norm_tri: [[totNorm]],
        col_tri: [[totCol]]
    };
}

const SuperposeResultsToJSArray = (superposeResults: libcootApi.SuperposeResultsT): libcootApi.SuperposeResultsJS => {
    const alignedPairsVec = superposeResults.aligned_pairs
    const alignedPairsVecSize = alignedPairsVec.size()
    let alignedPairsData: { reference: libcootApi.ValidationInformationJS, moving: libcootApi.ValidationInformationJS }[] = []
    
    for (let i = 0; i < alignedPairsVecSize; i++) {
        const alignedPairs = alignedPairsVec.get(i)
        const refResidueData = alignedPairs.first
        const refResidueSpec = refResidueData.residue_spec
        const movResidueData = alignedPairs.second
        const movResidueSpec = movResidueData.residue_spec
        const currentPairData = {
            reference: {
                chainId: refResidueSpec.chain_id,
                insCode: refResidueSpec.ins_code,
                seqNum: refResidueSpec.res_no,
                restype: "UNK",
                value: refResidueData.function_value,
                label: refResidueData.label
            },
            moving: {
                chainId: movResidueSpec.chain_id,
                insCode: movResidueSpec.ins_code,
                seqNum: movResidueSpec.res_no,
                restype: "UNK",
                value: movResidueData.function_value,
                label: movResidueData.label,
            }
        }
        
        movResidueData.delete()
        movResidueSpec.delete()
        refResidueData.delete()
        refResidueSpec.delete()

        alignedPairsData.push(currentPairData)
    }
    
    alignedPairsVec.delete()
    
    return {
        referenceSequence: superposeResults.alignment.first,
        movingSequence: superposeResults.alignment.second,
        supperposeInfo: superposeResults.superpose_info,
        alignedPairsData: alignedPairsData
    }
}

const colourRulesToJSArray = (colourRulesArray: emscriptem.vector<libcootApi.PairType<string, string>>) => {
    let returnResult: libcootApi.PairType<string, string>[] = []
    const colourRulesSize = colourRulesArray.size()
    for (let i = 0; i < colourRulesSize; i++) {
        const rule = colourRulesArray.get(i)
        returnResult.push(rule)
    }
    colourRulesArray.delete()
    return returnResult;
}

const floatArrayToJSArray = (floatArray: emscriptem.vector<number>) => {
    let returnResult: number[] = []
    const floatArraySize = floatArray.size()
    for (let i = 0; i < floatArraySize; i++) {
        const f = floatArray.get(i)
        returnResult.push(f);
    }
    floatArray.delete()
    return returnResult;
}

const mapMoleculeCentreInfoToJSObject = (mapMoleculeCentreInfo: libcootApi.MapMoleculeCentreInfo): libcootApi.MapMoleculeCentreInfoJS => {
    //Takes a coot::util::map_molecule_centre_info and returns a javascript object that resembles it
    //Disposes of the coordOrth
    const updatedCentre = mapMoleculeCentreInfo.updated_centre
    let returnResult = {
        updated_centre: [
            updatedCentre.x(),
            updatedCentre.y(),
            updatedCentre.z()
        ] as [number, number, number],
        success: mapMoleculeCentreInfo.success,
        suggested_radius: mapMoleculeCentreInfo.suggested_radius,
        suggested_contour_level: mapMoleculeCentreInfo.suggested_contour_level
    }
    updatedCentre.delete()
    return returnResult;
}

const intArrayToJSArray = (intArray: emscriptem.vector<number>) => {
    let returnResult: number[] = []
    const intArraySize = intArray.size()
    for (let i = 0; i < intArraySize; i++) {
        const f = intArray.get(i)
        returnResult.push(f);
    }
    intArray.delete()
    return returnResult;
}

const stringArrayToJSArray = (stringArray: emscriptem.vector<string>) => {
    let returnResult: string[] = []
    const stringArraySize = stringArray.size()
    for (let i = 0; i < stringArraySize; i++) {
        const s = stringArray.get(i)
        returnResult.push(s);
    }
    stringArray.delete()
    return returnResult;
}

const symmetryToJSData = (symmetryDataPair: libcootApi.PairType<libcootApi.SymmetryData, emscriptem.vector<number[][]>>) => {
    let result: { x: number; y: number; z: number; asString: string; isym: number; us: number; ws: number; vs: number; matrix: number[][]; }[]= []
    const symmetryData = symmetryDataPair.first
    const symmetryMatrices = symmetryDataPair.second
    const cell = symmetryData.cell
    const symm_trans = symmetryData.symm_trans
    const symmetrySize = symm_trans.size()

    for (let i = 0; i < symmetrySize; i++) {
        const currentSymmetry = symm_trans.get(i)
        const symTransT = currentSymmetry.first
        const cellTranslation = currentSymmetry.second
        const currentSymmMat = symmetryMatrices.get(i)

        result.push({
            x: symTransT.x(),
            y: symTransT.y(),
            z: symTransT.z(),
            asString: symTransT.symm_as_string,
            isym: symTransT.isym(),
            us: cellTranslation.us,
            ws: cellTranslation.ws,
            vs: cellTranslation.vs,
            matrix: currentSymmMat
        })
        symTransT.delete()
    }

    cell.delete()
    symm_trans.delete()
    symmetryMatrices.delete()
    symmetryData.delete()
    return result
}

const mmrrccStatsToJSArray = (mmrrccStats: libcootApi.PairType<emscriptem.map<libcootApi.DensityCorrelationStatsInfoT, libcootApi.ResidueSpecT>, emscriptem.map<libcootApi.DensityCorrelationStatsInfoT, libcootApi.ResidueSpecT>>): libcootApi.MMRCCStatsJS => {
    const parseStats = (stats: emscriptem.map<libcootApi.DensityCorrelationStatsInfoT, libcootApi.ResidueSpecT>) => {
        let result: {resNum: number; insCode: string; modelNumber: number; chainId: string; n: number; correlation: number; }[] = []
        const residueSpecs = stats.keys()
        const mapSize = residueSpecs.size()
        for (let i = 0; i < mapSize; i++) {
            const residueSpec = residueSpecs.get(i)
            const densityCorrStat = stats.get(residueSpec)
            result.push({
                resNum: residueSpec.res_no,
                insCode: residueSpec.ins_code,
                modelNumber: residueSpec.model_number,
                chainId: residueSpec.chain_id,
                n: densityCorrStat.n,
                correlation: densityCorrStat.correlation()
            })
            residueSpec.delete()
            densityCorrStat.delete()
        }
        residueSpecs.delete()
        return result
    }

    const first = mmrrccStats.first
    const second = mmrrccStats.second

    const returnResult = {
        "All atoms": parseStats(first),
        "Side-chains": parseStats(second)
    }

    first.delete()
    second.delete()
    return returnResult
}

const atomSpecToJSArray = (atomSpecs: emscriptem.vector<libcootApi.AtomSpecT>): libcootApi.AtomSpecJS[] => {
    let returnResult: libcootApi.AtomSpecJS[] = []
    const atomsSize = atomSpecs.size()
    for (let ic = 0; ic < atomsSize; ic++) {
        const atom = atomSpecs.get(ic)
        returnResult.push({
            chain_id: atom.chain_id,
            res_no: atom.res_no,
            ins_code: atom.ins_code,
            atom_name: atom.atom_name,
            alt_conf: atom.alt_conf,
            int_user_data: atom.int_user_data,
            float_user_data: atom.float_user_data,
            string_user_data: atom.string_user_data,
            model_number: atom.model_number
        })
        atom.delete()
    }
    atomSpecs.delete()
    return returnResult
}

const residueSpecToJSArray = (residueSpecs: emscriptem.vector<libcootApi.ResidueSpecT>): libcootApi.ResidueSpecJS[] => {
    let returnResult: { resNum: number; insCode: string; modelNumber: number; chainId: string; }[] = []
    const residuesSize = residueSpecs.size()
    for (let ic = 0; ic < residuesSize; ic++) {
        const residue = residueSpecs.get(ic)
        returnResult.push({
            resNum: residue.res_no,
            insCode: residue.ins_code,
            modelNumber: residue.model_number,
            chainId: residue.chain_id
        })
        residue.delete()
    }
    residueSpecs.delete()
    return returnResult
}

const validationDataToJSArray = (validationData: libcootApi.ValidationInformationT, chainID: string | null = null): libcootApi.ValidationInformationJS[] => {
    let returnResult: { chainId: string; insCode: string; seqNum: number; restype: string; value: number; }[] = []
    const cviv = validationData.cviv
    const chainSize = cviv.size()
    for (let chainIndex = 0; chainIndex < chainSize; chainIndex++) {
        const chain = cviv.get(chainIndex)
        if (chainID !== null && chain.chain_id !== chainID) {
            // pass
        } else {
            const resInfo = chain.rviv;
            const resInfoSize = resInfo.size()
            for (let ir = 0; ir < resInfoSize; ir++) {
                const residue = resInfo.get(ir)
                const residueSpec = residue.residue_spec
                returnResult.push({
                    chainId: residueSpec.chain_id,
                    insCode: residueSpec.ins_code,
                    seqNum: residueSpec.res_no,
                    restype: "UNK",
                    value: residue.function_value
                })
                residue.delete()
                residueSpec.delete()
            }
            resInfo.delete()
        }
        chain.delete()
    }
    cviv.delete()
    validationData.delete()
    return returnResult
}

const linesBoxToJSArray = (BoxData: libcootApi.Generic3dLinesBondsBoxT): libcootApi.Generic3dLinesBondsBoxJS[][] => {
    let envdata: {start: {x: number; y: number; z: number; }; end: {x: number; y: number; z: number; }; dist: number; }[][]= []
    const segments = BoxData.line_segments;
    const nSeg = segments.size()
    for (let i = 0; i < nSeg; i++) {
        let thisEnvdata: {start: {x: number; y: number; z: number; }; end: {x: number; y: number; z: number; }; dist: number; }[] = []
        const segsI = segments.get(i)
        const nSegI = segsI.size()
        for (let j = 0; j < nSegI; j++) {
            const seg = segsI.get(j)
            const start = seg.getStart()
            const end = seg.getFinish()
            const ampl = seg.amplitude()
            const startJS = { x: start.x(), y: start.y(), z: start.z() }
            const endJS = { x: end.x(), y: end.y(), z: end.z() }
            thisEnvdata.push({
                start: startJS,
                end: endJS,
                dist: ampl,
            })
            start.delete()
            end.delete()
            seg.delete()
        }
        segsI.delete()
        envdata.push(thisEnvdata)
    }
    segments.delete()
    BoxData.delete()

    return envdata
}

const vectorHBondToJSArray = (HBondData: emscriptem.vector<libcootApi.MoorhenHBond>): libcootApi.HBondJS[] => {
    let hbdata: libcootApi.HBondJS[] = []
    const hbondDataSize = HBondData.size()
    for (let ib = 0; ib < hbondDataSize; ib++) {
        const hb = HBondData.get(ib)
        hbdata.push({
            hb_hydrogen: hb.hb_hydrogen,
            donor: hb.donor,
            acceptor: hb.acceptor,
            donor_neigh: hb.donor_neigh,
            acceptor_neigh: hb.acceptor_neigh,
            angle_1: hb.angle_1,
            angle_2: hb.angle_2,
            angle_3: hb.angle_3,
            dist: hb.dist,
            ligand_atom_is_donor: hb.ligand_atom_is_donor,
            hydrogen_is_ligand_atom: hb.hydrogen_is_ligand_atom,
            bond_has_hydrogen_flag: hb.bond_has_hydrogen_flag,
        })
    }
    HBondData.delete()
    return hbdata
}

const interestingPlaceDataToJSArray = (interestingPlaceData: emscriptem.vector<libcootApi.InterestingPlaceT>): libcootApi.InterestingPlaceDataJS[] => {
    let returnResult: { 
        modelNumber: number;
        chainId: string;
        insCode: string;
        resNum: number;
        featureType: string;
        featureValue: number;
        buttonLabel: string;
        badness: number;
        coordX: number;
        coordY: number;
        coordZ: number;
     }[] = [];

    const interestingPlaceDataSize = interestingPlaceData.size()
    for (let ir = 0; ir < interestingPlaceDataSize; ir++) {
        const residue = interestingPlaceData.get(ir)
        const residueSpec = residue.residue_spec
        returnResult.push({
            modelNumber: residueSpec.model_number,
            chainId: residueSpec.chain_id,
            insCode: residueSpec.ins_code,
            resNum: residueSpec.res_no,
            featureType: residue.feature_type,
            featureValue: residue.feature_value,
            buttonLabel: residue.button_label,
            badness: residue.badness,
            coordX: residue.x,
            coordY: residue.y,
            coordZ: residue.z
        })
        residue.delete()
        residueSpec.delete()
    }
    interestingPlaceData.delete()
    return returnResult
}

const histogramInfoToJSData = (histogramInfo: libcootApi.HistogramInfo): libcootApi.HistogramInfoJS => {
    const histogramCounts = histogramInfo.counts
    const counts = intArrayToJSArray(histogramCounts)
    const result: libcootApi.HistogramInfoJS = {
        counts: counts,
        bin_width: histogramInfo.bin_width,
        base: histogramInfo.base
    }
    return result
}

const autoReadMtzInfoVectToJSArray = (autoReadMtzInfoArray: emscriptem.vector<libcootApi.AutoReadMtzInfo>): libcootApi.AutoReadMtzInfoJS[] => {
    let returnResult: {idx: number;
        F: string;
        phi: string;
        w: string;
        F_obs: string;
        sigF_obs: string;
        Rfree: string;
        weights_used: boolean;
    }[] = []

    const autoReadMtzInfoArraySize = autoReadMtzInfoArray.size()
    for(let i = 0; i < autoReadMtzInfoArraySize; i++) {
        const autoReadMtzInfo = autoReadMtzInfoArray.get(i)
        returnResult.push({
            idx: autoReadMtzInfo.idx,
            F: autoReadMtzInfo.F,
            phi: autoReadMtzInfo.phi,
            w: autoReadMtzInfo.w,
            Rfree: autoReadMtzInfo.Rfree,
            F_obs: autoReadMtzInfo.F_obs,
            sigF_obs: autoReadMtzInfo.sigF_obs,
            weights_used: autoReadMtzInfo.weights_used,
        })
        autoReadMtzInfo.delete()
    }
    autoReadMtzInfoArray.delete()
    
    return returnResult
}

const ramachandranDataToJSArray = (ramachandraData: emscriptem.vector<libcootApi.CootPhiPsiProbT>, chainID: string): libcootApi.RamaDataJS[] => {
    let returnResult: { chainId: string; insCode: string; seqNum: number; restype: string; isOutlier: boolean; phi: number; psi: number; is_pre_pro: boolean; }[] = [];
    const ramachandraDataSize = ramachandraData.size()
    for (let ir = 0; ir < ramachandraDataSize; ir++) {
        const residue = ramachandraData.get(ir)
        const phiPsi = residue.phi_psi
        if (phiPsi.chain_id === chainID) {
            returnResult.push({
                chainId: phiPsi.chain_id,
                insCode: phiPsi.ins_code,
                seqNum: phiPsi.residue_number,
                restype: residue.residue_name(),
                isOutlier: !residue.is_allowed_flag,
                phi: phiPsi.phi(),
                psi: phiPsi.psi(),
                is_pre_pro: residue.residue_name() === 'PRO'
            })
        }
        residue.delete()
        phiPsi.delete()
    }
    ramachandraData.delete()
    return returnResult
}

const simpleMeshToLineMeshData = (simpleMesh: libcootApi.SimpleMeshT, normalLighting: boolean): libcootApi.SimpleMeshJS => {
    const vertices = simpleMesh.vertices;
    const triangles = simpleMesh.triangles;
    let totIdxs: number[] = [];
    let totPos: number[] = [];
    let totNorm: number[] = [];
    let totCol: number[] = [];

    const trianglesSize = triangles.size()
    for (let i = 0; i < trianglesSize; i++) {
        const triangle = triangles.get(i)
        const idxs = triangle.point_id;
        totIdxs.push(...[idxs[0], idxs[1], idxs[0], idxs[2], idxs[1], idxs[2]]);
    }
    triangles.delete()

    const verticesSize = vertices.size()
    for (let i = 0; i < verticesSize; i++) {
        const vert = vertices.get(i);
        totPos.push(...vert.pos);
        totNorm.push(...vert.normal);
        totCol.push(...vert.color);
        vert.delete()
    }
    vertices.delete()

    simpleMesh.delete()

    if (normalLighting)
        return { prim_types: [["NORMALLINES"]], useIndices: [[true]], idx_tri: [[totIdxs]], vert_tri: [[totPos]], additional_norm_tri: [[totNorm]], norm_tri: [[totNorm]], col_tri: [[totCol]] };
    else
        return { prim_types: [["LINES"]], useIndices: [[true]], idx_tri: [[totIdxs]], vert_tri: [[totPos]], norm_tri: [[totNorm]], col_tri: [[totCol]] };

}

const auto_open_mtz = (mtzData: ArrayBufferLike) => {
    const theGuid = guid()
    const asUint8Array = new Uint8Array(mtzData)
    cootModule.FS_createDataFile(".", `${theGuid}.mtz`, asUint8Array, true, true);
    const tempFilename = `./${theGuid}.mtz`
    const result = molecules_container.auto_read_mtz(tempFilename)
    cootModule.FS_unlink(tempFilename)
    return result
}

const replace_map_by_mtz_from_file = (imol: number, mtzData: ArrayBufferLike, selectedColumns: { F: string; PHI: string; }) => {
    const theGuid = guid()
    const tempFilename = `./${theGuid}.mtz`
    const asUint8Array = new Uint8Array(mtzData)
    cootModule.FS_createDataFile(".", tempFilename, asUint8Array, true, true);
    const readMtzArgs: [number, string, string, string, string, boolean] = [imol, tempFilename, selectedColumns.F, selectedColumns.PHI, "", false]
    const result = molecules_container.replace_map_by_mtz_from_file(...readMtzArgs)
    cootModule.FS_unlink(tempFilename)
    return result
}

const new_positions_for_residue_atoms = (molToUpDate: number, residues: libcootApi.AtomInfo[][]) => {
    let success = 0
    const movedResidueVector  = new cootModule.Vectormoved_residue_t()
    residues.forEach(atoms => {
        if (atoms.length > 0) {
            const cidFields = atoms[0].label.split('/')
            let [resNoStr, insCode] = cidFields[3].split(".")
            insCode = insCode ? insCode : ""
            const movedResidue = new cootModule.moved_residue_t(cidFields[2], parseInt(resNoStr), insCode)
            atoms.forEach(atom => {
                const movedAtom = new cootModule.moved_atom_t(atom.name, atom.alt_loc, atom.x, atom.y, atom.z, -1)
                movedResidue.add_atom(movedAtom)
                movedAtom.delete()
            })
            movedResidueVector.push_back(movedResidue)
            movedResidue.delete()
        }
    })
    const thisSuccess = molecules_container.new_positions_for_atoms_in_residues(molToUpDate, movedResidueVector)
    success += thisSuccess
    movedResidueVector.delete()
    return success
}

const read_mtz = (mapData: ArrayBufferLike, name: string, selectedColumns: { F: string; PHI: string; isDifference: boolean; }) => {
    const theGuid = guid()
    const asUint8Array = new Uint8Array(mapData)
    cootModule.FS_createDataFile(".", `${theGuid}.mtz`, asUint8Array, true, true);
    const tempFilename = `./${theGuid}.mtz`
    const read_mtz_args: [string, string, string, string, boolean, boolean] = [tempFilename, selectedColumns.F,
        selectedColumns.PHI, "", false, selectedColumns.isDifference]
    const molNo = molecules_container.read_mtz(...read_mtz_args)
    cootModule.FS_unlink(tempFilename)
    return molNo
}

const associate_data_mtz_file_with_map = (iMol: number, mtzData: { data: ArrayBufferLike; fileName: string; }, F: string, SIGF: string, FREE: string) => {
    const asUint8Array = new Uint8Array(mtzData.data)
    cootModule.FS_createDataFile(".", `${mtzData.fileName}.mtz`, asUint8Array, true, true);
    const mtzFilename = `./${mtzData.fileName}.mtz`
    const args: [number, string, string, string, string] = [iMol, mtzFilename, F, SIGF, FREE]
    molecules_container.associate_data_mtz_file_with_map(...args)
    return mtzFilename
}

const read_ccp4_map = (mapData: ArrayBufferLike, name: string, isDiffMap: boolean) => {
    const theGuid = guid()
    const asUint8Array = new Uint8Array(mapData)
    let fileExtension = ''
    if (name.endsWith('.map.gz')) {
        fileExtension = 'map.gz'
    } else if (name.endsWith('.mrc.gz')) {
        fileExtension = 'mrc.gz'
    }
    cootModule.FS_createDataFile(".", `${theGuid}${fileExtension}`, asUint8Array, true, true);
    const tempFilename = `./${theGuid}${fileExtension}`
    const read_map_args: [string, boolean] = [tempFilename, isDiffMap]
    const molNo = molecules_container.read_ccp4_map(...read_map_args)
    cootModule.FS_unlink(tempFilename)
    return molNo
}

const setUserDefinedBondColours = (imol: number, colours: { cid: string; rgb: [number, number, number] }[], applyColourToNonCarbonAtoms: boolean = false) => {
    let colourMap = new cootModule.MapIntFloat3()
    let indexedResiduesVec = new cootModule.VectorStringUInt_pair()
    
    colours.forEach((colour, index) => {
        colourMap.set(index + 51, colour.rgb)
        const i = { first: colour.cid, second: index + 51 }
        indexedResiduesVec.push_back(i)
    })

    molecules_container.set_user_defined_bond_colours(imol, colourMap)
    molecules_container.set_user_defined_atom_colour_by_selection(imol, indexedResiduesVec, applyColourToNonCarbonAtoms)

    indexedResiduesVec.delete()
    colourMap.delete()
}

const doCootCommand = (messageData: { 
    myTimeStamp: number;
    chainID?: string;
    messageId?: string;
    message: string;
    returnType: string;
    command: string;
    commandArgs: any[];
 }) => {

    const { returnType, command, commandArgs, messageId, myTimeStamp, message } = messageData
    
    try {

        let cootResult
        switch (command) {
            case 'shim_new_positions_for_residue_atoms':
                cootResult = new_positions_for_residue_atoms(...commandArgs as [number, libcootApi.AtomInfo[][]])
                break
            case 'shim_read_mtz':
                cootResult = read_mtz(...commandArgs as [ArrayBufferLike, string, { F: string; PHI: string; isDifference: boolean; }])
                break
            case 'shim_auto_open_mtz':
                cootResult = auto_open_mtz(...commandArgs as [ArrayBuffer])
                break
            case 'shim_read_ccp4_map':
                cootResult = read_ccp4_map(...commandArgs as [ArrayBuffer, string, boolean])
                break
            case 'shim_associate_data_mtz_file_with_map':
                cootResult = associate_data_mtz_file_with_map(...commandArgs as [number, { data: ArrayBufferLike; fileName: string; }, string, string, string])
                break
            case 'shim_replace_map_by_mtz_from_file':
                cootResult = replace_map_by_mtz_from_file(...commandArgs as [number, ArrayBufferLike, { F: string; PHI: string; }])
                break
            case 'shim_set_bond_colours':
                cootResult = setUserDefinedBondColours(...commandArgs as [number, { cid: string; rgb: [number, number, number] }[], boolean])
                break
            default:
                cootResult = molecules_container[command](...commandArgs)
                break
        }

        let returnResult;
        switch (returnType) {
            case 'instanced_mesh_perm':
                returnResult = instancedMeshToMeshData(cootResult, true)
                break;
            case 'histogram_info_t':
                returnResult = histogramInfoToJSData(cootResult)
                break
            case 'symmetry':
                returnResult = symmetryToJSData(cootResult)
                break;
            case 'mmrrcc_stats':
                returnResult = mmrrccStatsToJSArray(cootResult)
                break;
            case 'colour_rules':
                returnResult = colourRulesToJSArray(cootResult)
                break;
            case 'instanced_mesh_perfect_spheres':
                returnResult = instancedMeshToMeshData(cootResult, false, true)
                break;
            case 'instanced_mesh':
                returnResult = instancedMeshToMeshData(cootResult, false, false, 5)
                break;
            case 'mesh_perm':
                returnResult = simpleMeshToMeshData(cootResult, true)
                break;
            case 'mesh':
                returnResult = simpleMeshToMeshData(cootResult)
                break;
            case 'lit_lines_mesh':
                returnResult = simpleMeshToLineMeshData(cootResult, true)
                break;
            case 'lines_mesh':
                returnResult = simpleMeshToLineMeshData(cootResult, false)
                break;
            case 'float_array':
                returnResult = floatArrayToJSArray(cootResult)
                break;
            case 'int_array':
                returnResult = intArrayToJSArray(cootResult)
                break;
            case 'auto_read_mtz_info_array':
                returnResult = autoReadMtzInfoVectToJSArray(cootResult)
                break;
            case 'map_molecule_centre_info_t':
                returnResult = mapMoleculeCentreInfoToJSObject(cootResult)
                break;
            case 'string_array':
                returnResult = stringArrayToJSArray(cootResult)
                break;
            case 'residue_specs':
                returnResult = residueSpecToJSArray(cootResult)
                break;
            case 'atom_specs':
                returnResult = atomSpecToJSArray(cootResult)
                break;
            case 'ramachandran_data':
                returnResult = ramachandranDataToJSArray(cootResult, messageData.chainID as string)
                break;
            case 'validation_data':
                returnResult = validationDataToJSArray(cootResult, messageData.chainID)
                break;
            case 'interesting_places_data':
                returnResult = interestingPlaceDataToJSArray(cootResult)
                break;
            case 'superpose_results':
                returnResult = SuperposeResultsToJSArray(cootResult)
                break
            case 'generic_3d_lines_bonds_box':
                returnResult = linesBoxToJSArray(cootResult)
                break;
            case 'vector_hbond':
                returnResult = vectorHBondToJSArray(cootResult)
                break;
            case 'status_instanced_mesh_pair':
                returnResult = { status: cootResult.first, mesh: instancedMeshToMeshData(cootResult.second, false, false, 5) }
                break;
            case 'void':
                returnResult = null
                break
            case 'status':
            default:
                returnResult = cootResult
                break;
        }

        return {
            messageId, messageSendTime: Date.now(), command: command,
            consoleMessage: `Completed ${command} in ${Date.now() - myTimeStamp} ms`,
            result: { status: 'Completed', result: returnResult }
        }
    }

    catch (err) {
        console.log(err)
        return {
            messageId, myTimeStamp, message,
            consoleMessage: `EXCEPTION RAISED IN ${command}, ${err}`,
            result: { status: 'Exception' }
        }
    }
}

onmessage = function (e) {
    if (e.data.message === 'CootInitialize') {
        createRSRModule({
            locateFile: (file) => `./wasm/${file}`,
            onRuntimeInitialized: () => { },
            mainScriptUrlOrBlob: "moorhen.js",
            print: print,
            printErr: print,
        })
            .then((returnedModule) => {
                postMessage({ consoleMessage: 'Initialized molecules_container', message: e.data.message, messageId: e.data.messageId })
                cootModule = returnedModule;
                molecules_container = new cootModule.molecules_container_js(false)
                molecules_container.set_show_timings(false)
                molecules_container.fill_rotamer_probability_tables()
                molecules_container.set_map_sampling_rate(1.7)
                cootModule.FS.mkdir("COOT_BACKUP");
            })
            .catch((e) => {
                console.log(e)
                print(e);
            });
        
        createCCP4Module({
            locateFile: (file) => `./wasm/${file}`,
            onRuntimeInitialized: () => { },
            mainScriptUrlOrBlob: "web_example.js",
            print: print,
            printErr: print,
        })
            .then((returnedModule) => {
                ccp4Module = returnedModule;
            })
            .catch((e) => {
                console.log(e)
                print(e);
            });
    }

    else if (e.data.message === 'get_mtz_data') {
        const mtzData = cootModule.FS.readFile(e.data.fileName, { encoding: 'binary' });
        postMessage({
            messageId: e.data.messageId,
            myTimeStamp: e.data.myTimeStamp,
            consoleMessage: `Fetched mtz data for map ${e.data.molNo}`,
            message: e.data.message,
            result: { molNo: e.data.molNo, mtzData: mtzData }
        })
    }

    else if (e.data.message === 'get_map') {
        const theGuid = guid()
        const tempFilename = `./${theGuid}.map`
        molecules_container.writeCCP4Map(e.data.molNo, tempFilename)

        const mapData = cootModule.FS.readFile(tempFilename, { encoding: 'binary' }) as Uint8Array;
        cootModule.FS_unlink(tempFilename)
        postMessage({
            messageId: e.data.messageId,
            myTimeStamp: e.data.myTimeStamp,
            consoleMessage: `Fetched map of map ${e.data.molNo}`,
            message: e.data.message,
            result: { molNo: e.data.molNo, mapData: mapData.buffer }
        })
    }

    else if (e.data.message === 'delete_file_name') {
        const result = cootModule.FS_unlink(e.data.fileName)

        postMessage({
            messageId: e.data.messageId,
            myTimeStamp: e.data.myTimeStamp,
            messageTag: "result",
            result: result,
        })

    } else if (e.data.message === 'coot_command_list') {
        const resultList = e.data.commandList.map(command => doCootCommand({...e.data, ...command}))
        postMessage({
            messageId: e.data.messageId, resultList
        })
    }

    if (e.data.message === 'coot_command') {
        const result = doCootCommand(e.data)
        postMessage(result)
    }
}
