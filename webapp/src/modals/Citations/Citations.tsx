import React, { type Dispatch, type SetStateAction } from 'react';
import Modal from 'react-modal';

const customStyles = {
    content: {
        top: '50%',
        left: '50%',
        right: 'auto',
        bottom: 'auto',
        // // marginRight: '-50%',
        transform: 'translate(-50%, -50%)',
        color: 'black',
        borderRadius: 10,
        overflowY: 'auto',
        maxHeight: '100vh',
    },
};

export default function CitationModal(props: {
    modalOpen: boolean;
    setModalOpen: Dispatch<SetStateAction<boolean>>;
}) {
    const handleModalClose = () => {
        props.setModalOpen(false);
    };

    return (
        <div className="max-h-screen overflow-auto">
            <Modal
                isOpen={props.modalOpen}
                onRequestClose={handleModalClose}
                contentLabel="Example Modal"
                style={customStyles}
            >
                <div className="flex flex-row w-full space-between ">
                    <h2 className="mr-16 sm:mr-96">Citations</h2>
                    <button
                        type="button"
                        onClick={handleModalClose}
                        className="text-gray-400 bg-transparent hover:scale-105 rounded-lg text-sm w-8 h-8 ms-auto inline-flex justify-center items-center "
                    >
                        <svg
                            className="w-3 h-3"
                            aria-hidden="true"
                            xmlns="http://www.w3.org/2000/svg"
                            fill="none"
                            viewBox="0 0 14 14"
                        >
                            <path
                                stroke="currentColor"
                                strokeLinecap="round"
                                strokeLinejoin="round"
                                strokeWidth="2"
                                d="m1 1 6 6m0 0 6 6M7 7l6-6M7 7l-6 6"
                            />
                        </svg>
                        <span className="sr-only">Close modal</span>
                    </button>
                </div>

                <div className="mt-4 flex flex-col space-y-2  text-sm sm:text-base">
                    <h3>
                        If you have used the{' '}
                        <b>
                            <i>Privateer Web App</i>
                        </b>
                        , please cite:{' '}
                    </h3>
                    <strong>
                        Dialpuri, J. S., Bagdonas, H., Schofield, L. C., Pham,
                        P. T., Holland, L., Bond, P. S., Sánchez Rodríguez, F.,
                        McNicholas, S. J. & Agirre, J. (2024). Online
                        carbohydrate 3D structure validation with the Privateer
                        web app. Acta Cryst. F80.
                        <p><a href="https://doi.org/10.1107/S2053230X24000359" target="_blank">https://doi.org/10.1107/S2053230X24000359</a></p>
                    </strong>
                </div>

                {/* <div className="mt-4 flex flex-col space-y-2"> */}
                {/*    <h3>If you have used the <b><i>Privateer Database</i></b>, please cite: </h3> */}
                {/*    <span>Dialpuri, J. S., Bagdonas, H., Schofield, L. C., Pham, P. T., Holland, L., Bond, P. S., Sánchez Rodríguez, F., McNicholas, S. J. & Agirre, J. (2024). Online carbohydrate 3D structure validation with the Privateer web app. Acta Cryst. F80. https://doi.org/10.1107/S2053230X24000359</span> */}
                {/* </div> */}

                <div className="mt-4 flex flex-col space-y-2 font text-sm sm:text-base">
                    <h3>
                        The <i>Privateer Web App</i> is possible thanks to these
                        foundational studies:
                    </h3>
                    <span>
                        Dialpuri, J. S., Bagdonas, H., Atanasova, M., Schofield,
                        L. C., Hekkelman, M. L., Joosten, R. P., & Agirre, J.
                        (2023). Analysis and validation of overall N-glycan
                        conformation in Privateer. Acta Crystallographica
                        Section D: Structural Biology, 79(6).
                    </span>
                    <span>
                        Bagdonas, H., Ungar, D., & Agirre, J. (2020). Leveraging
                        glycomics data in glycoprotein 3D structure validation
                        with Privateer. Beilstein Journal of Organic Chemistry,
                        16(1), 2523-2533.
                    </span>
                    <span>
                        McNicholas, S., & Agirre, J. (2017). Glycoblocks: a
                        schematic three-dimensional representation for glycans
                        and their interactions. Acta Crystallographica Section
                        D: Structural Biology, 73(2), 187-194.
                    </span>
                    <span>
                        Atanasova, M., Nicholls, R. A., Joosten, R. P., &
                        Agirre, J. (2022). Updated restraint dictionaries for
                        carbohydrates in the pyranose form. Acta
                        Crystallographica Section D: Structural Biology, 78(4),
                        455-465.
                    </span>
                    <span>
                        Agirre, J., Iglesias-Fernández, J., Rovira, C., Davies,
                        G. J., Wilson, K. S., & Cowtan, K. D. (2015). Privateer:
                        software for the conformational validation of
                        carbohydrate structures. Nature structural & molecular
                        biology, 22(11), 833-834.
                    </span>
                    <span>
                        Agirre, J., Davies, G., Wilson, K., & Cowtan, K. (2015).
                        Carbohydrate anomalies in the PDB. Nature chemical
                        biology, 11(5), 303-303.
                    </span>
                    <span>
                        Atanasova, M., Bagdonas, H., & Agirre, J. (2020).
                        Structural glycobiology in the age of electron
                        cryo-microscopy. Current Opinion in Structural Biology,
                        62, 70-78.
                    </span>
                </div>
            </Modal>
        </div>
    );
}
