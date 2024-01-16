import React, { lazy, type ReactElement, useEffect, useState } from 'react';
import { Information } from '../../shared/Information/Information.tsx';
import {type StatisticsHeaderProps} from '../../interfaces/types';
import StatisticsHeader from "../../layouts/StatisticsHeader.tsx"
const Footer = lazy(async () => await import('../../layouts/Footer.tsx'));
const BorderElement = lazy(
    async () => await import('../../layouts/BorderElement.tsx')
);

export default function DatabaseSection(): ReactElement {
    const [resetApp, setResetApp] = useState<boolean>(false);

    const mainProps: StatisticsHeaderProps = {
        resetApp: resetApp, setResetApp: setResetApp
    }
    return (
        <>
            <StatisticsHeader {...mainProps} />
            <BorderElement
                topColor={'#D6D9E5'}
                bottomColor={'#F4F9FF'}
                reverse={false}
            ></BorderElement>
            <Information />
            <BorderElement
                topColor={'#F4F9FF'}
                bottomColor={'#D6D9E5'}
                reverse={true}
            ></BorderElement>
            <Footer></Footer>
        </>
    );
}
