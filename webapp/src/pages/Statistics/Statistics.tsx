import React, { lazy, type ReactElement } from 'react';
import { Information } from '../../shared/Information/Information.tsx';
import StatisticsHeader from '../../layouts/StatisticsHeader.tsx';
const Footer = lazy(async () => await import('../../layouts/Footer.tsx'));
const BorderElement = lazy(
    async () => await import('../../layouts/BorderElement.tsx')
);

export default function DatabaseSection(): ReactElement {
    return (
        <>
            <StatisticsHeader />
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
