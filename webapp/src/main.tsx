import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import './index.css';
import { BrowserRouter } from 'react-router-dom';
import ReactGA from "react-ga4";
import {onCLS, onFID, onLCP} from 'web-vitals';

ReactGA.initialize("G-PGPMR0MEYT");
ReactDOM.createRoot(document.getElementById('root')).render(
    // <React.StrictMode>
    <BrowserRouter>
        <App />
    </BrowserRouter>

    // </React.StrictMode>,
);
const SendAnalytics = ()=> {
    ReactGA.send({
        hitType: "pageview",
        page: window.location.pathname,
    });
}

try {
    onCLS(SendAnalytics);
    onLCP(SendAnalytics);
    onFID(SendAnalytics);
}
catch (err) {
    console.error(err)
}
