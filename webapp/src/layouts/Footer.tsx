import React, { type ReactElement } from 'react';
export default function Footer(): ReactElement {
    return (
        <footer className="bg-gray w-full mt-auto">
            <div className="w-full mx-auto max-w-screen-xl p-4 flex items-center justify-center">
                <span className="text-sm text-primary text-center justify-center dark:text-gray-400">
                    York Structural Biology Laboratory
                    <br />
                    © 2023-2024 University of York. All Rights Reserved.
                    <br />
                    Email: support@privateer.york.ac.uk
                    <br />
                    <a href="https://www.york.ac.uk/about/legal-statements/">
                        Disclaimer
                    </a>
                    <div className="flex space-x-5 flex-wrap justify-center">
                        <img
                            className="w-48 h-auto object-contain mx-auto"
                            src="bbsrc_logo.png"
                        />
                        <img
                            className="w-48 h-auto object-contain mx-auto"
                            src="ccp4_logo.png"
                        />
                        <img
                            className="w-48 h-auto object-contain mx-auto"
                            src="uoy_logo.png"
                        />
                        <img
                            className="w-48 h-auto object-contain mx-auto"
                            src="ysbl_logo.png"
                        />
                        <img
                            className="w-36 h-auto object-contain mx-auto"
                            src="royal_society_logo.png"
                        />
                    </div>
                </span>
            </div>
        </footer>
    );
}
