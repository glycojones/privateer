const puppeteer = require('puppeteer');
describe('Header Tests', () => {
    let browser;
    let page;

    beforeAll(async () => {
        jest.setTimeout(35000);
        browser = await puppeteer.launch({ headless: 'new' });
        page = await browser.newPage();
        page.on('console', (msg) => {
            console.log('PAGE LOG:', msg.text());
        });
    });

    it('File Upload Present?', async () => {
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#upload');
        const text = await page.$eval('#upload', (e) => {
            return e.textContent;
        });
        expect(text).toContain('Choose a file');
        expect(text).toContain('PDB, mmCIF or MTZ');
        expect(text).toContain('Files will never be sent externally');
    });

    it('PDB Fetch Present?', async () => {
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#PDBFetch');
        const text = await page.$eval('#PDBFetch', (e) => {
            return e.textContent;
        });
        expect(text).toContain('Fetch from PDB');
        expect(text).toContain('Fetch');

        const searchBar = await page.$('#PDBFetch input');
        expect(searchBar).toBeTruthy();

        const fetchButton = await page.$('#PDBFetch button');
        expect(fetchButton).toBeTruthy();
    });

    it('PDB Fetch Works?', async () => {
        await page.setRequestInterception(true);
        page.on('request', (interceptedRequest) => {
            if (interceptedRequest.url().endsWith('.pdb')) {
                interceptedRequest.abort();
            } else interceptedRequest.continue();
        });
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#PDBFetch');
        await page.type('#PDBFetch #code', '2h6o');
        await page.click('#PDBFetch #fetch');

        // after PDBFetch click, upload and PDB should disappear
        // and be replaced with loading
        const uploadBox = await page.$('#upload');
        expect(uploadBox).toBeFalsy();

        const PDBFetch = await page.$('#PDBFetch');
        expect(PDBFetch).toBeFalsy();
    });

    afterAll(() => browser.close());
});
