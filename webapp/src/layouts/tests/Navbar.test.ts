const puppeteer = require('puppeteer');
describe('NavBar Tests', () => {
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

    it('NavBar Text Loading', async () => {
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#navbar');
        const text = await page.$eval('#navbar', (e) => {
            return e.textContent;
        });
        expect(text).toContain('Validate your carbohydrates online with');
        expect(text).toContain('Privateer');
        expect(text).toContain(
            'The Swiss Army knife for carbohydrate structure validation, refinement and analysis'
        );
    });

    it('Buttons Present', async () => {
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#navbar');
        const ids = await page.$$eval('#navbar a', (e) => {
            return e.map((a) => a.id);
        });
        expect(ids).toContain('database');
        expect(ids).toContain('mainpage');
        expect(ids).toContain('citation');
        expect(ids).toContain('github');
    });

    afterAll(() => browser.close());
});
