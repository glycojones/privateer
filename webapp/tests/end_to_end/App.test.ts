const puppeteer = require('puppeteer');
describe('NavBar Tests', () => {
    let browser;
    let page;

    beforeAll(async () => {
        jest.setTimeout(35000);
        browser = await puppeteer.launch({ headless: 'new' });
        page = await browser.newPage();
    });

    it('NavBar Text loading', async () => {
        await page.goto('http://localhost:5173/');
        await page.waitForSelector('#main');
        const text = await page.$eval('#main', (e) => {
            return e.textContent;
        });
        expect(text).toBeTruthy();
    });

    afterAll(() => browser.close());
});
