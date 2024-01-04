module.exports = {
    launch: {
        headless: true,
        slowMo: 50,
        devtools: true,
        product: 'chrome',
        args: [
            `--no-sandbox`,
            `--disable-setuid-sandbox`,
            '--window-size=1920,1080',
        ],
    },
    server: {
        command: 'npx vite --port 5173 --mode development',
        launchTimeout: 10000,
    },
};
