/** @type {import('tailwindcss').Config} */

const colors = require("tailwindcss/colors");

export default {
  content: ["./index.html", "./src/**/*.{js,ts,jsx,tsx}"],
  theme: {
    colors: {
      ...colors,
      gray: "#D6D9E5",
      primary: "#12161E",
      secondary: "#4a5568",
      tertiary: "#F4F9FF",
      hover: "#e0e7ff",
      darkgray: "#B6BAC6",
    },
    extend: {
      fontFamily: {
        body: ["OpenSans", "sans-serif "],
        title: ["Poppins", "cursive"],
      },
    },
  },
  plugins: [],
};
