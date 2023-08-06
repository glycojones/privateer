/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    colors: { 
      gray: "#eef2ff",
      primary: '#1b202b',
      secondary: '#4a5568',
      tertiary: '#fafafa',
      hover: "#e0e7ff",

    },
    extend: { 
      fontFamily:{
      body: ['OpenSans', 'sans-serif '],
      title: ['Poppins', 'cursive'], 
    }},
  },
  plugins: [],
}