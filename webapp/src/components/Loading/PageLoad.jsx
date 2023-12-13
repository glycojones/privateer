import ClipLoader from "react-spinners/ClipLoader";

export default function PageLoad() {
  return (
    <div className="min-h-screen flex flex-col space-y-6 text-center w-screen h-screen bg-gray justify-center items-center text-primary text-m ">
      {/* <h1>Loading Privateer Web App</h1> */}
      <ClipLoader
        loading={true}
        size={100}
        aria-label="Loading Spinner"
        data-testid="loader"
      />
    </div>
  );
}
