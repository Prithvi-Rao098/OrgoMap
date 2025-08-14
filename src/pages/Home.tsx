import React from 'react';
import { Link } from 'react-router-dom';
import chemBg from "../imgs/chem-bg.png";

const Home = () => {
  return (
    <div className="flex flex-col">

      {/* Full-screen Hero with Background Image */}
      <section className="relative h-screen w-full flex items-center justify-center overflow-hidden">
        <img
          className="absolute inset-0 w-full h-full object-cover"
          src={chemBg}
          alt="Chemistry background"
        />
        <div className="absolute inset-0 bg-black/50 z-0" />
        <div className="z-10 text-center px-4 text-white">
          <h1 className="text-5xl md:text-7xl font-bold mb-6">OrgoMap</h1>
          <p className="text-2xl md:text-3xl text-emerald-200 mb-8">
            AI-Powered Organic Chemistry Pathways
          </p>
          <div className="flex justify-center">
            <Link
              to="/demo"
              className="bg-emerald-600 hover:bg-emerald-700 text-white px-8 py-3 rounded-md font-semibold transition-colors"
            >
              Try It Out
            </Link>
          </div>
        </div>
      </section>

      {/* How It Works Section */}
      <section className="py-20 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl font-bold mb-6">How It Works</h2>
          <p className="text-gray-600 max-w-2xl mx-auto">
            Upload your organic reaction as SMILES or by drawing it.
            Our AI analyzes the electron flow and intermediates,
            then shows you the full step-by-step mechanism with visualizations.
          </p>
        </div>
      </section>

      {/* Interactive Demo Section */}
      <section className="py-20 bg-black-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-3xl font-bold mb-8">Try Our Interactive Demo</h2>
          <div className="bg-white p-8 rounded-lg shadow-lg max-w-4xl mx-auto">
            <div className="aspect-w-16 aspect-h-9 bg-gray-100 rounded-lg mb-6 flex items-center justify-center">
              <p className="text-gray-500">[Interactive Reaction Demo Embed]</p>
            </div>
            <Link
              to="/demo"
              className="inline-block bg-emerald-600 text-white px-6 py-3 rounded-md font-semibold hover:bg-emerald-700 transition-colors"
            >
              Launch Full Demo
            </Link>
          </div>
        </div>
      </section>
    </div>
  );
};

export default Home;
