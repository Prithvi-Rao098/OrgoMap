import React from 'react';

import { Link } from 'react-router-dom';
import { FlaskConical, BookOpen, Atom, Lightbulb, TrendingUp, ScanEye, User, MessageSquare, Info} from 'lucide-react';
import chemBg from "../imgs/chem-bg.png";

// Place your image in: public/images/hero-bg.jpg
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
          <div className="flex flex-col sm:flex-row justify-center gap-4">
            <Link to="/demo" className="bg-emerald-600 hover:bg-emerald-700 text-white px-8 py-3 rounded-md font-semibold transition-colors">
              Try It Out
            </Link>
            <Link to="/demo" className="bg-transparent border-2 border-white text-white px-8 py-3 rounded-md font-semibold hover:bg-white/10 transition-colors">
              Live Demo
            </Link>
          </div>
        </div>
      </section>


      {/* Features Section */}
      <section className="py-20 bg-white">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <h2 className="text-3xl font-bold text-center mb-12">How It Works</h2>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
            <FeatureCard
              icon={<ScanEye className="h-8 w-8" />}
              title="Upload Reaction"
              description="Paste SMILES or draw your organic reaction"
            />
            <FeatureCard
              icon={<FlaskConical className="h-8 w-8" />}
              title="AI Analysis"
              description="Our engine maps electron flow and intermediates"
            />
            <FeatureCard
              icon={<BookOpen className="h-8 w-8" />}
              title="Learn Mechanisms"
              description="Step-by-step explanations with 3D visualizations"
            />
          </div>
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

// Reusable Feature Card Component
const FeatureCard = ({ icon, title, description }) => (
  <div className="p-6 bg-white rounded-lg shadow-md border border-gray-100 hover:shadow-lg transition-shadow">
    <div className="text-emerald-600 mb-4">{icon}</div>
    <h3 className="text-xl font-semibold mb-2">{title}</h3>
    <p className="text-gray-600">{description}</p>
  </div>
);

export default Home;