import React, { useEffect, useState } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate, useLocation } from 'react-router-dom';
import { Session } from '@supabase/supabase-js';
import { supabase } from './lib/supabase';
import Navbar from './components/Navbar';
import Footer from './components/Footer';
import Home from './pages/Home';
import Login from './pages/Login';
import AIChat from './pages/AIChat';
import Profile from './pages/Profile';
import Signup from './pages/Signup';
import Test from './pages/test';

function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    // Check existing session
    supabase.auth.getSession().then(({ data: { session } }) => {
      setSession(session);
      setLoading(false);
    });

    // Listen for auth changes
    const { data: { subscription } } = supabase.auth.onAuthStateChange((event, session) => {
      console.log(`Supabase auth event: ${event}`);
      setSession(session);
    });

    return () => subscription.unsubscribe();
  }, []);

  // Wrapper component to handle auth redirects
  const AuthWrapper = ({ children }: { children: JSX.Element }) => {
    const location = useLocation();
    
    if (loading) {
      return (
        <div className="min-h-screen flex items-center justify-center bg-gray-50">
          <div className="flex flex-col items-center gap-4">
            <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-emerald-600"></div>
            <p className="text-gray-600">Loading...</p>
          </div>
        </div>
      );
    }

    if (!session) {
      return <Navigate to="/login" state={{ from: location }} replace />;
    }

    return children;
  };

  if (loading) {
    return (
      <div className="min-h-screen flex items-center justify-center bg-gray-50">
        <div className="flex flex-col items-center gap-4">
          <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-emerald-600"></div>
          <p className="text-gray-600">Loading...</p>
        </div>
      </div>
    );
  }

  return (
    <Router>
      <div className="flex flex-col min-h-screen">
        <Navbar session={session} />
        <main className="flex-grow">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route 
              path="/login" 
              element={session ? <Navigate to="/" replace /> : <Login />} 
            />
            <Route 
              path="/signup" 
              element={session ? <Navigate to="/" replace /> : <Signup />} 
            />
            <Route 
              path="/profile" 
              element={<AuthWrapper><Profile /></AuthWrapper>} 
            />
            <Route 
              path="/ai-chat" 
              //element={<AuthWrapper><AIChat /></AuthWrapper>} 
              element={<Test />} 
            />
            <Route path="*" element={<Navigate to="/" replace />} />
          </Routes>
        </main>
        <Footer />
      </div>
    </Router>
  );
}

export default App;