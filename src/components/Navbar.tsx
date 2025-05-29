import React from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { FlaskConical, Menu, X, User, BookOpen, MessageSquare, Atom } from 'lucide-react';
import { supabase } from '../lib/supabase';

interface NavbarProps {
  session: Session | null;
}

const Navbar = ({ session }: NavbarProps) => {
  const [isMenuOpen, setIsMenuOpen] = React.useState(false);
  const navigate = useNavigate();

  const handleSignOut = async () => {
    try {
      await supabase.auth.signOut();
      navigate('/');
    } catch (error) {
      console.error('Error signing out:', error);
    }
  };

  return (
    <nav className="bg-emerald-700 text-white shadow-md">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="flex justify-between h-16">
          <div className="flex items-center">
            <Link to="/" className="flex items-center">
              <FlaskConical className="h-8 w-8 text-emerald-200" />
              <span className="ml-2 text-xl font-bold">OrgoMap</span>
            </Link>
          </div>

          {/* Desktop Navigation */}
          <div className="hidden md:flex md:items-center md:space-x-6">
            <Link to="/lessons" className="flex items-center hover:text-emerald-200 gap-1">
              <BookOpen className="h-5 w-5" /> Lessons
            </Link>
            <Link to="/ai-chat" className="flex items-center hover:text-emerald-200 gap-1">
              <MessageSquare className="h-5 w-5" /> AI Chat
            </Link>
            <Link to="/mechanisms" className="flex items-center hover:text-emerald-200 gap-1">
              <Atom className="h-5 w-5" /> Mechanisms
            </Link>
            
            {session ? (
              <div className="flex items-center space-x-4 ml-4">
                <Link to="/profile" className="hover:text-emerald-200">
                  <User className="h-6 w-6" />
                </Link>
                <button
                  onClick={handleSignOut}
                  className="bg-white/10 hover:bg-white/20 px-4 py-2 rounded-md transition-colors"
                >
                  Sign Out
                </button>
              </div>
            ) : (
              <div className="flex items-center space-x-4 ml-4">
                <Link 
                  to="/login" 
                  className="bg-white/10 hover:bg-white/20 px-4 py-2 rounded-md flex items-center gap-1 transition-colors"
                >
                  <User size={18} /> Login
                </Link>
                <Link 
                  to="/signup" 
                  className="bg-emerald-500 hover:bg-emerald-600 px-4 py-2 rounded-md font-medium transition-colors"
                >
                  Sign Up
                </Link>
              </div>
            )}
          </div>

          {/* Mobile menu button */}
          <div className="flex items-center md:hidden">
            <button
              onClick={() => setIsMenuOpen(!isMenuOpen)}
              className="text-white hover:text-emerald-200"
            >
              {isMenuOpen ? <X className="h-6 w-6" /> : <Menu className="h-6 w-6" />}
            </button>
          </div>
        </div>
      </div>

      {/* Mobile Navigation */}
      {isMenuOpen && (
        <div className="md:hidden bg-emerald-800">
          <div className="px-2 pt-2 pb-3 space-y-1">
            <Link 
              to="/lessons" 
              className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
              onClick={() => setIsMenuOpen(false)}
            >
              <BookOpen className="h-5 w-5" /> Lessons
            </Link>
            <Link 
              to="/ai-chat" 
              className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
              onClick={() => setIsMenuOpen(false)}
            >
              <MessageSquare className="h-5 w-5" /> AI Chat
            </Link>
            <Link 
              to="/mechanisms" 
              className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
              onClick={() => setIsMenuOpen(false)}
            >
              <Atom className="h-5 w-5" /> Mechanisms
            </Link>
            
            {session ? (
              <>
                <Link 
                  to="/profile" 
                  className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
                  onClick={() => setIsMenuOpen(false)}
                >
                  <User className="h-5 w-5" /> Profile
                </Link>
                <button
                  onClick={() => {
                    handleSignOut();
                    setIsMenuOpen(false);
                  }}
                  className="w-full text-left flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
                >
                  Sign Out
                </button>
              </>
            ) : (
              <>
                <Link 
                  to="/login" 
                  className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
                  onClick={() => setIsMenuOpen(false)}
                >
                  <User className="h-5 w-5" /> Login
                </Link>
                <Link 
                  to="/signup" 
                  className="flex items-center gap-2 px-3 py-2 hover:bg-emerald-700 rounded-md"
                  onClick={() => setIsMenuOpen(false)}
                >
                  Sign Up
                </Link>
              </>
            )}
          </div>
        </div>
      )}
    </nav>
  );
};

export default Navbar;