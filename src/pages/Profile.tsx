import React, { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { User, LogOut, Settings, BookOpen, Brain } from 'lucide-react';
import { supabase } from '../lib/supabase';

interface UserProfile {
  email: string;
  created_at: string;
}

const Profile = () => {
  const navigate = useNavigate();
  const [profile, setProfile] = useState<UserProfile | null>(null);
  const [loading, setLoading] = useState(true);
  const [stats, setStats] = useState({
    flashcards: 0,
    studyTime: '0',
    lessonsCompleted: 0,
    streak: 0,
  });

  useEffect(() => {
    fetchProfile();
    fetchStats();
  }, []);

  const fetchProfile = async () => {
    try {
      const { data: { user } } = await supabase.auth.getUser();
      if (user) {
        setProfile({
          email: user.email || '',
          created_at: user.created_at,
        });
      }
    } catch (error) {
      console.error('Error fetching profile:', error);
    } finally {
      setLoading(false);
    }
  };

  const fetchStats = async () => {
    try {
      const { data: { user } } = await supabase.auth.getUser();
      if (user) {
        // Fetch flashcards count
        const { count } = await supabase
          .from('flashcards')
          .select('*', { count: 'exact' })
          .eq('user_id', user.id);

        setStats(prev => ({
          ...prev,
          flashcards: count || 0,
          // Other stats would be fetched from their respective tables
          studyTime: '2.5', // Placeholder
          lessonsCompleted: 5, // Placeholder
          streak: 3, // Placeholder
        }));
      }
    } catch (error) {
      console.error('Error fetching stats:', error);
    }
  };

  const handleSignOut = async () => {
    try {
      await supabase.auth.signOut();
      navigate('/');
    } catch (error) {
      console.error('Error signing out:', error);
    }
  };

  if (loading) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="animate-spin rounded-full h-12 w-12 border-t-2 border-b-2 border-indigo-600"></div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50 py-12">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
        <div className="bg-white rounded-lg shadow-lg overflow-hidden">
          {/* Profile Header */}
          <div className="bg-gradient-to-r from-indigo-600 to-purple-600 px-6 py-8">
            <div className="flex items-center justify-between">
              <div className="flex items-center">
                <div className="bg-white p-3 rounded-full">
                  <User className="h-8 w-8 text-indigo-600" />
                </div>
                <div className="ml-4 text-white">
                  <h1 className="text-2xl font-bold">{profile?.email}</h1>
                  <p className="text-indigo-100">
                    Member since {new Date(profile?.created_at || '').toLocaleDateString()}
                  </p>
                </div>
              </div>
              <div className="flex gap-4">
                <button
                  onClick={() => navigate('/settings')}
                  className="bg-white/10 text-white px-4 py-2 rounded-md hover:bg-white/20 transition-colors flex items-center gap-2"
                >
                  <Settings className="h-5 w-5" />
                  Settings
                </button>
                <button
                  onClick={handleSignOut}
                  className="bg-white text-indigo-600 px-4 py-2 rounded-md hover:bg-indigo-50 transition-colors flex items-center gap-2"
                >
                  <LogOut className="h-5 w-5" />
                  Sign Out
                </button>
              </div>
            </div>
          </div>

          {/* Stats Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6 p-6">
            <div className="bg-gray-50 rounded-lg p-6">
              <div className="flex items-center gap-4">
                <div className="bg-indigo-100 p-3 rounded-full">
                  <Brain className="h-6 w-6 text-indigo-600" />
                </div>
                <div>
                  <p className="text-gray-600">Flashcards Created</p>
                  <h3 className="text-2xl font-bold">{stats.flashcards}</h3>
                </div>
              </div>
            </div>
            <div className="bg-gray-50 rounded-lg p-6">
              <div className="flex items-center gap-4">
                <div className="bg-purple-100 p-3 rounded-full">
                  <BookOpen className="h-6 w-6 text-purple-600" />
                </div>
                <div>
                  <p className="text-gray-600">Study Time</p>
                  <h3 className="text-2xl font-bold">{stats.studyTime}h</h3>
                </div>
              </div>
            </div>
            <div className="bg-gray-50 rounded-lg p-6">
              <div className="flex items-center gap-4">
                <div className="bg-green-100 p-3 rounded-full">
                  <BookOpen className="h-6 w-6 text-green-600" />
                </div>
                <div>
                  <p className="text-gray-600">Lessons Completed</p>
                  <h3 className="text-2xl font-bold">{stats.lessonsCompleted}</h3>
                </div>
              </div>
            </div>
            <div className="bg-gray-50 rounded-lg p-6">
              <div className="flex items-center gap-4">
                <div className="bg-orange-100 p-3 rounded-full">
                  <BookOpen className="h-6 w-6 text-orange-600" />
                </div>
                <div>
                  <p className="text-gray-600">Study Streak</p>
                  <h3 className="text-2xl font-bold">{stats.streak} days</h3>
                </div>
              </div>
            </div>
          </div>

          {/* Recent Activity */}
          <div className="p-6 border-t">
            <h2 className="text-xl font-bold mb-4">Recent Activity</h2>
            <div className="space-y-4">
              {/* Placeholder activities */}
              <div className="flex items-center gap-4 p-4 bg-gray-50 rounded-lg">
                <div className="bg-indigo-100 p-2 rounded-full">
                  <Brain className="h-5 w-5 text-indigo-600" />
                </div>
                <div>
                  <p className="font-medium">Created new flashcard set</p>
                  <p className="text-sm text-gray-600">2 hours ago</p>
                </div>
              </div>
              <div className="flex items-center gap-4 p-4 bg-gray-50 rounded-lg">
                <div className="bg-purple-100 p-2 rounded-full">
                  <BookOpen className="h-5 w-5 text-purple-600" />
                </div>
                <div>
                  <p className="font-medium">Completed Biochemistry lesson</p>
                  <p className="text-sm text-gray-600">Yesterday</p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default Profile;