import { createClient } from '@supabase/supabase-js';

// Use import.meta.env for Vite projects (not process.env)
const supabaseUrl = 'https://fkalnmkdehbpilnabevc.supabase.co';
const supabaseKey = import.meta.env.VITE_SUPABASE_ANON_KEY; // Changed from process.env

if (!supabaseKey) {
  throw new Error('Supabase key is missing. Please check your environment variables.');
}

// Create and export the Supabase client
export const supabase = createClient(supabaseUrl, supabaseKey);