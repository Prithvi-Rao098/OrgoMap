/*
  # Create flashcards table

  1. New Tables
    - `flashcards`
      - `id` (uuid, primary key)
      - `user_id` (uuid, references auth.users)
      - `question` (text)
      - `answer` (text)
      - `created_at` (timestamp)

  2. Security
    - Enable RLS on `flashcards` table
    - Add policies for users to:
      - Read their own flashcards
      - Create their own flashcards
      - Update their own flashcards
      - Delete their own flashcards
*/

CREATE TABLE IF NOT EXISTS flashcards (
  id uuid PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id uuid REFERENCES auth.users NOT NULL,
  question text NOT NULL,
  answer text NOT NULL,
  created_at timestamptz DEFAULT now()
);

ALTER TABLE flashcards ENABLE ROW LEVEL SECURITY;

-- Policy for users to read their own flashcards
CREATE POLICY "Users can read own flashcards"
  ON flashcards
  FOR SELECT
  TO authenticated
  USING (auth.uid() = user_id);

-- Policy for users to create their own flashcards
CREATE POLICY "Users can create own flashcards"
  ON flashcards
  FOR INSERT
  TO authenticated
  WITH CHECK (auth.uid() = user_id);

-- Policy for users to update their own flashcards
CREATE POLICY "Users can update own flashcards"
  ON flashcards
  FOR UPDATE
  TO authenticated
  USING (auth.uid() = user_id)
  WITH CHECK (auth.uid() = user_id);

-- Policy for users to delete their own flashcards
CREATE POLICY "Users can delete own flashcards"
  ON flashcards
  FOR DELETE
  TO authenticated
  USING (auth.uid() = user_id);