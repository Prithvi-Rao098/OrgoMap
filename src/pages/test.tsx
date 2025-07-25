import React, { useState } from 'react';
import { Send, MessageSquare, Plus, Bot, FlaskConical, Box } from 'lucide-react';

interface Message {
  id: string;
  text: string;
  isUser: boolean;
  timestamp: Date;
  smileImageUrl?: string;
  sdfUrl?: string;
}

interface Chat {
  id: string;
  title: string;
  lastMessage: Date;
}

const AIChat = () => {
  const [messages, setMessages] = useState<Message[]>([{
    id: '1',
    text: "Hello! I'm your OrgoMap AI Tutor. How can I assist you today?",
    isUser: false,
    timestamp: new Date(),
  }]);
  const [chats, setChats] = useState<Chat[]>([{
    id: '1', title: "First Chat", lastMessage: new Date()
  }]);
  const [activeChat, setActiveChat] = useState('1');
  const [newMessage, setNewMessage] = useState('');
  const [loading, setLoading] = useState(false);
  const [show3DModal, setShow3DModal] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    const userInput = newMessage.trim();
    if (!userInput) return;

    const userMessage: Message = {
      id: Date.now().toString(),
      text: userInput,
      isUser: true,
      timestamp: new Date(),
    };

    setMessages(prev => [...prev, userMessage]);
    setNewMessage('');
    setLoading(true);

    try {
      const response = await fetch("http://localhost:8000/chat", {
        method: "POST",
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ message: userInput })
      });

      if (!response.ok) throw new Error("Failed to get response from server");

      const data = await response.json();
      console.log("Response data:", data);
      const aiResponse: Message = {
        id: (Date.now() + 1).toString(),
        text: data.response || "No response from AI.",
        isUser: false,
        timestamp: new Date(),
        ...(data.image_url ? { smileImageUrl: data.image_url } : {}),
        ...(data.sdf_url ? { sdfUrl: data.sdf_url } : {})
      };
      console.log(aiResponse)
      setMessages(prev => [...prev, aiResponse]);

      if (messages.length === 1) {
        const firstFewWords = userInput.split(' ').slice(0, 5).join(' ');
        setChats(chats.map(chat => chat.id === activeChat ? { ...chat, title: firstFewWords, lastMessage: new Date() } : chat));
      } else {
        setChats(chats.map(chat => chat.id === activeChat ? { ...chat, lastMessage: new Date() } : chat));
      }

    } catch (error) {
      const errorMessage: Message = {
        id: (Date.now() + 1).toString(),
        text: 'An error occurred while processing your request.',
        isUser: false,
        timestamp: new Date(),
      };
      setMessages(prev => [...prev, errorMessage]);
    } finally {
      setLoading(false);
    }
  };

  const startNewChat = () => {
    const newChatId = Date.now().toString();
    setChats([...chats, {
      id: newChatId,
      title: `Chat ${chats.length + 1}`,
      lastMessage: new Date()
    }]);
    setActiveChat(newChatId);
    setMessages([{
      id: '1',
      text: "Hello! I'm your OrgoMap AI Tutor. How can I assist you today?",
      isUser: false,
      timestamp: new Date(),
    }]);
  };

  return (
    <div className="flex h-screen bg-gray-50">
      {/* Sidebar */}
      <div className="w-64 bg-white border-r border-gray-200 flex flex-col">
        <div className="p-4">
          <button onClick={startNewChat} className="flex items-center justify-center w-full gap-2 p-3 text-sm font-medium text-white bg-blue-600 rounded-lg hover:bg-blue-700">
            <Plus className="h-4 w-4" /> New Chat
          </button>
        </div>
        <div className="flex-1 overflow-y-auto">
          <div className="px-2 py-2">
            <h3 className="px-3 mb-1 text-xs font-semibold text-gray-500 uppercase">Recent chats</h3>
            <div className="space-y-1">
              {chats.map(chat => (
                <button key={chat.id} onClick={() => setActiveChat(chat.id)} className={`flex items-center w-full gap-2 p-3 text-sm rounded-lg ${activeChat === chat.id ? 'bg-blue-100 text-blue-800' : 'hover:bg-gray-100 text-gray-700'}`}>
                  <MessageSquare className="h-4 w-4" />
                  <span className="truncate">{chat.title}</span>
                </button>
              ))}
            </div>
          </div>
        </div>
        <div className="p-4 border-t border-gray-200">
          <div className="flex items-center gap-2 text-sm text-gray-500">
            <FlaskConical className="h-4 w-4" /> <span>OrgoMap AI Tutor</span>
          </div>
        </div>
      </div>

      {/* Main Chat Area */}
      <div className="flex-1 flex flex-col overflow-hidden bg-white">
        <div className="flex-1 overflow-y-auto p-4">
          <div className="max-w-3xl mx-auto space-y-6">
            {messages.map((message) => (
              <div key={message.id} className={`flex ${message.isUser ? 'justify-end' : 'justify-start'}`}>
                <div className={`max-w-[80%] rounded-lg p-4 ${message.isUser ? 'bg-blue-600 text-white' : 'bg-gray-100 text-gray-900'}`}>
                  <div className="flex items-center gap-2 mb-1">
                    {!message.isUser && <Bot className="h-4 w-4 text-blue-600" />}
                    <span className="text-xs font-medium">{message.isUser ? 'You' : 'OrgoMap AI'}</span>
                  </div>
                  <p className="whitespace-pre-wrap">{message.text}</p>
                  {message.smileImageUrl && (
                    <>
                      <img
                        //key={message.smileImageUrl}
                        src={`http://localhost:8000${message.smileImageUrl}`}
                        alt="SMILES 2D"
                        className="mt-2 rounded border w-48"
                      />
                      {console.log("smileImageUrl for message", message.id, ":", message.smileImageUrl)}
                    </>
                  )}
                  {message.sdfUrl && (
                    <button
                      className="mt-2 text-sm text-blue-600 underline flex items-center gap-1"
                      onClick={() => setShow3DModal(true)}
                    >
                      <Box className="w-4 h-4" /> View 3D Molecule
                    </button>
                  )}
                  <span className="text-xs opacity-75 mt-1 block">{message.timestamp.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}</span>
                </div>
              </div>
            ))}
            {loading && (
              <div className="flex justify-start">
                <div className="bg-gray-100 text-gray-900 rounded-lg p-4 max-w-[80%]">
                  <div className="flex items-center gap-2">
                    <Bot className="h-4 w-4 text-blue-600" />
                    <span className="text-xs font-medium">OrgoMap AI</span>
                  </div>
                  <div className="flex gap-1 mt-2">
                    <div className="w-2 h-2 rounded-full bg-blue-400 animate-bounce" />
                    <div className="w-2 h-2 rounded-full bg-blue-400 animate-bounce delay-100" />
                    <div className="w-2 h-2 rounded-full bg-blue-400 animate-bounce delay-200" />
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>

        {/* Input */}
        <div className="p-4 border-t border-gray-200 bg-white">
          <div className="max-w-3xl mx-auto">
            <form onSubmit={handleSubmit}>
              <div className="flex gap-2">
                <input
                  type="text"
                  value={newMessage}
                  onChange={(e) => setNewMessage(e.target.value)}
                  placeholder="Ask about reaction mechanisms, synthesis, or molecular structures..."
                  className="flex-grow p-3 border rounded-lg focus:outline-none focus:ring-2 focus:ring-blue-500"
                  disabled={loading}
                />
                <button type="submit" className="bg-blue-600 text-white p-3 rounded-lg hover:bg-blue-700 disabled:bg-blue-400" disabled={loading || !newMessage.trim()}>
                  {loading ? (
                    <div className="w-5 h-5 flex justify-center items-center">
                      <div className="w-3 h-3 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    </div>
                  ) : <Send className="h-5 w-5" />}
                </button>
              </div>
            </form>
            <p className="text-xs text-gray-500 mt-2 text-center">
              OrgoMap AI can make mistakes. Consider checking important information.
            </p>
          </div>
        </div>
      </div>

      {/* 3D Molecule Modal */}
      {show3DModal && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
          <div className="bg-white rounded-lg shadow-lg w-full max-w-4xl p-4 relative">
            <button onClick={() => setShow3DModal(false)} className="absolute top-2 right-3 text-gray-500 hover:text-black">✖</button>
            <iframe
              src="training\static\3d.html"
              title="3D Molecule Viewer"
              width="100%"
              height="500px"
              className="rounded"
            ></iframe>
          </div>
        </div>
      )}
    </div>
  );
};

export default AIChat;
