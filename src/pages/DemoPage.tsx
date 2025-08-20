import React, { useEffect, useRef, useState } from "react";

type Msg = { id: string; role: "user" | "assistant"; text?: string; imageUrl?: string | null };

// Use environment variables
const API_URL = import.meta.env.VITE_CHEMCHAT_URL || "https://prithvirao-orgomap-backend.hf.space/api/chemchat";
const API_ORIGIN = import.meta.env.VITE_CHEMCHAT_ORIGIN || "https://prithvirao-orgomap-backend.hf.space";

export default function ChemChat() {
  const [messages, setMessages] = useState<Msg[]>([]);
  const [input, setInput] = useState("");
  const [busy, setBusy] = useState(false);
  const [allowDiagrams, setAllowDiagrams] = useState(true);
  const endRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => { 
    endRef.current?.scrollIntoView({ behavior: "smooth" }); 
  }, [messages]);

  async function send() {
    const msg = input.trim();
    if (!msg || busy) return;
    
    setBusy(true);
    setMessages((m) => [...m, { id: crypto.randomUUID(), role: "user", text: msg }]);
    setInput("");

    try {
      const res = await fetch(API_URL, {
        method: "POST",
        headers: { 
          "Content-Type": "application/json",
          "Accept": "application/json",
        },
        body: JSON.stringify({ 
          message: msg, 
          allow_diagrams: allowDiagrams 
        }),
      });
      
      if (!res.ok) {
        const errorText = await res.text();
        throw new Error(`HTTP ${res.status}: ${errorText}`);
      }
      
      const data = await res.json();

      // Handle image URL construction
      let imageUrl = null;
      if (data.image_url) {
        imageUrl = data.image_url.startsWith('http') 
          ? data.image_url 
          : `${API_ORIGIN}${data.image_url}`;
      }

      setMessages((m) => [
        ...m,
        {
          id: crypto.randomUUID(),
          role: "assistant",
          text: data.text,
          imageUrl: imageUrl,
        },
      ]);
    } catch (e: any) {
      console.error("API Error:", e);
      setMessages((m) => [...m, { 
        id: crypto.randomUUID(), 
        role: "assistant", 
        text: `Error: ${e.message || 'Failed to get response. Please try again.'}` 
      }]);
    } finally {
      setBusy(false);
    }
  }

  function onKey(e: React.KeyboardEvent<HTMLTextAreaElement>) {
    if (e.key === "Enter" && !e.shiftKey) { 
      e.preventDefault(); 
      send(); 
    }
  }

  const clearChat = () => {
    setMessages([]);
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-emerald-50 text-gray-900 flex flex-col">
      {/* Header */}
      <header className="sticky top-0 border-b border-gray-200 bg-white/80 backdrop-blur-sm z-20 shadow-sm">
        <div className="max-w-4xl mx-auto px-4 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 bg-gradient-to-r from-emerald-600 to-blue-600 rounded-lg flex items-center justify-center">
              <svg className="w-6 h-6 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
              </svg>
            </div>
            <div>
              <h1 className="text-xl font-bold text-gray-800">ChemChat</h1>
              <p className="text-sm text-gray-600">AI-Powered Chemistry Assistant</p>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            <label className="flex items-center gap-2 text-sm bg-gray-100 px-3 py-1 rounded-full">
              <div className="relative">
                <input 
                  type="checkbox" 
                  className="sr-only"
                  checked={allowDiagrams} 
                  onChange={(e) => setAllowDiagrams(e.target.checked)} 
                />
                <div className={`w-10 h-6 rounded-full transition-colors ${
                  allowDiagrams ? 'bg-emerald-600' : 'bg-gray-300'
                }`} />
                <div className={`absolute top-0.5 left-0.5 w-5 h-5 bg-white rounded-full transition-transform ${
                  allowDiagrams ? 'transform translate-x-4' : ''
                }`} />
              </div>
              <span className="text-gray-700">Diagrams</span>
            </label>
            
            {messages.length > 0 && (
              <button
                onClick={clearChat}
                className="text-sm text-gray-500 hover:text-gray-700 px-3 py-1 rounded-full hover:bg-gray-100 transition-colors"
              >
                Clear Chat
              </button>
            )}
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="flex-1 overflow-hidden">
        <div className="max-w-4xl mx-auto px-4 py-6 space-y-6 h-full flex flex-col">
          {messages.length === 0 && (
            <div className="text-center py-12">
              <div className="w-24 h-24 mx-auto mb-6 bg-gradient-to-r from-emerald-500 to-blue-500 rounded-2xl flex items-center justify-center">
                <svg className="w-12 h-12 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z" />
                </svg>
              </div>
              <h2 className="text-2xl font-bold text-gray-800 mb-4">Welcome to ChemChat!</h2>
              <p className="text-gray-600 mb-8">Ask me anything about organic chemistry mechanisms and reactions.</p>
              
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 max-w-2xl mx-auto">
                <div className="bg-white rounded-xl p-4 shadow-sm border border-gray-100 hover:shadow-md transition-shadow">
                  <div className="w-8 h-8 bg-emerald-100 rounded-lg flex items-center justify-center mb-2">
                    <svg className="w-4 h-4 text-emerald-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                    </svg>
                  </div>
                  <h3 className="font-semibold text-gray-800 mb-1">Mechanisms</h3>
                  <p className="text-sm text-gray-600">"Show SN1 vs SN2 mechanism"</p>
                </div>
                
                <div className="bg-white rounded-xl p-4 shadow-sm border border-gray-100 hover:shadow-md transition-shadow">
                  <div className="w-8 h-8 bg-blue-100 rounded-lg flex items-center justify-center mb-2">
                    <svg className="w-4 h-4 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 8h10M7 12h4m1 8l-4-4H5a2 2 0 01-2-2V6a2 2 0 012-2h14a2 2 0 012 2v8a2 2 0 01-2 2h-3l-4 4z" />
                    </svg>
                  </div>
                  <h3 className="font-semibold text-gray-800 mb-1">Reactions</h3>
                  <p className="text-sm text-gray-600">"HBr addition to propene"</p>
                </div>
                
                <div className="bg-white rounded-xl p-4 shadow-sm border border-gray-100 hover:shadow-md transition-shadow">
                  <div className="w-8 h-8 bg-purple-100 rounded-lg flex items-center justify-center mb-2">
                    <svg className="w-4 h-4 text-purple-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
                    </svg>
                  </div>
                  <h3 className="font-semibold text-gray-800 mb-1">Explanations</h3>
                  <p className="text-sm text-gray-600">"Explain Markovnikov's rule"</p>
                </div>
              </div>
            </div>
          )}

          {/* Messages */}
          <div className="space-y-4 flex-1 overflow-y-auto">
            {messages.map((m) => (
              <div key={m.id} className={`flex gap-4 ${m.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                {m.role === 'assistant' && (
                  <div className="w-8 h-8 bg-gradient-to-r from-emerald-500 to-blue-500 rounded-full flex items-center justify-center flex-shrink-0">
                    <svg className="w-4 h-4 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
                    </svg>
                  </div>
                )}
                
                <div className={`max-w-[80%] ${m.role === 'user' ? 'order-first' : ''}`}>
                  <div className={`rounded-2xl px-4 py-3 ${
                    m.role === 'user' 
                      ? 'bg-gradient-to-r from-emerald-600 to-emerald-500 text-white' 
                      : 'bg-white border border-gray-200 shadow-sm'
                  }`}>
                    {m.text && (
                      <p className="whitespace-pre-wrap text-sm leading-relaxed">
                        {m.text}
                      </p>
                    )}
                    {m.imageUrl && (
                      <div className="mt-3">
                        <div className="bg-gray-100 rounded-xl p-3 border border-gray-200">
                          <img 
                            src={m.imageUrl} 
                            alt="Reaction diagram" 
                            className="rounded-lg max-w-full mx-auto" 
                            loading="lazy" 
                          />
                        </div>
                        <p className="text-xs text-gray-500 mt-2 text-center">Reaction diagram</p>
                      </div>
                    )}
                  </div>
                  
                  {m.role === 'user' && (
                    <div className="w-8 h-8 bg-gray-200 rounded-full flex items-center justify-center flex-shrink-0 ml-2">
                      <span className="text-xs font-medium text-gray-700">You</span>
                    </div>
                  )}
                </div>
              </div>
            ))}
          </div>
          
          <div ref={endRef} />
        </div>
      </main>

      {/* Footer with Input */}
      <footer className="border-t border-gray-200 bg-white/80 backdrop-blur-sm">
        <div className="max-w-4xl mx-auto px-4 py-4">
          <div className="bg-white rounded-2xl border border-gray-200 shadow-sm p-3">
            <div className="flex items-end gap-3">
              <textarea
                value={input} 
                onChange={(e) => setInput(e.target.value)} 
                onKeyDown={onKey}
                placeholder="Ask any chemistry question… (Shift+Enter for newline)"
                className="flex-1 resize-none bg-transparent outline-none p-3 text-sm min-h-[60px] max-h-40 border-0 focus:ring-0"
                rows={1}
                disabled={busy}
              />
              <button
                onClick={send} 
                disabled={busy || !input.trim()}
                className="px-6 py-3 text-sm font-medium text-white bg-gradient-to-r from-emerald-600 to-emerald-500 rounded-xl hover:from-emerald-700 hover:to-emerald-600 disabled:opacity-50 disabled:cursor-not-allowed transition-all duration-200 transform hover:scale-105 disabled:transform-none"
              >
                {busy ? (
                  <div className="flex items-center gap-2">
                    <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    Thinking…
                  </div>
                ) : (
                  <div className="flex items-center gap-2">
                    <svg className="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 5l7 7-7 7M5 5l7 7-7 7" />
                    </svg>
                    Send
                  </div>
                )}
              </button>
            </div>
            
            <div className="flex items-center justify-between mt-2 px-1">
              <p className="text-xs text-gray-500">
                Demo only. Verify mechanisms before use in academic or professional settings.
              </p>
              <div className="flex items-center gap-4">
                <span className="text-xs text-gray-500">
                  {messages.length} message{messages.length !== 1 ? 's' : ''}
                </span>
              </div>
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
}