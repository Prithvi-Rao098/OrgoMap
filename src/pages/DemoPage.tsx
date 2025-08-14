import React, { useEffect, useRef, useState } from "react";

type Msg = { id: string; role: "user" | "assistant"; text?: string; imageUrl?: string | null };

const API = import.meta.env.VITE_CHEMCHAT_URL || "http://127.0.0.1:8000/api/chemchat";
const API_ORIGIN = import.meta.env.VITE_CHEMCHAT_ORIGIN || "http://127.0.0.1:8000";

export default function ChemChat() {
  const [messages, setMessages] = useState<Msg[]>([]);
  const [input, setInput] = useState("");
  const [busy, setBusy] = useState(false);
  const [allowDiagrams, setAllowDiagrams] = useState(true);
  const endRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => { endRef.current?.scrollIntoView({ behavior: "smooth" }); }, [messages]);

  async function send() {
    const msg = input.trim();
    if (!msg || busy) return;
    setBusy(true);
    setMessages((m) => [...m, { id: crypto.randomUUID(), role: "user", text: msg }]);
    setInput("");

    try {
      const res = await fetch(API, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ message: msg, allow_diagrams: allowDiagrams }),
      });
      if (!res.ok) throw new Error(await res.text());
      const data = await res.json();

      setMessages((m) => [
        ...m,
        {
          id: crypto.randomUUID(),
          role: "assistant",
          text: data.text,
          imageUrl: data.image_url ? `${API_ORIGIN}${data.image_url}` : null,
        },
      ]);
    } catch (e: any) {
      setMessages((m) => [...m, { id: crypto.randomUUID(), role: "assistant", text: `Error: ${e.message || e}` }]);
    } finally {
      setBusy(false);
    }
  }

  function onKey(e: React.KeyboardEvent<HTMLTextAreaElement>) {
    if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); send(); }
  }

  return (
    <div className="min-h-screen bg-neutral-50 text-neutral-900 flex flex-col">
      <header className="sticky top-0 border-b bg-white z-10">
        <div className="max-w-3xl mx-auto px-4 py-3 flex items-center justify-between">
          <h1 className="text-lg font-semibold">Chem Chat + Mechanisms</h1>
          <label className="flex items-center gap-2 text-sm">
            <input type="checkbox" className="h-4 w-4"
              checked={allowDiagrams} onChange={(e) => setAllowDiagrams(e.target.checked)} />
            Enable diagrams
          </label>
        </div>
      </header>

      <main className="flex-1">
        <div className="max-w-3xl mx-auto px-4 py-6 space-y-4">
          {messages.length === 0 && (
            <div className="rounded-2xl border bg-white p-4 text-sm text-neutral-600">
              <p className="font-medium">Try:</p>
              <ul className="list-disc ml-5 mt-1 space-y-1">
                <li>“Explain SN1 vs SN2 briefly”</li>
                <li>“Show the mechanism for HBr addition to propene”</li>
                <li>“Nitration of bromobenzene mechanism”</li>
              </ul>
            </div>
          )}
          {messages.map((m) => (
            <div key={m.id} className="flex gap-3">
              <div className={`h-8 w-8 rounded-full flex items-center justify-center text-xs font-bold
                ${m.role === "user" ? "bg-emerald-600 text-white" : "bg-gray-200 text-gray-800"}`}>
                {m.role === "user" ? "U" : "A"}
              </div>
              <div className="flex-1 space-y-3">
                {m.text && <p className="whitespace-pre-wrap text-sm">{m.text}</p>}
                {m.imageUrl && (
                  <a href={m.imageUrl} target="_blank" rel="noreferrer">
                    <img src={m.imageUrl} alt="Reaction diagram" className="rounded-xl border max-w-full" loading="lazy" />
                  </a>
                )}
              </div>
            </div>
          ))}
          <div ref={endRef} />
        </div>
      </main>

      <footer className="border-t bg-white">
        <div className="max-w-3xl mx-auto px-4 py-3">
          <div className="rounded-2xl border bg-white p-2 flex items-end gap-2">
            <textarea
              value={input} onChange={(e) => setInput(e.target.value)} onKeyDown={onKey}
              placeholder="Ask any chem question… (Shift+Enter for newline)"
              className="flex-1 resize-none bg-transparent outline-none p-3 text-sm min-h-[52px] max-h-40"
            />
            <button
              onClick={send} disabled={busy || !input.trim()}
              className="rounded-xl px-4 py-2 text-sm font-medium text-white bg-emerald-600 disabled:opacity-40">
              {busy ? "Thinking…" : "Send"}
            </button>
          </div>
          <p className="mt-2 text-[11px] text-neutral-500">
            Demo only. Verify mechanisms before use.
          </p>
        </div>
      </footer>
    </div>
  );
}
