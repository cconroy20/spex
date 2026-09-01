const list = await (await fetch('http://localhost:9333/json/list')).json();
const ws = new WebSocket(list.find(t => t.type === 'page').webSocketDebuggerUrl);
await new Promise(r => ws.onopen = r);
let id = 0; const w = new Map();
ws.onmessage = e => { const m = JSON.parse(e.data); if (m.id && w.has(m.id)) { w.get(m.id)(m.result || m); w.delete(m.id); } };
const send = (me, p = {}) => new Promise(r => { const i = ++id; w.set(i, r); ws.send(JSON.stringify({ id: i, method: me, params: p })); });
const ev = async e => (await send('Runtime.evaluate', { expression: e, returnByValue: true, awaitPromise: true })).result?.value;
await send('Page.enable');
await send('Page.navigate', { url: 'http://localhost:8731/?r=240#bb=1' });
await new Promise(r => setTimeout(r, 8000));
console.log(await ev(`(() => {
  const sb=document.querySelector('#sidebar'), sp=document.querySelector('#species');
  const sep=sp.querySelector('.chipsep');
  const first=[...sp.children].indexOf(sep);
  const mols=[...document.querySelectorAll('.chip')].filter(c=>c.classList.contains('mol')).length;
  return JSON.stringify({ sepPresent: !!sep, chipsBefore: first, molChips: mols,
    listH: sp.scrollHeight, railVisible: sb.clientHeight, scrolls: sp.scrollHeight>sb.clientHeight,
    err:(document.querySelector('#err')||{}).textContent||'' });
})()`));
const shot = await send('Page.captureScreenshot', { format:'png' });
(await import('fs')).writeFileSync('sep.png', Buffer.from(shot.data,'base64'));
ws.close(); process.exit(0);
