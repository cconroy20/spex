const list = await (await fetch('http://localhost:9333/json/list')).json();
const ws = new WebSocket(list.find(t => t.type === 'page').webSocketDebuggerUrl);
await new Promise(r => ws.onopen = r);
let id=0; const w=new Map(); const errs=[];
ws.onmessage=e=>{const m=JSON.parse(e.data);
  if(m.method==='Runtime.exceptionThrown') errs.push((m.params.exceptionDetails.exception||{}).description);
  if(m.id&&w.has(m.id)){w.get(m.id)(m.result||m);w.delete(m.id);}};
const send=(me,p={})=>new Promise(r=>{const i=++id;w.set(i,r);ws.send(JSON.stringify({id:i,method:me,params:p}));});
const ev=async e=>(await send('Runtime.evaluate',{expression:e,returnByValue:true,awaitPromise:true})).result?.value;
const sleep=ms=>new Promise(r=>setTimeout(r,ms));
await send('Runtime.enable'); await send('Page.enable');
await send('Page.navigate',{url:'http://localhost:8731/?r=610#w=5160.00-5200.00&s=Fe%20I&R=50000&snr=50,3'}); await sleep(9000);
console.log('buttons:', await ev(`[...document.querySelectorAll('.seg.pp')].map(b=>b.dataset.p).join(' ')`));
for (const p of [1,2,3,4,5]) {
  await ev(`[...document.querySelectorAll('.seg.pp')].find(b=>b.dataset.p==='${p}').click()`);
  await sleep(1600);
  console.log(`  p=${p}: ppre=${await ev('state.ppre')}  step=${(await ev('pixelStep(state.R,state.ppre)')).toFixed(2)}  usable=${await ev('pixelGridOK(state.R,state.ppre)')}  label="${await ev(`document.querySelector('#snrlabel').textContent`)}"`);
}
// p=5 at the highest R must be refused
await ev(`(() => { const e=document.querySelector('#res'); e.value=rToPos(300000);
  e.dispatchEvent(new Event('input',{bubbles:true})); })()`); await sleep(2200);
console.log('R=300k,p=5:', await ev(`document.querySelector('#snrlabel').textContent`));
console.log('errors:', errs.length?errs.slice(0,2):'(none)');
ws.close(); process.exit(0);
