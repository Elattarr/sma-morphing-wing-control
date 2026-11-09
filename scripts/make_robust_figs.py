"""Charts and an animated page for the heat-only robustness results.
python scripts/make_robust_figs.py   -> docs/figures/heatonly_robust_bars.png, heatonly_windows.png, docs/heatonly_animation.html  (simulation only)"""
import json, os, re, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import evaluate_heatonly as eh, evaluate_controller as ev
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]

# 1. bars from the committed results table
rows = []
for ln in (ROOT / "docs/results/eval_heatonly_robust.txt").read_text().splitlines()[2:]:
    m = re.match(r"(.+?)\s+(\d) of 4\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s+(\d+)\s*$", ln)
    if m:
        rows.append((m.group(1).strip(), int(m.group(2)), float(m.group(3)), float(m.group(5)), float(m.group(7))))
names = [r[0].replace(" (severe, NOT the data)", " (severe)") for r in rows][::-1]
f, ax = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
col = ["#2a7" if r[1] == 4 else "#c33" for r in rows][::-1]
for a, i, t in zip(ax, (2, 3, 4), ("max hold error [% of stroke]", "peak wire temperature [C]", "energy [J]")):
    a.barh(names, [r[i] for r in rows][::-1], color=col); a.set_title(t); a.grid(axis="x", alpha=.3)
ax[1].axvline(73, color="k", ls="--", lw=1); ax[1].text(73.3, 0, "73 C limit", fontsize=8)
f.suptitle("Heat-only controller, parameter errors (simulation). Green = 4 of 4 targets, red = fewer", y=.99)
f.tight_layout(); f.savefig(ROOT / "docs/figures/heatonly_robust_bars.png", dpi=130); plt.close(f)

# 2. time histories, wide windows
cases = [("nominal (A_s 51, A_f 54)", {}), ("18 K window", dict(p_true=eh.wing_window())), ("30 K window", dict(p_true=eh.wing_window(W=30.0))),
         ("25 K, A_f 73", dict(p_true=eh.wing_window(W=25.0, Af_target=73.0)))]
data = []
f, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
for name, kw in cases:
    L, seg, home = ev.sequence(targets=eh.UP, hold_s=10.0, t_max=1500.0, label=name, **kw)
    L = np.asarray(L); stride = max(1, len(L) // 400)
    t, T, r, tg, I = L[:, 0], L[:, 1], L[:, 2], L[:, 3], L[:, 4]
    ax[0].plot(t, T, label=name); ax[1].plot(t, r * 100); ax[2].plot(t, I)
    data.append(dict(name=name, t=[round(x, 2) for x in t[::stride]], T=[round(x, 2) for x in T[::stride]], r=[round(x, 4) for x in r[::stride]],
                     tg=[round(x, 3) for x in tg[::stride]], I=[round(x, 2) for x in I[::stride]]))
ax[1].step(L[:, 0], L[:, 3] * 100, "k--", lw=1, where="post", label="target")
ax[0].set_ylabel("wire T [C]"); ax[1].set_ylabel("stroke [%]"); ax[2].set_ylabel("current [A]"); ax[2].set_xlabel("time [s]")
ax[0].axhline(73, color="k", ls=":", lw=1); ax[0].legend(fontsize=8); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
f.suptitle("Heat-only controller on four transformation windows (simulation, controller knows only the nominal wire)")
f.tight_layout(); f.savefig(ROOT / "docs/figures/heatonly_windows.png", dpi=130); plt.close(f)

# 3. animated page
html = """<!doctype html><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1"><title>Heat-only controller, animation</title>
<style>body{font:15px system-ui,sans-serif;margin:16px;background:#101418;color:#e6e6e6}select,button{font:inherit;padding:4px 8px}canvas{width:100%;max-width:900px;background:#181e24;border-radius:6px;display:block;margin:8px 0}
.n{color:#aab;max-width:900px}</style>
<h2>Heat-only SMA controller: one run, animated</h2>
<p class=n>Simulation only, nothing measured. The controller knows the nominal wire; the plant has the chosen transformation window. Top: wire temperature. Middle: stroke against target. Bottom: current.</p>
<select id=c></select> <button id=p>pause</button> <input type=range id=s min=0 max=1000 value=0 style="width:300px">
<canvas id=a width=900 height=520></canvas><div id=o class=n></div>
<script>const D=__DATA__;let k=0,run=true,cs=0;const sel=document.getElementById('c'),cv=document.getElementById('a'),g=cv.getContext('2d'),sl=document.getElementById('s'),o=document.getElementById('o');
D.forEach((d,i)=>sel.add(new Option(d.name,i)));sel.onchange=()=>{cs=+sel.value;k=0};document.getElementById('p').onclick=e=>{run=!run;e.target.textContent=run?'pause':'play'};
sl.oninput=()=>{k=Math.floor(sl.value/1000*(D[cs].t.length-1))};
function panel(y0,h,ser,lo,hi,col,lab,tg){g.strokeStyle='#345';g.strokeRect(50,y0,820,h);g.fillStyle='#9ab';g.font='12px sans-serif';g.fillText(lab,54,y0+13);g.fillText(hi,8,y0+12);g.fillText(lo,8,y0+h);
const d=D[cs],n=d.t.length,T=d.t[n-1];const X=i=>50+820*d.t[i]/T,Y=v=>y0+h-(v-lo)/(hi-lo)*h;
if(tg){g.strokeStyle='#888';g.setLineDash([4,3]);g.beginPath();for(let i=0;i<=k;i++){i?g.lineTo(X(i),Y(tg[i])):g.moveTo(X(i),Y(tg[i]))}g.stroke();g.setLineDash([])}
g.strokeStyle=col;g.lineWidth=2;g.beginPath();for(let i=0;i<=k;i++){i?g.lineTo(X(i),Y(ser[i])):g.moveTo(X(i),Y(ser[i]))}g.stroke();g.lineWidth=1;g.fillStyle=col;g.beginPath();g.arc(X(k),Y(ser[k]),4,0,7);g.fill()}
function draw(){const d=D[cs];g.clearRect(0,0,900,520);panel(10,150,d.T,20,80,'#f84','wire temperature [C]');
g.strokeStyle='#fff3';g.beginPath();g.moveTo(50,10+150-(73-20)/60*150);g.lineTo(870,10+150-(73-20)/60*150);g.stroke();
panel(180,150,d.r.map(v=>v*100),0,105,'#4af','stroke [%]',d.tg.map(v=>v*100));panel(350,150,d.I,0,Math.max(...d.I)*1.1||1,'#8d5','current [A]');
o.textContent=`t = ${d.t[k].toFixed(0)} s   T = ${d.T[k].toFixed(1)} C   stroke = ${(d.r[k]*100).toFixed(1)} %   target = ${(d.tg[k]*100).toFixed(0)} %   I = ${d.I[k].toFixed(1)} A`;sl.value=k/(d.t.length-1)*1000}
setInterval(()=>{if(run){const n=D[cs].t.length;k=(k+1)%n}draw()},40);draw();</script>"""
(ROOT / "docs/heatonly_animation.html").write_text(html.replace("__DATA__", json.dumps(data)))
print("ok", [len(d["t"]) for d in data])
