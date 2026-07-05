#!/usr/bin/env python3
# Three hero banners for "The Price of an Answer" (Journey II).
#   v1  "receipt"  — editorial serif title; a solver's runaway bill vs 7 solves
#   v2  "curve"    — 7 nodes light up the whole parametric curve, glowing
#   v3  "ledger"   — the price list: forward / inverse / dial, three glowing rows
import math, pathlib

OUT = pathlib.Path(__file__).resolve().parent
W, H = 1400, 740
GREEN="#19C39A"; AMBER="#FFB020"; RED="#FF5A4D"; BLUE="#4C8DFF"
INK0="#0A0F1E"; INK1="#111A31"; MUT="rgba(255,255,255,.62)"; FAINT="rgba(255,255,255,.40)"

BASE=f"""
*{{margin:0;padding:0;box-sizing:border-box;}}
html,body{{width:{W}px;height:{H}px;overflow:hidden;background:{INK0};
 font-family:'Lato',sans-serif;-webkit-font-smoothing:antialiased;text-rendering:geometricPrecision;}}
.stage{{position:relative;width:{W}px;height:{H}px;overflow:hidden;
 background:linear-gradient(150deg,{INK0} 0%,{INK1} 100%);}}
.dots{{position:absolute;inset:0;opacity:.05;
 background-image:radial-gradient(rgba(255,255,255,.9) 1px, transparent 1px);background-size:28px 28px;}}
.mark{{position:absolute;right:64px;bottom:36px;font-size:14.5px;color:{FAINT};z-index:9;}}
.mark b{{color:rgba(255,255,255,.85);font-weight:800;}}
.kicker{{font-size:14.5px;font-weight:700;letter-spacing:.34em;text-transform:uppercase;color:{FAINT};}}
"""
def page(name,css,body):
    html=(f"<!doctype html><html><head><meta charset='utf-8'><style>{BASE}{css}</style>"
          f"</head><body><div class='stage'>{body}</div></body></html>")
    (OUT/f"{name}.html").write_text(html,encoding="utf-8")
    print("wrote",name)

# ══════════════════════════════════════════════════════════════════════════
# V1 — RECEIPT: serif title; runaway meter vs seven nodes
# ══════════════════════════════════════════════════════════════════════════
def v1():
    # background: a rising "cost" area on the right that a flat green line undercuts
    x0,x1=760,1360; base=560
    area=[]
    for i in range(120):
        t=i/119; x=x0+(x1-x0)*t
        y=base-(30+300*t*t)          # solver cost explodes
        area.append((x,y))
    apath="M "+" L ".join(f"{x:.1f} {y:.1f}" for x,y in area)
    solver=(f"<path d='{apath} L {x1} {base} L {x0} {base} Z' fill='rgba(76,141,255,.16)'/>"
            f"<path d='{apath}' fill='none' stroke='{BLUE}' stroke-width='3'/>"
            f"<text x='{x1}' y='{base-330}' text-anchor='end' font-size='19' "
            f"font-weight='700' fill='{BLUE}'>the solver's bill &mdash; per answer</text>")
    # flat green line low, 7 nodes
    gy=base-14; nodes=""
    for j in range(7):
        gx=x0+ (x1-x0)*(0.06+0.88*j/6)
        nodes+=(f"<circle cx='{gx:.1f}' cy='{gy}' r='9' fill='{GREEN}'/>")
    green=(f"<line x1='{x0}' y1='{gy}' x2='{x1}' y2='{gy}' stroke='{GREEN}' stroke-width='3'/>"
           f"{nodes}<text x='{x0}' y='{gy+34}' font-size='19' font-weight='800' "
           f"fill='{GREEN}'>7 solves &mdash; then microseconds, forever</text>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}' style='position:absolute;inset:0'>"
         f"{solver}{green}</svg>")
    css=f"""
.glow{{position:absolute;left:120px;top:150px;width:520px;height:300px;border-radius:50%;
 background:rgba(25,195,154,.16);filter:blur(100px);}}
.wrap{{position:absolute;left:80px;top:70px;width:640px;z-index:5;}}
.h1{{margin-top:26px;font-family:'EB Garamond',serif;font-size:84px;line-height:1.02;
 color:#F4F6FB;letter-spacing:-.005em;}}
.h1 em{{font-style:italic;background:linear-gradient(95deg,{GREEN},{AMBER} 60%,{RED});
 -webkit-background-clip:text;background-clip:text;color:transparent;}}
.sub{{margin-top:24px;font-size:21px;font-weight:300;line-height:1.55;color:{MUT};max-width:600px;}}
.sub b{{color:#fff;font-weight:700;}}
"""
    body=(f"<div class='glow'></div>{svg}"
          f"<div class='wrap'>"
          f"<div class='kicker'>Spectra&nbsp;Without&nbsp;Matrices &middot; Journey&nbsp;II</div>"
          f"<div class='h1'>The Price of<br><em>an Answer</em></div>"
          f"<div class='sub'>What does one answer actually cost? For 35 nonlinear PDEs, "
          f"the honest price was <b>seven solves</b> &mdash; and microseconds ever after.</div>"
          f"</div>"
          f"<div class='mark'><b>lastsolve</b> &nbsp;&middot;&nbsp; built on resona</div>")
    page("hero-v1",css,body)

# ══════════════════════════════════════════════════════════════════════════
# V2 — CURVE: seven glowing nodes light up the whole parametric curve
# ══════════════════════════════════════════════════════════════════════════
def v2():
    x0,x1=90,1310; ymid=520; amp=150
    def fy(t): return ymid-amp*math.sin(2.4*t+0.6)*math.exp(-0.35*t)
    pts=[(x0+(x1-x0)*i/259, fy(i/259)) for i in range(260)]
    cpath="M "+" L ".join(f"{x:.1f} {y:.1f}" for x,y in pts)
    curve=(f"<path d='{cpath}' fill='none' stroke='url(#cg)' stroke-width='4' "
           f"filter='url(#glow)'/>")
    nodes=""
    for j in range(7):
        t=0.5-0.5*math.cos((2*j+1)*math.pi/14)
        x=x0+(x1-x0)*t; y=fy(t)
        nodes+=(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='13' fill='{GREEN}' "
                f"filter='url(#glow)'/><circle cx='{x:.1f}' cy='{y:.1f}' r='6' fill='#fff'/>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}' style='position:absolute;inset:0'>"
         f"<defs><linearGradient id='cg' x1='0' y1='0' x2='1' y2='0'>"
         f"<stop offset='0' stop-color='{GREEN}'/><stop offset='.55' stop-color='{AMBER}'/>"
         f"<stop offset='1' stop-color='{RED}'/></linearGradient>"
         f"<filter id='glow' x='-40%' y='-40%' width='180%' height='180%'>"
         f"<feDropShadow dx='0' dy='0' stdDeviation='6' flood-color='{GREEN}' flood-opacity='.5'/>"
         f"</filter></defs>{curve}{nodes}</svg>")
    css=f"""
.wrap{{position:absolute;left:80px;top:66px;right:80px;z-index:5;}}
.h1{{margin-top:24px;font-size:70px;font-weight:900;line-height:1.04;letter-spacing:-.02em;color:#fff;}}
.h1 span{{background:linear-gradient(95deg,{GREEN},{AMBER} 55%,{RED});
 -webkit-background-clip:text;background-clip:text;color:transparent;}}
.sub{{margin-top:22px;font-size:22px;font-weight:300;line-height:1.5;color:{MUT};max-width:820px;}}
.sub b{{color:#fff;font-weight:700;}}
.cap{{position:absolute;left:80px;bottom:150px;font-size:18px;color:{FAINT};z-index:6;}}
.cap b{{color:{GREEN};font-weight:800;}}
"""
    body=(f"{svg}"
          f"<div class='wrap'>"
          f"<div class='kicker'>Spectra&nbsp;Without&nbsp;Matrices &middot; Journey&nbsp;II</div>"
          f"<div class='h1'>The Price of <span>an Answer</span></div>"
          f"<div class='sub'>Seven solves light up the entire parametric curve at machine "
          f"precision. The nonlinear solver was <b>redundant</b> at query time all along.</div>"
          f"</div>"
          f"<div class='cap'><b>7 solves</b> &nbsp;&rarr;&nbsp; every answer in between, "
          f"~5&times;10&#8315;&#185;&#8309;, in microseconds</div>"
          f"<div class='mark'><b>lastsolve</b> &nbsp;&middot;&nbsp; built on resona</div>")
    page("hero-v2",css,body)

# ══════════════════════════════════════════════════════════════════════════
# V3 — LEDGER: the price list, three glowing rows
# ══════════════════════════════════════════════════════════════════════════
def v3():
    rows=[("Forward","u(T; k) for any k","~5&times;10&#8315;&#185;&#8309;, 576&times; faster",GREEN),
          ("Inverse","k&#770; from one observation","Cram&eacute;r&ndash;Rao optimal, 0 extra solves",AMBER),
          ("The dial","&Phi;&#8321; &mdash; how much structure","refuses honestly at the wall",RED)]
    ry0=300; rh=112; rw=700; rx=640
    cards=""
    for i,(a,b,c,col) in enumerate(rows):
        y=ry0+i*rh
        cards+=(f"<div style='position:absolute;left:{rx}px;top:{y}px;width:{rw}px;height:88px;"
                f"border-radius:16px;background:rgba(255,255,255,.04);"
                f"border-left:5px solid {col};padding:16px 26px;'>"
                f"<div style='font-size:24px;font-weight:900;color:#fff'>{a}"
                f"<span style='font-size:18px;font-weight:400;color:{MUT};margin-left:12px'>{b}</span></div>"
                f"<div style='font-size:18px;font-weight:700;color:{col};margin-top:8px'>{c}</div></div>")
    css=f"""
.glow1{{position:absolute;left:80px;top:180px;width:420px;height:360px;border-radius:50%;
 background:rgba(25,195,154,.14);filter:blur(110px);}}
.wrap{{position:absolute;left:80px;top:200px;width:480px;z-index:5;}}
.h1{{margin-top:22px;font-size:66px;font-weight:900;line-height:1.05;letter-spacing:-.02em;color:#fff;}}
.h1 span{{background:linear-gradient(95deg,{GREEN},{AMBER} 55%,{RED});
 -webkit-background-clip:text;background-clip:text;color:transparent;}}
.sub{{margin-top:22px;font-size:20px;font-weight:300;line-height:1.5;color:{MUT};}}
.lead{{position:absolute;left:640px;top:236px;font-size:15px;font-weight:800;letter-spacing:.18em;
 text-transform:uppercase;color:{FAINT};z-index:6;}}
"""
    body=(f"<div class='dots'></div><div class='glow1'></div>"
          f"<div class='wrap'>"
          f"<div class='kicker'>Spectra&nbsp;Without&nbsp;Matrices &middot; Journey&nbsp;II</div>"
          f"<div class='h1'>The Price of <span>an Answer</span></div>"
          f"<div class='sub'>Seven solves buy the whole price list &mdash; measured first, "
          f"then paid.</div></div>"
          f"<div class='lead'>what seven solves buy</div>{cards}"
          f"<div class='mark'><b>lastsolve</b> &nbsp;&middot;&nbsp; built on resona</div>")
    page("hero-v3",css,body)

v1(); v2(); v3()
