#!/usr/bin/env python3
# Concept illustrations for "The Price of an Answer" (Journey II).
# HTML/SVG -> headless Chrome PNG @2x. Same visual family as Journey I:
# green→amber→red dial, blue = solver/data/price. White cards, big text.
import math, pathlib

OUT = pathlib.Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

INK="#1A2233"; MUTED="#6B7689"; FAINT="#9AA3B2"; HAIR="#EAECEF"
GREEN="#0CA678"; GREEN_BG="#E6FCF5"; RED="#E8453C"; RED_BG="#FFF1F0"
BLUE="#2B6CF6"; BLUE_BG="#EAF0FE"; AMBER="#F59F00"; AMBER_BG="#FFF6E6"
TRACK="#F1F3F5"; INDIGO="#3B4CCA"

CSS=f"""
*{{margin:0;padding:0;box-sizing:border-box;}}
html,body{{font-family:'Lato',sans-serif;background:#fff;color:{INK};
 -webkit-font-smoothing:antialiased;text-rendering:geometricPrecision;}}
.card{{position:relative;width:1200px;padding:48px 60px 40px;}}
.title{{font-size:40px;font-weight:900;letter-spacing:-.01em;}}
"""
def frame(inner,h):
    return (f"<!doctype html><html><head><meta charset='utf-8'><style>{CSS}</style>"
            f"</head><body><div class='card' style='height:{h}px'>{inner}</div>"
            f"</body></html>")

def grad_defs():
    return (f"<defs><linearGradient id='dial' x1='0' y1='0' x2='1' y2='0'>"
            f"<stop offset='0' stop-color='{GREEN}'/>"
            f"<stop offset='.5' stop-color='{AMBER}'/>"
            f"<stop offset='1' stop-color='{RED}'/></linearGradient></defs>")

# ══════════════════════════════════════════════════════════════════════════
# 1 — THE PRICE: solver charges per answer vs 7 solves then free
# ══════════════════════════════════════════════════════════════════════════
def concept_price():
    W,H=1080,480
    # left: the solver, a coin stack per query — expensive forever
    lcx=230
    coins=""
    for i in range(9):
        y=420-i*26
        coins+=(f"<ellipse cx='{lcx}' cy='{y}' rx='60' ry='16' fill='{BLUE_BG}' "
                f"stroke='{BLUE}' stroke-width='2.5'/>")
    left=(f"<text x='{lcx}' y='110' text-anchor='middle' font-size='24' "
          f"font-weight='900' fill='{INK}'>The solver</text>"
          f"<text x='{lcx}' y='140' text-anchor='middle' font-size='19' "
          f"fill='{MUTED}'>charges per answer</text>"
          f"{coins}"
          f"<text x='{lcx}' y='470' text-anchor='middle' font-size='21' "
          f"font-weight='800' fill='{BLUE}'>1000s of solves</text>")
    mid=(f"<text x='540' y='300' text-anchor='middle' font-size='38' "
         f"fill='{FAINT}'>vs</text>")
    rcx=810; rx=660
    dots=""
    for i in range(7):
        x=rx+i*50
        dots+=(f"<circle cx='{x}' cy='320' r='10' fill='{GREEN}'/>")
    flat=(f"<line x1='{rx}' y1='320' x2='{rx+300}' y2='320' stroke='{GREEN}' "
          f"stroke-width='3'/>")
    right=(f"<text x='{rcx}' y='110' text-anchor='middle' font-size='24' "
           f"font-weight='900' fill='{INK}'>lastsolve</text>"
           f"<text x='{rcx}' y='140' text-anchor='middle' font-size='19' "
           f"fill='{MUTED}'>pays once, then answers free</text>"
           f"<text x='{rcx}' y='270' text-anchor='middle' font-size='20' "
           f"font-weight='800' fill='{GREEN}'>~7 solves</text>"
           f"{flat}{dots}"
           f"<text x='{rcx}' y='400' text-anchor='middle' font-size='22' "
           f"font-weight='800' fill='{GREEN}'>then microseconds &mdash; forever</text>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}'>{left}{mid}{right}</svg>")
    inner=(f"<div class='title'>What does one answer actually cost?</div>"
           f"<div style='margin-top:24px'>{svg}</div>")
    return frame(inner,650)

# ══════════════════════════════════════════════════════════════════════════
# 2 — SEVEN SOLVES: Chebyshev nodes → the whole curve at machine precision
# ══════════════════════════════════════════════════════════════════════════
def concept_seven():
    W,H=1080,470
    x0,x1=110,970; y0,ymid=120,300; amp=140
    # the true curve u(T;k) — a smooth wave
    pts=[]
    for i in range(200):
        t=i/199; x=x0+(x1-x0)*t
        y=ymid-amp*math.sin(2.6*t+0.5)*math.exp(-0.4*t)
        pts.append((x,y))
    path="M "+" L ".join(f"{x:.1f} {y:.1f}" for x,y in pts)
    curve=f"<path d='{path}' fill='none' stroke='{BLUE}' stroke-width='3.5'/>"
    # 7 Chebyshev nodes ON the curve
    nodes=""
    for j in range(7):
        t=0.5-0.5*math.cos((2*j+1)*math.pi/14)   # chebyshev in [0,1]
        x=x0+(x1-x0)*t
        y=ymid-amp*math.sin(2.6*t+0.5)*math.exp(-0.4*t)
        nodes+=(f"<line x1='{x:.1f}' y1='{y:.1f}' x2='{x:.1f}' y2='440' "
                f"stroke='{GREEN}' stroke-width='1.5' opacity='.4'/>"
                f"<circle cx='{x:.1f}' cy='{y:.1f}' r='11' fill='#fff' "
                f"stroke='{GREEN}' stroke-width='4'/>")
    axis=(f"<line x1='{x0}' y1='440' x2='{x1}' y2='440' stroke='{HAIR}' "
          f"stroke-width='2'/>"
          f"<text x='{x1}' y='430' text-anchor='end' font-size='18' "
          f"fill='{FAINT}'>parameter k &rarr;</text>"
          f"<text x='{x0-8}' y='150' font-size='19' font-weight='700' "
          f"fill='{BLUE}'>u(T; k)</text>")
    label=(f"<text x='{(x0+x1)/2:.0f}' y='92' text-anchor='middle' font-size='21' "
           f"font-weight='800' fill='{GREEN}'>7 solves &mdash; the green nodes</text>"
           f"<text x='{(x0+x1)/2:.0f}' y='488' text-anchor='middle' font-size='19' "
           f"font-weight='700' fill='{BLUE}'>&hellip;and the blue curve is exact "
           f"between them &mdash; ~5&times;10&#8315;&#185;&#8309; over the whole range</text>")
    svg=f"<svg width='{W}' height='510' viewBox='0 0 {W} 510'>{axis}{curve}{nodes}{label}</svg>"
    inner=(f"<div class='title'>Seven solves buy the whole curve</div>"
           f"<div style='margin-top:22px'>{svg}</div>")
    return frame(inner,600)

# ══════════════════════════════════════════════════════════════════════════
# 3 — THE DIAL reads the family: 35 PDEs, almost all at Φ₁ = 1
# ══════════════════════════════════════════════════════════════════════════
def concept_manifold():
    W,H=1080,470
    bx0,bx1=90,990; by=250; bh=30
    grad=grad_defs()
    bar=f"<rect x='{bx0}' y='{by}' width='{bx1-bx0}' height='{bh}' rx='{bh/2}' fill='url(#dial)'/>"
    def at(phi): return bx0+(bx1-bx0)*min(phi,3.0)/3.0
    # cluster of 33 dots at Φ₁≈1, one outlier at 2.7
    cluster=""
    import random as _r; rr=_r.Random(3)
    for i in range(33):
        x=at(1.0)+rr.uniform(-16,16); y=by-40+rr.uniform(-14,14)
        cluster+=f"<circle cx='{x:.1f}' cy='{y:.1f}' r='6' fill='{GREEN}' opacity='.85'/>"
    soliton=(f"<circle cx='{at(2.7):.1f}' cy='{by-40}' r='10' fill='{RED}'/>"
             f"<text x='{at(2.7):.1f}' y='{by-58}' text-anchor='middle' font-size='19' "
             f"font-weight='800' fill='{RED}'>soliton (2 of 35)</text>")
    ticks=""
    for phi in (1,2):
        x=at(phi)
        ticks+=(f"<line x1='{x:.1f}' y1='{by+bh}' x2='{x:.1f}' y2='{by+bh+10}' "
                f"stroke='{FAINT}' stroke-width='2'/>"
                f"<text x='{x:.1f}' y='{by+bh+34}' text-anchor='middle' font-size='18' "
                f"fill='{MUTED}'>&Phi;&#8321; = {phi}</text>")
    callout=(f"<rect x='{at(1.0)-160:.0f}' y='118' width='320' height='46' rx='23' "
             f"fill='{GREEN_BG}'/>"
             f"<text x='{at(1.0):.0f}' y='148' text-anchor='middle' font-size='21' "
             f"font-weight='800' fill='{GREEN}'>33 of 35 PDEs sit here</text>"
             f"<line x1='{at(1.0):.0f}' y1='164' x2='{at(1.0):.0f}' y2='{by-56}' "
             f"stroke='{GREEN}' stroke-width='1.5' opacity='.5'/>")
    ends=(f"<text x='{bx0}' y='{by+bh+34}' font-size='18' font-weight='700' "
          f"fill='{GREEN}'>one direction &mdash; cheap</text>"
          f"<text x='{bx1}' y='{by+bh+34}' text-anchor='end' font-size='18' "
          f"font-weight='700' fill='{RED}'>a genuine wall</text>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}'>{grad}{callout}{bar}"
         f"{cluster}{soliton}{ticks}{ends}</svg>")
    inner=(f"<div class='title'>The dial reads the whole family at once</div>"
           f"<div style='margin-top:26px'>{svg}</div>")
    return frame(inner,560)

# ══════════════════════════════════════════════════════════════════════════
# 4 — READ BACKWARDS: same kernel, forward and inverse, at the CRB floor
# ══════════════════════════════════════════════════════════════════════════
def concept_backwards():
    W,H=1080,430
    cx=540
    # forward row
    fk=(f"<rect x='120' y='90' width='150' height='70' rx='14' fill='{BLUE_BG}'/>"
        f"<text x='195' y='133' text-anchor='middle' font-size='26' font-weight='900' "
        f"fill='{BLUE}'>k</text>"
        f"<text x='195' y='185' text-anchor='middle' font-size='17' fill='{MUTED}'>parameter</text>")
    fu=(f"<rect x='810' y='90' width='150' height='70' rx='14' fill='{GREEN_BG}'/>"
        f"<text x='885' y='133' text-anchor='middle' font-size='26' font-weight='900' "
        f"fill='{GREEN}'>u</text>"
        f"<text x='885' y='185' text-anchor='middle' font-size='17' fill='{MUTED}'>the field</text>")
    fwd=(f"<line x1='285' y1='125' x2='795' y2='125' stroke='{INK}' stroke-width='3'/>"
         f"<path d='M 795 125 l -16 -9 v 18 z' fill='{INK}'/>"
         f"<text x='540' y='110' text-anchor='middle' font-size='19' font-weight='700' "
         f"fill='{INK}'>forward &mdash; 7 solves</text>")
    # inverse row
    inv=(f"<line x1='795' y1='270' x2='285' y2='270' stroke='{AMBER}' stroke-width='3'/>"
         f"<path d='M 285 270 l 16 -9 v 18 z' fill='{AMBER}'/>"
         f"<text x='540' y='255' text-anchor='middle' font-size='19' font-weight='700' "
         f"fill='{AMBER}'>inverse &mdash; zero extra solves</text>"
         f"<text x='540' y='300' text-anchor='middle' font-size='18' fill='{MUTED}'>"
         f"k&#770; &plusmn; Cram&eacute;r&ndash;Rao bar &mdash; the floor no method can beat</text>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}'>{fk}{fu}{fwd}{inv}</svg>")
    inner=(f"<div class='title'>The same kernel, read backwards</div>"
           f"<div style='margin-top:20px'>{svg}</div>")
    return frame(inner,540)

# ══════════════════════════════════════════════════════════════════════════
# 5 — THE WALL: soliton, the dial saw it first (convergence stalls)
# ══════════════════════════════════════════════════════════════════════════
def concept_wall():
    W,H=1080,500
    L,R,T,B=110,60,70,96
    def X(t): return L+(W-L-R)*t
    def Y(v): return (H-B)-(v)*((H-B)-T)   # v in [0,1], log-ish already normalized
    grid=""
    for gy,lab in [(1.0,"1"),(0.66,"1e-5"),(0.33,"1e-10"),(0.0,"1e-15")]:
        y=Y(gy)
        grid+=(f"<line x1='{L}' y1='{y:.1f}' x2='{W-R}' y2='{y:.1f}' stroke='{HAIR}' "
               f"stroke-width='1'/><text x='{L-12}' y='{y+5:.1f}' text-anchor='end' "
               f"font-size='16' fill='{FAINT}'>{lab}</text>")
    # healthy: plunges to machine precision
    hp=[(X(t),Y(1.0-1.0*min(t*3,1.0))) for t in [0,.12,.24,.36,.5,.7,1.0]]
    hpath="M "+" L ".join(f"{x:.1f} {y:.1f}" for x,y in hp)
    healthy=(f"<path d='{hpath}' fill='none' stroke='{GREEN}' stroke-width='3.5'/>")
    # soliton: stalls near the top (n-width wall)
    sp=[(X(t),Y(1.0-0.30*t)) for t in [0,.15,.3,.5,.7,.85,1.0]]
    spath="M "+" L ".join(f"{x:.1f} {y:.1f}" for x,y in sp)
    sol=(f"<path d='{spath}' fill='none' stroke='{RED}' stroke-width='3.5'/>"
         f"<circle cx='{X(1.0):.1f}' cy='{Y(0.70):.1f}' r='8' fill='{RED}'/>"
         f"<text x='{X(0.60):.0f}' y='{Y(0.92):.0f}' text-anchor='middle' font-size='19' "
         f"font-weight='800' fill='{RED}'>soliton &mdash; &Phi;&#8321; = 2.7, stalls at 7&times;10&#8315;&#179;</text>")
    hlab=(f"<text x='{X(0.60):.0f}' y='{Y(0.055):.0f}' text-anchor='middle' font-size='19' "
          f"font-weight='800' fill='{GREEN}'>healthy PDE &mdash; &Phi;&#8321; = 1.0, machine precision</text>")
    xlab=(f"<text x='{(L+W-R)/2:.0f}' y='{H-26}' text-anchor='middle' font-size='17' "
          f"fill='{MUTED}'>nodes spent &rarr; (7 &hellip; 200 solves)</text>")
    svg=(f"<svg width='{W}' height='{H}' viewBox='0 0 {W} {H}'>{grid}{healthy}{sol}{hlab}{xlab}</svg>")
    inner=(f"<div class='title'>The one wall &mdash; and the dial called it first</div>"
           f"<div style='margin-top:22px'>{svg}</div>")
    return frame(inner,560)

for name,fn in [("concept-price",concept_price),
                ("concept-seven",concept_seven),
                ("concept-manifold",concept_manifold),
                ("concept-backwards",concept_backwards),
                ("concept-wall",concept_wall)]:
    (OUT/f"{name}.html").write_text(fn(),encoding="utf-8")
    print("wrote",name)
