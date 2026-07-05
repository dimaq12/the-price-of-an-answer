#!/usr/bin/env python3
# Render assets/*.html cards to PNG @2x via headless Chrome.
# NOTE: headless Chrome paints ~90px less than --window-size vertically,
# so we render with headroom and crop to the card's declared height.
import re, subprocess, sys, pathlib
from PIL import Image, ImageChops

ASSETS = pathlib.Path(__file__).resolve().parent

def render(name):
    p = ASSETS / f"{name}.html"
    html = p.read_text(encoding="utf-8")
    h = int(re.search(r"class='card' style='(?:width:\d+px;)?height:(\d+)px'", html).group(1))
    png = str(ASSETS / f"{name}.png")
    subprocess.run(
        ["google-chrome", "--headless=new", "--disable-gpu", "--hide-scrollbars",
         "--force-device-scale-factor=2", f"--window-size=1200,{h+160}",
         f"--screenshot={png}", p.as_uri()],
        check=True, capture_output=True)
    im = Image.open(png).crop((0, 0, 2400, h * 2))
    # trim surplus bottom whitespace to a uniform 96px (48 css) margin
    bbox = ImageChops.difference(im.convert("RGB"),
                                 Image.new("RGB", im.size, "#FFFFFF")).getbbox()
    bottom = min(h * 2, (bbox[3] if bbox else h * 2) + 96)
    im.crop((0, 0, 2400, bottom)).save(png)
    print(f"{name}.png  2400x{bottom}")

for name in sys.argv[1:]:
    render(name)
