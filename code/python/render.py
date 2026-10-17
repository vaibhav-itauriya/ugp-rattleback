import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sys, numpy as np
from model import rotmat, geom, a1, a2, a3, b
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from scipy.spatial import ConvexHull
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

SS = 2
FONT = next((f for f in ['/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', '/Library/Fonts/DejaVuSans.ttf', 'C:/Windows/Fonts/DejaVuSans.ttf', 'DejaVuSans.ttf'] if os.path.exists(f)), 'DejaVuSans.ttf')
def font(px): return ImageFont.truetype(FONT, int(px*SS))

C_INK=(31,41,51); C_MID=(100,116,132); C_LIGHT=(148,163,184)
COL = {'w1':(0,114,178),'w2':(213,94,0),'w3':(0,158,115),'th':(31,41,51)}

# ---------------------------------------------------------------- mesh ---
NU, NV = 120, 36
u = np.linspace(0, 2*np.pi, NU, endpoint=False)
v = np.linspace(-np.pi/2, 0, NV)
U, V = np.meshgrid(u, v)
VERT = np.stack([a1*np.cos(V)*np.cos(U), a2*np.cos(V)*np.sin(U), a3*np.sin(V)+b], -1)
FACES = []
for j in range(NV-1):
    for i in range(NU):
        i2 = (i+1) % NU
        FACES.append((j*NU+i, j*NU+i2, (j+1)*NU+i2, (j+1)*NU+i))
FACES = np.array(FACES)
VFLAT = VERT.reshape(-1, 3)
CENT = VFLAT[FACES].mean(1)
NB = np.stack([CENT[:,0]/a1**2, CENT[:,1]/a2**2, (CENT[:,2]-b)/a3**2], -1)
NB /= np.linalg.norm(NB, axis=1, keepdims=True)
RIM = VERT[-1]                                   # v = 0 ring, the flat cap
STRIPE = np.array([[-0.92*a1, 0, b], [0.92*a1, 0, b]])
TIP = np.array([0.80*a1, 0, b])

LIGHT = np.array([-0.35, -0.55, 0.80]); LIGHT /= np.linalg.norm(LIGHT)
UNDER = np.array([170, 104, 56]); CAP = np.array([212, 160, 104])

class Cam:
    def __init__(self, W, H, az=-62, el=13, dist=0.80, f=1700, target=(0, 0, 0.010), cy=0.60):
        az, el = np.radians(az), np.radians(el)
        self.T = np.array(target)
        fwd = -np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])
        self.E = self.T - dist*fwd
        self.r = np.cross(fwd, [0, 0, 1]); self.r /= np.linalg.norm(self.r)
        self.u = np.cross(self.r, fwd); self.fwd = fwd
        self.W, self.H, self.f, self.cy = W, H, f*SS*(W/(760*SS)), cy
    def proj(self, P):
        d = np.atleast_2d(P) - self.E
        z = d @ self.fwd
        x = self.f*(d @ self.r)/z + self.W/2
        y = -self.f*(d @ self.u)/z + self.H*self.cy
        return np.stack([x, y], -1), z

def static_layer(cam):
    W, H = cam.W, cam.H
    g = np.linspace(0, 1, H)[:, None]
    top, bot = np.array([252, 252, 253]), np.array([236, 240, 244])
    bg = (top*(1-g) + bot*g).astype(np.uint8)
    img = Image.fromarray(np.repeat(bg[:, None, :], W, 1).reshape(H, W, 3))
    d = ImageDraw.Draw(img)
    S = 0.9
    corners = np.array([[-S,-S,0],[S,-S,0],[S,S,0],[-S,S,0]])
    p, _ = cam.proj(corners)
    d.polygon([tuple(q) for q in p], fill=(229, 233, 238))
    for x in np.arange(-S, S+1e-9, 0.05):
        for seg in ([[x,-S,0],[x,S,0]], [[-S,x,0],[S,x,0]]):
            q, _ = cam.proj(np.array(seg))
            d.line([tuple(q[0]), tuple(q[1])], fill=(214, 220, 228), width=SS)
    return img

def draw_body(img, cam, ang, w, show_contact=False):
    Q = rotmat(ang)
    _, _, r, _, h = geom(ang, w)
    o = np.array([0, 0, h])
    VW = VFLAT @ Q.T + o
    # soft shadow: project along the light onto the table
    sh = VW - (VW[:, 2:3]/LIGHT[2])*LIGHT
    hull = ConvexHull(sh[:, :2]).vertices
    ps, _ = cam.proj(sh[hull])
    small = Image.new('L', (cam.W//4, cam.H//4), 0)
    ImageDraw.Draw(small).polygon([tuple(q/4) for q in ps], fill=150)
    small = small.filter(ImageFilter.GaussianBlur(3.5*SS/2))
    mask = small.resize((cam.W, cam.H), Image.BILINEAR)
    img.paste(Image.new('RGB', img.size, (120, 128, 140)), (0, 0), mask)

    d = ImageDraw.Draw(img)
    P2, _ = cam.proj(VW)
    NW = NB @ Q.T
    CW = CENT @ Q.T + o
    vis = np.einsum('ij,ij->i', NW, cam.E - CW) > 0
    lam = np.clip(NW @ LIGHT, 0, 1)
    Hh = LIGHT + (cam.E - CW)/np.linalg.norm(cam.E - CW, axis=1, keepdims=True)
    Hh /= np.linalg.norm(Hh, axis=1, keepdims=True)
    spec = np.clip(np.einsum('ij,ij->i', NW, Hh), 0, 1)**40
    vdir = (cam.E - CW)/np.linalg.norm(cam.E - CW, axis=1, keepdims=True)
    fill = np.clip(np.einsum('ij,ij->i', NW, vdir), 0, 1)
    col = np.clip(UNDER*(0.30 + 0.50*lam[:, None] + 0.40*fill[:, None]) + 60*spec[:, None], 0, 255).astype(int)
    for k in np.where(vis)[0]:
        poly = [tuple(P2[i]) for i in FACES[k]]
        c = tuple(col[k]); d.polygon(poly, fill=c, outline=c)
    # flat cap, stripe along e1, and a mark at +e1 so the spin sense is visible
    ncap = Q @ np.array([0, 0, 1.0])
    rimw = RIM @ Q.T + o
    if ncap @ (cam.E - rimw.mean(0)) > 0:
        lc = 0.42 + 0.58*max(0.0, ncap @ LIGHT)
        cc = tuple(np.clip(CAP*lc, 0, 255).astype(int))
        pr, _ = cam.proj(rimw); d.polygon([tuple(q) for q in pr], fill=cc, outline=cc)
        ps2, _ = cam.proj(STRIPE @ Q.T + o + ncap*1e-4)
        d.line([tuple(ps2[0]), tuple(ps2[1])], fill=(122, 74, 38), width=3*SS)
        pt, _ = cam.proj(TIP @ Q.T + o + ncap*1e-4)
        x, y = pt[0]; R = 4.5*SS
        d.ellipse([x-R, y-R, x+R, y+R], fill=(250, 246, 240), outline=(122, 74, 38), width=SS)
    if show_contact:
        pc, _ = cam.proj(Q @ r + o)
        x, y = pc[0]; R = 4*SS
        d.ellipse([x-R-SS, y-R-SS, x+R+SS, y+R+SS], fill=(255, 255, 255))
        d.ellipse([x-R, y-R, x+R, y+R], fill=(211, 63, 73))
    return img

# ----------------------------------------------------- trace panel ---
class Panel:
    """matplotlib draws the static axes once; PIL overlays the elapsed part."""
    def __init__(self, W, H, t, series, rows, tmax):
        self.W, self.H = W, H
        dpi = 100*SS
        fig = plt.figure(figsize=(W/dpi, H/dpi), dpi=dpi, facecolor='white')
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5,
                             'axes.edgecolor': '#94a3b8', 'axes.linewidth': 0.7,
                             'xtick.color': '#64748b', 'ytick.color': '#64748b',
                             'axes.labelcolor': '#334155'})
        n = len(rows)
        top, bot, gap = 0.94, 0.08, 0.10
        hh = (top - bot - gap*(n-1))/n
        self.axes, self.lines = [], []
        sel = t <= tmax
        for k, (keys, ylabel, ylim) in enumerate(rows):
            ax = fig.add_axes([0.17, top - (k+1)*hh - k*gap, 0.78, hh])
            for key in keys:
                ax.plot(t[sel], series[key][sel], color='#d5dce4',
                        lw=1.6 if key in ('w3', 'th') else 0.7)
            ax.set_xlim(0, tmax); ax.set_ylim(*ylim)
            ax.set_ylabel(ylabel)
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(color='#e2e8f0', lw=0.6); ax.set_axisbelow(True)
            if k == n-1: ax.set_xlabel('t (s)')
            else: ax.tick_params(labelbottom=False)
            self.axes.append((ax, keys))
        names = {'w1': 'ω₁', 'w2': 'ω₂', 'w3': 'ω₃', 'th': 'θ'}
        ax0, keys0 = self.axes[0]
        hs = [plt.Line2D([], [], color=np.array(COL[k])/255, lw=2 if k=='w3' else 1.2) for k in keys0]
        ax0.legend(hs, [names[k] for k in keys0], loc='upper left', ncol=len(keys0),
                   frameon=False, fontsize=8.5, handlelength=1.4, columnspacing=1.0)
        fig.canvas.draw()
        self.bg = Image.frombuffer('RGBA', fig.canvas.get_width_height(),
                                   fig.canvas.buffer_rgba()).convert('RGB')
        self.px = {}
        Hpx = self.bg.size[1]
        ts = t[sel]; stride = max(1, len(ts)//2400)
        idx = np.arange(0, len(ts), stride)
        self.ts = ts[idx]
        for ax, keys in self.axes:
            for key in keys:
                xy = ax.transData.transform(np.stack([self.ts, series[key][sel][idx]], -1))
                xy[:, 1] = Hpx - xy[:, 1]
                self.px[key] = xy
            self.yb = ax.transData.transform([[0, ax.get_ylim()[0]], [0, ax.get_ylim()[1]]])
        self.axpix = []
        for ax, keys in self.axes:
            yy = ax.transData.transform([[0, ax.get_ylim()[0]], [0, ax.get_ylim()[1]]])[:, 1]
            self.axpix.append((Hpx - yy[1], Hpx - yy[0], keys))
        plt.close(fig)
    def frame(self, tnow):
        img = self.bg.copy(); d = ImageDraw.Draw(img)
        k = int(np.searchsorted(self.ts, tnow))
        for (ytop, ybot, keys) in self.axpix:
            for key in keys:
                pts = self.px[key][:max(k, 2)].copy()
                pts[:, 1] = np.clip(pts[:, 1], ytop, ybot)
                d.line([tuple(p) for p in pts], fill=COL[key],
                       width=int((2.0 if key in ('w3', 'th') else 1.1)*SS))
            if k > 0:
                x = self.px[keys[0]][min(k, len(self.ts)-1), 0]
                d.line([(x, ytop), (x, ybot)], fill=(148, 163, 184), width=SS)
                for key in keys:
                    p = self.px[key][min(k, len(self.ts)-1)]; R = 3.2*SS
                    d.ellipse([p[0]-R, p[1]-R, p[0]+R, p[1]+R], fill=COL[key], outline=(255,255,255), width=SS)
        return img

def text(d, xy, s, px, fill=C_INK, anchor='la'):
    s = s.replace('-', '−') if ('= ' in s and 'non-' not in s) else s
    d.text((xy[0]*SS, xy[1]*SS), s, font=font(px), fill=fill, anchor=anchor)
