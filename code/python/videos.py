import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sys, os, subprocess, numpy as np
from render import *

OUT = os.path.join(HERE, 'movies'); os.makedirs(OUT, exist_ok=True)
FPS, SPEED, TMAX = 40, 0.4, 6.0          # 0.4x real time keeps the 11 Hz rattle resolved

META = {
 'rolling_nonpref':      ('Rolling without slipping', 'non-preferred sense, ω₀ = −5 rad/s'),
 'rolling_pref':         ('Rolling without slipping', 'preferred sense, ω₀ = +5 rad/s'),
 'coulomb_nonpref':      ('Coulomb friction, μ = 0.15', 'non-preferred sense, ω₀ = −5 rad/s'),
 'coulomb_spin_nonpref': ('Coulomb + drilling friction, ρ = 0.5 mm', 'non-preferred sense, ω₀ = −5 rad/s'),
 'viscous_nonpref':      ('Viscous friction, c = 500 kg/s', 'non-preferred sense, ω₀ = −5 rad/s'),
}

def load(name):
    alias = {'coulomb_nonpref': 'coulomb_non', 'coulomb_spin_nonpref': 'spin05_non'}
    D = np.load(os.path.join(HERE, 'runs', alias.get(name, name) + '.npz')); return D['t'], D['Y']

def series(Y):
    return {'w1': Y[:,0], 'w2': Y[:,1], 'w3': Y[:,2], 'th': np.degrees(Y[:,4])}

def lim(a, pad=0.12, floor=None):
    lo, hi = float(a.min()), float(a.max())
    if floor is not None: lo, hi = min(lo, -floor), max(hi, floor)
    sp = hi - lo; return lo - pad*sp, hi + pad*sp

def encoder(path, W, H):
    return subprocess.Popen(['ffmpeg','-loglevel','error','-y','-f','rawvideo','-pix_fmt','rgb24',
        '-s',f'{W}x{H}','-r',str(FPS),'-i','-','-c:v','libx264','-preset','slow','-crf','20',
        '-pix_fmt','yuv420p','-movflags','+faststart',path], stdin=subprocess.PIPE)

def frame_times():
    n = int(TMAX/SPEED*FPS) + 1
    return np.linspace(0, TMAX, n)

def single(name, still_at=None):
    t, Y = load(name); S = series(Y); title, sub = META[name]
    W, H, LW = 1280, 720, 780
    cam = Cam(LW*SS, H*SS); st = static_layer(cam)
    rows = [(['w1','w2','w3'], 'ω (rad/s)', lim(np.r_[S['w1'][:int(TMAX*1000)],S['w2'][:int(TMAX*1000)],S['w3'][:int(TMAX*1000)]])),
            (['th'], 'θ (deg)', lim(S['th'][:int(TMAX*1000)], floor=0.05))]
    pan = Panel((W-LW)*SS, H*SS, t, S, rows, TMAX)
    enc = encoder(f'{OUT}/{name}.mp4', W, H)
    stills = {}
    for j, tn in enumerate(frame_times()):
        k = min(int(round(tn*1000)), len(t)-1)
        left = st.copy(); draw_body(left, cam, Y[k,3:6], Y[k,0:3])
        d = ImageDraw.Draw(left)
        text(d, (32, 28), title, 22)
        text(d, (32, 60), sub, 15, C_MID)
        text(d, (32, H-44), f't = {t[k]:4.2f} s', 17)
        text(d, (190, H-44), f'ω₃ = {Y[k,2]:+5.2f} rad/s', 17, COL['w3'])
        text(d, (LW-24, H-44), f'{SPEED:g}× real time', 13, C_LIGHT, anchor='ra')
        full = Image.new('RGB', (W*SS, H*SS), 'white')
        full.paste(left, (0, 0)); full.paste(pan.frame(t[k]), (LW*SS, 0))
        fr = full.resize((W, H), Image.LANCZOS)
        enc.stdin.write(fr.tobytes())
        if still_at is not None and abs(tn - still_at) < 0.5*TMAX/len(frame_times()): fr.save(f'{OUT}/{name}_still.png')
    enc.stdin.close(); enc.wait(); print('done', name, flush=True)

def grid4():
    names = ['rolling_nonpref','coulomb_nonpref','coulomb_spin_nonpref','viscous_nonpref']
    labels = ['Rolling without slipping','Coulomb, μ = 0.15','Coulomb + drilling, ρ = 0.5 mm','Viscous, c = 500 kg/s']
    data = [load(n) for n in names]
    W, H, tw, th = 1280, 720, 640, 360
    cam = Cam(tw*SS, th*SS, f=1700*0.92); st = static_layer(cam)
    enc = encoder(f'{OUT}/comparison_all_models.mp4', W, H)
    for tn in frame_times():
        full = Image.new('RGB', (W*SS, H*SS), 'white')
        for q, ((t, Y), lab) in enumerate(zip(data, labels)):
            k = min(int(round(tn*1000)), len(t)-1)
            tile = st.copy(); draw_body(tile, cam, Y[k,3:6], Y[k,0:3])
            d = ImageDraw.Draw(tile)
            text(d, (22, 18), lab, 18)
            text(d, (22, th-34), f'ω₃ = {Y[k,2]:+5.2f} rad/s', 15, COL['w3'])
            full.paste(tile, ((q%2)*tw*SS, (q//2)*th*SS))
        d = ImageDraw.Draw(full)
        d.line([(tw*SS, 0), (tw*SS, H*SS)], fill=(203,213,225), width=SS)
        d.line([(0, th*SS), (W*SS, th*SS)], fill=(203,213,225), width=SS)
        text(d, (W-18, th-34), f't = {tn:4.2f} s', 15, C_MID, anchor='ra')
        text(d, (W-18, H-34), f'non-preferred sense, {SPEED:g}× real time', 13, C_LIGHT, anchor='ra')
        enc.stdin.write(full.resize((W, H), Image.LANCZOS).tobytes())
    enc.stdin.close(); enc.wait(); print('done grid', flush=True)

def senses():
    names = ['rolling_pref', 'rolling_nonpref']
    labels = ['Preferred sense, ω₀ = +5 rad/s', 'Non-preferred sense, ω₀ = −5 rad/s']
    data = [load(n) for n in names]
    W, H, tw, th3 = 1280, 720, 640, 470
    cam = Cam(tw*SS, th3*SS, f=1700*0.95); st = static_layer(cam)
    pans = []
    for (t, Y) in data:
        S = series(Y)
        pans.append(Panel(tw*SS, (H-th3)*SS, t, S, [(['w3'], 'ω₃ (rad/s)', (-6.2, 6.2))], TMAX))
    enc = encoder(f'{OUT}/rolling_preferred_vs_nonpreferred.mp4', W, H)
    for tn in frame_times():
        full = Image.new('RGB', (W*SS, H*SS), 'white')
        for q, ((t, Y), lab) in enumerate(zip(data, labels)):
            k = min(int(round(tn*1000)), len(t)-1)
            tile = st.copy(); draw_body(tile, cam, Y[k,3:6], Y[k,0:3])
            d = ImageDraw.Draw(tile)
            text(d, (24, 20), lab, 19)
            text(d, (24, 50), 'rolling without slipping', 14, C_MID)
            full.paste(tile, (q*tw*SS, 0)); full.paste(pans[q].frame(t[k]), (q*tw*SS, th3*SS))
        d = ImageDraw.Draw(full)
        d.line([(tw*SS, 0), (tw*SS, H*SS)], fill=(203,213,225), width=SS)
        text(d, (W-20, th3-34), f't = {tn:4.2f} s   ·   {SPEED:g}× real time', 13, C_LIGHT, anchor='ra')
        enc.stdin.write(full.resize((W, H), Image.LANCZOS).tobytes())
    enc.stdin.close(); enc.wait(); print('done senses', flush=True)

if __name__ == '__main__':
    job = sys.argv[1]
    if job == 'grid': grid4()
    elif job == 'senses': senses()
    else:
        stills = {'rolling_nonpref': 2.3, 'coulomb_nonpref': 2.3, 'coulomb_spin_nonpref': 2.3, 'viscous_nonpref': 2.3, 'rolling_pref': 2.0}
        single(job, still_at=stills.get(job))
