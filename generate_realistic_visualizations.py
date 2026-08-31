#!/usr/bin/env python3
"""
ns3-ntn-toolkit — realistic visualisation generator
===================================================

Outputs three artefacts:
  1. ns3_ntn_toolkit_architecture.png   — clean layered architecture
  2. ntn_realistic_mobility.gif         — 14 UEs (7 classes) over a LEO pass
  3. ntn_handover_realistic.gif         — same constellation showing
                                          per-class HO behaviour

Each UE class gets a hand-drawn matplotlib icon (no emoji fonts):
  static handheld -> phone shape
  pedestrian      -> stick figure
  vehicular       -> car silhouette
  HST             -> train capsule
  maritime        -> ship hull + mast
  aviation        -> aircraft outline
  IoT             -> dish antenna

Author: Muhammad Uzair
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import (
    FancyBboxPatch, Circle, Polygon, Rectangle, Wedge, FancyArrowPatch, Arc,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       '..', 'visualization')
os.makedirs(OUT_DIR, exist_ok=True)

# =============================================================================
# Colour palette (picked for high contrast on dark navy background)
# =============================================================================
PALETTE = {
    'static':     '#9e9e9e',
    'pedestrian': '#00e676',
    'vehicular':  '#ff9100',
    'hst':        '#e91e63',
    'maritime':   '#00b8d4',
    'aviation':   '#b388ff',
    'iot':        '#ffc107',
    'sat':        '#00d4ff',
    'beam':       '#1976d2',
    'gnd_link':   '#ffffff',
    'feeder':     '#ffc10780',
    'isl':        '#448aff80',
}

# =============================================================================
# UE icon library — each draws a small symbol at (cx, cy) of size s
# =============================================================================

def _phone(ax, cx, cy, s, c):
    body = FancyBboxPatch((cx-0.45*s, cy-0.7*s), 0.9*s, 1.4*s,
                          boxstyle='round,pad=0.0,rounding_size=0.15',
                          fc=c, ec='black', lw=0.8, zorder=10)
    ax.add_patch(body)
    ax.add_patch(Rectangle((cx-0.30*s, cy-0.45*s), 0.6*s, 0.9*s,
                            fc='white', ec='black', lw=0.4, zorder=11))


def _pedestrian(ax, cx, cy, s, c):
    # head
    ax.add_patch(Circle((cx, cy+0.55*s), 0.20*s, fc=c, ec='black',
                        lw=0.8, zorder=10))
    # body
    ax.plot([cx, cx], [cy-0.35*s, cy+0.35*s], color=c, lw=2.2, zorder=10,
            solid_capstyle='round')
    # arms
    ax.plot([cx-0.30*s, cx+0.30*s], [cy+0.15*s, cy+0.05*s],
            color=c, lw=2.0, zorder=10, solid_capstyle='round')
    # legs
    ax.plot([cx, cx-0.25*s], [cy-0.35*s, cy-0.85*s],
            color=c, lw=2.0, zorder=10, solid_capstyle='round')
    ax.plot([cx, cx+0.25*s], [cy-0.35*s, cy-0.80*s],
            color=c, lw=2.0, zorder=10, solid_capstyle='round')


def _car(ax, cx, cy, s, c):
    body = FancyBboxPatch((cx-0.85*s, cy-0.30*s), 1.7*s, 0.55*s,
                          boxstyle='round,pad=0.0,rounding_size=0.10',
                          fc=c, ec='black', lw=0.8, zorder=10)
    ax.add_patch(body)
    # roof
    pts = [(cx-0.55*s, cy+0.25*s), (cx+0.55*s, cy+0.25*s),
           (cx+0.40*s, cy+0.55*s), (cx-0.40*s, cy+0.55*s)]
    ax.add_patch(Polygon(pts, fc=c, ec='black', lw=0.7, zorder=10))
    # wheels
    ax.add_patch(Circle((cx-0.55*s, cy-0.35*s), 0.16*s,
                        fc='#1c1c1c', ec='black', lw=0.5, zorder=11))
    ax.add_patch(Circle((cx+0.55*s, cy-0.35*s), 0.16*s,
                        fc='#1c1c1c', ec='black', lw=0.5, zorder=11))


def _train(ax, cx, cy, s, c):
    # nose + body capsule
    body = FancyBboxPatch((cx-1.10*s, cy-0.30*s), 2.2*s, 0.55*s,
                          boxstyle='round,pad=0.0,rounding_size=0.30',
                          fc=c, ec='black', lw=0.8, zorder=10)
    ax.add_patch(body)
    # window strip
    ax.add_patch(Rectangle((cx-0.85*s, cy-0.05*s), 1.7*s, 0.20*s,
                            fc='white', ec='black', lw=0.3, zorder=11))
    # rail tick
    ax.plot([cx-1.20*s, cx+1.20*s], [cy-0.45*s, cy-0.45*s],
            color='#37474f', lw=0.8, zorder=9)


def _ship(ax, cx, cy, s, c):
    # hull (trapezoid)
    pts = [(cx-1.0*s, cy-0.05*s), (cx+1.0*s, cy-0.05*s),
           (cx+0.7*s, cy-0.45*s), (cx-0.7*s, cy-0.45*s)]
    ax.add_patch(Polygon(pts, fc=c, ec='black', lw=0.8, zorder=10))
    # superstructure
    ax.add_patch(Rectangle((cx-0.30*s, cy-0.05*s), 0.6*s, 0.40*s,
                            fc='white', ec='black', lw=0.6, zorder=11))
    # mast
    ax.plot([cx, cx], [cy+0.35*s, cy+0.85*s], color='black', lw=1.0, zorder=11)
    # waterline
    ax.plot([cx-1.10*s, cx+1.10*s], [cy-0.50*s, cy-0.50*s],
            color='#0288d1', lw=0.6, zorder=9)


def _plane(ax, cx, cy, s, c):
    # fuselage
    ax.add_patch(FancyBboxPatch((cx-0.85*s, cy-0.10*s), 1.7*s, 0.20*s,
                                 boxstyle='round,pad=0.0,rounding_size=0.10',
                                 fc=c, ec='black', lw=0.6, zorder=10))
    # wings
    pts_wing = [(cx-0.10*s, cy+0.05*s), (cx+0.30*s, cy+0.05*s),
                (cx+0.05*s, cy+0.65*s), (cx-0.30*s, cy+0.65*s)]
    ax.add_patch(Polygon(pts_wing, fc=c, ec='black', lw=0.5, zorder=10))
    pts_wing2 = [(cx-0.10*s, cy-0.05*s), (cx+0.30*s, cy-0.05*s),
                 (cx+0.05*s, cy-0.65*s), (cx-0.30*s, cy-0.65*s)]
    ax.add_patch(Polygon(pts_wing2, fc=c, ec='black', lw=0.5, zorder=10))
    # tail
    pts_tail = [(cx-0.85*s, cy+0.05*s), (cx-0.55*s, cy+0.05*s),
                (cx-0.70*s, cy+0.40*s), (cx-0.85*s, cy+0.30*s)]
    ax.add_patch(Polygon(pts_tail, fc=c, ec='black', lw=0.5, zorder=10))
    # nose dot
    ax.add_patch(Circle((cx+0.85*s, cy), 0.06*s, fc='white',
                        ec='black', lw=0.4, zorder=11))


def _iot(ax, cx, cy, s, c):
    # base box
    ax.add_patch(Rectangle((cx-0.30*s, cy-0.55*s), 0.6*s, 0.45*s,
                            fc=c, ec='black', lw=0.6, zorder=10))
    # antenna mast
    ax.plot([cx, cx], [cy-0.10*s, cy+0.30*s], color='black', lw=1.0, zorder=11)
    # dish arc
    ax.add_patch(Wedge((cx, cy+0.30*s), 0.45*s, 30, 150,
                        fc=c, ec='black', lw=0.6, zorder=11))


ICONS = {
    'static':     _phone,
    'pedestrian': _pedestrian,
    'vehicular':  _car,
    'hst':        _train,
    'maritime':   _ship,
    'aviation':   _plane,
    'iot':        _iot,
}


# =============================================================================
# 1. Toolkit architecture diagram
# =============================================================================

def generate_toolkit_architecture():
    print("[1/3] toolkit architecture diagram...")

    fig, ax = plt.subplots(1, 1, figsize=(22, 15.5), facecolor='white')
    ax.set_xlim(0, 20); ax.set_ylim(0, 13.5)
    ax.axis('off')

    def box(x, y, w, h, title, color, sub='', tc='white', body_fontsize=11.0,
            title_fontsize=14.5):
        r = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.05,rounding_size=0.20',
                           fc=color, ec='#37474f', lw=1.6, alpha=0.92)
        ax.add_patch(r)
        ax.text(x+w/2, y+h-0.36, title, ha='center', va='top',
                fontsize=title_fontsize, fontweight='bold', color=tc)
        if sub:
            ax.text(x+w/2, y+h-0.92, sub, ha='center', va='top',
                    fontsize=body_fontsize, color=tc, alpha=0.95,
                    linespacing=1.35)

    def arrow(x1, y1, x2, y2, label='', color='#37474f', lw=1.8, label_dy=0.22):
        a = FancyArrowPatch((x1, y1), (x2, y2),
                             arrowstyle='-|>', mutation_scale=18,
                             color=color, lw=lw, zorder=4)
        ax.add_patch(a)
        if label:
            mx, my = (x1+x2)/2, (y1+y2)/2 + label_dy
            ax.text(mx, my, label, fontsize=10, color='#263238',
                    ha='center', fontstyle='italic',
                    bbox=dict(boxstyle='round,pad=0.22', fc='white',
                              ec='#cfd8dc', lw=0.7, alpha=0.95))

    # ── Title block ──
    ax.text(10, 12.85,
            'ns3-ntn-toolkit — Integrated 6G Non-Terrestrial Network Simulation Platform',
            fontsize=22, fontweight='bold', ha='center', color='#0d47a1')
    ax.text(10, 12.35,
            'ns-3.43  •  3GPP Rel-17 NTN  •  O-RAN  •  HITRAN/ITU-R THz  •  Realistic UE Mobility (TR 38.811)',
            fontsize=13, ha='center', color='#37474f', fontstyle='italic')

    # ── Tier 4: ns-3 core (full width at bottom) ──
    box(0.5, 0.6, 19.0, 1.2, 'ns-3.43 core',
        '#37474f', tc='white',
        sub='Simulator engine  •  events  •  channels  •  helpers  •  attribute system  •  tracing',
        body_fontsize=12.0, title_fontsize=15)

    # ── Tier 3: upstream contrib ──
    box(0.5, 2.3, 4.5, 1.4, 'satellite (SNS3)', '#01579b',
        sub='SGP4 orbit propagator\n3GPP TR 38.811 NTN channel\nLoo / Markov fading',
        body_fontsize=9.5, title_fontsize=13)
    box(5.3, 2.3, 4.5, 1.4, 'mmwave', '#1565c0',
        sub='5G NR PHY/MAC\nDual connectivity\nBeamforming',
        body_fontsize=9.5, title_fontsize=13)
    box(10.1, 2.3, 4.5, 1.4, 'lte (patched)', '#283593',
        sub='X2 PDCP/RLC\nInter-RAT HO\nMcEnbPdcp / McUePdcp',
        body_fontsize=9.5, title_fontsize=13)
    box(14.9, 2.3, 4.6, 1.4, 'traffic + magister-stats', '#3949ab',
        sub='NTN traffic generators\nStatistics collectors',
        body_fontsize=9.5, title_fontsize=13)

    # ── Tier 2: custom contrib (this work) ──
    box(0.5, 4.4, 4.0, 2.2, 'contrib/ntn-cho', '#bf360c',
        sub='3GPP Rel-17 CHO state machine\nTTE estimator (binary search)\nN(0,σ²) shadow fading\nITU-R P.676 atmosphere\nntn-realistic-mobility helper\n7 UE classes / 5 scenarios',
        body_fontsize=10.0, title_fontsize=14)
    box(4.8, 4.4, 4.0, 2.2, 'contrib/oran-ntn', '#1b5e20',
        sub='Space O-RAN architecture\nE2 / A1 / O1 interfaces\n13 xApps + Near-RT RIC\n28 E2SM-RC actions\n11 A1 policies\n5 conflict strategies',
        body_fontsize=10.0, title_fontsize=14)
    box(9.1, 4.4, 4.0, 2.2, 'contrib/thz-ntn', '#4a148c',
        sub='HITRAN-2020 line-by-line\nITU-R P.835 6-layer atm.\nUM-MIMO array factor\nRIS / ISAC / EKF tracker\n9 example sweeps\nVan Vleck-Weisskopf line shape',
        body_fontsize=10.0, title_fontsize=14)
    box(13.4, 4.4, 4.0, 2.2, 'contrib/ai (ns3-ai fork)', '#311b92',
        sub='Gymnasium API\nShared-memory IPC\nLibTorch / pybind11\nns-3.43 + Py 3.13 + NumPy 2 fixes\n4 RL agents (DQN / Dueling /\nLSTM / FedDQN)',
        body_fontsize=10.0, title_fontsize=14)
    box(17.7, 4.4, 1.8, 2.2, 'contrib/\ntraffic-ntn', '#01579b',
        sub='NTN-aware\nflows', body_fontsize=10, title_fontsize=12.5)

    # ── Tier 1: applications / outputs ──
    box(0.5, 7.3, 4.5, 1.6, 'Reference Scenarios', '#0277bd',
        sub='Walker-Star 6×11=66 sat\n10-seed × 600-s Monte-Carlo\nMulti-xApp coexistence  •  THz LEO pass',
        body_fontsize=10.0, title_fontsize=13.5)
    box(5.3, 7.3, 4.5, 1.6, 'Analysis & Figures', '#00695c',
        sub='Python build_figures.py\n300-DPI PDF exports  •  mc_table.tex\nReproducible anti-overlap layout',
        body_fontsize=10.0, title_fontsize=13.5)
    box(10.1, 7.3, 4.5, 1.6, 'Reinforcement-Learning Agents', '#5e35b1',
        sub='Gymnasium environments\nReward shaping  •  Federated learning\n68-feature observation space',
        body_fontsize=10.0, title_fontsize=13.5)
    box(14.9, 7.3, 4.6, 1.6, '3-D Visualisation Dashboard', '#c62828',
        sub='CesiumJS web viewer\nConstellation animation\nPer-UE HO trace  •  KPM heat map',
        body_fontsize=10.0, title_fontsize=13.5)

    # ── Side legend: standards alignment ──
    box(0.5, 9.3, 19.0, 2.7, '', '#eceff1', tc='#263238',
        sub='')
    standards = [
        ('3GPP TR 38.811 v15.4',  'NTN UE classes & channel'),
        ('3GPP TS 38.331 Rel-17', 'CHO state machine + T304'),
        ('3GPP TR 38.821 v16.2',  'NR-NTN solutions'),
        ('3GPP TR 38.901 v17',    'Channel 0.5–100 GHz / HST'),
        ('ITU-R P.676-13',        'Gas absorption reference'),
        ('ITU-R P.618-14',        'Earth-space scintillation'),
        ('O-RAN.WG3.E2AP / A1AP', 'O-RAN interface profiles'),
        ('IMO Res. A.857(20)',    'Maritime TSS lanes'),
    ]
    ax.text(10, 11.65, 'Standards alignment',
            ha='center', fontsize=14.5, fontweight='bold', color='#0d47a1')
    for i, (ref, desc) in enumerate(standards):
        col = i % 4
        row = i // 4
        x = 0.95 + col * 4.55
        y = 11.00 - row * 0.65
        ax.text(x, y, '•', fontsize=15, color='#01579b', va='center')
        ax.text(x+0.22, y+0.10, ref, fontsize=10.5, fontweight='bold',
                color='#263238', va='center')
        ax.text(x+0.22, y-0.20, desc, fontsize=9.0,
                color='#546e7a', fontstyle='italic', va='center')

    # ── Connection arrows (vertical, between tiers) ──
    arrow(2.7, 4.4, 2.7, 3.7, 'SGP4 / TR 38.811',  '#bf360c')
    arrow(6.8, 4.4, 6.8, 3.7, 'NR PHY / MAC',     '#1b5e20')
    arrow(11.1, 4.4, 11.1, 3.7, 'L1/L2 split',     '#4a148c')
    arrow(15.4, 4.4, 15.4, 3.7, 'Gym IPC',         '#311b92')

    arrow(2.7, 2.3, 2.7, 1.8, '', '#01579b')
    arrow(7.5, 2.3, 7.5, 1.8, '', '#1565c0')
    arrow(12.4, 2.3, 12.4, 1.8, '', '#283593')
    arrow(17.2, 2.3, 17.2, 1.8, '', '#3949ab')

    # tier-2 to tier-1
    arrow(2.7, 7.3, 2.7, 6.6, '', '#0277bd')
    arrow(7.5, 7.3, 7.5, 6.6, '', '#00695c')
    arrow(12.4, 7.3, 12.4, 6.6, '', '#5e35b1')
    arrow(17.2, 7.3, 17.2, 6.6, '', '#c62828')

    # ── Bottom info ──
    ax.text(10, 0.25,
            'Build:  ./ns3 configure --enable-tests --enable-examples  →  ./ns3 build  →  ./ns3 run "<example>"',
            fontsize=11.5, ha='center', color='#37474f', fontstyle='italic')

    path = os.path.join(OUT_DIR, 'ns3_ntn_toolkit_architecture.png')
    fig.savefig(path, dpi=180, bbox_inches='tight', facecolor='white',
                pad_inches=0.25)
    plt.close(fig)
    print(f"  → {path} ({os.path.getsize(path)/1024:.0f} KB)")


# =============================================================================
# Helper: realistic UE generator on a Europe / North-Atlantic map
# =============================================================================

class RealisticUe:
    def __init__(self, uid, cls, lat, lon, vN, vE, color, alt_m=0):
        self.id = uid
        self.cls = cls
        self.lat = lat
        self.lon = lon
        self.vN = vN
        self.vE = vE
        self.alt_m = alt_m
        self.color = color
        self.heading = np.arctan2(vE, vN)
        self.speed = np.hypot(vN, vE)
        self.next_change = 8.0  # for pedestrian / vehicular


def step_ue(ue, dt, rng):
    if ue.cls in ('static', 'iot'):
        return
    if ue.cls == 'pedestrian':
        ue.next_change -= dt
        if ue.next_change <= 0:
            ue.heading = rng.uniform(0, 2*np.pi)
            ue.speed = rng.normal(1.2, 0.4)
            ue.speed = float(np.clip(ue.speed, 0.3, 3.0))
            ue.next_change = rng.uniform(5, 15)
            ue.vE = ue.speed * np.sin(ue.heading)
            ue.vN = ue.speed * np.cos(ue.heading)
    elif ue.cls == 'vehicular':
        # Gauss-Markov heading
        alpha = np.exp(-dt/60.0)
        sigma = 0.2
        ue.heading = alpha*ue.heading + sigma*np.sqrt(1-alpha*alpha)*rng.normal()
        ue.vE = ue.speed * np.sin(ue.heading)
        ue.vN = ue.speed * np.cos(ue.heading)

    M = 111320.0
    ue.lat += ue.vN * dt / M
    ue.lon += ue.vE * dt / (M * max(np.cos(np.radians(ue.lat)), 0.05))


def make_population(rng):
    ues = []
    # 2 static handhelds in major cities
    ues.append(RealisticUe(0, 'static', 51.5, -0.1, 0, 0, PALETTE['static']))      # London
    ues.append(RealisticUe(1, 'static', 48.85, 2.35, 0, 0, PALETTE['static']))     # Paris
    # 2 pedestrians
    ues.append(RealisticUe(2, 'pedestrian', 52.5, 13.4, 0.8, 0.9, PALETTE['pedestrian'])) # Berlin
    ues.append(RealisticUe(3, 'pedestrian', 41.9, 12.5,-0.5, 1.1, PALETTE['pedestrian'])) # Rome
    # 2 vehicular (highway)
    ues.append(RealisticUe(4, 'vehicular', 50.0, 8.0, 22.0, -10.0, PALETTE['vehicular']))
    ues.append(RealisticUe(5, 'vehicular', 47.0, 19.0, 12.0,  18.0, PALETTE['vehicular']))
    # 2 HST (linear track)
    ues.append(RealisticUe(6, 'hst', 45.5, -1.0, 0.0,  82.0, PALETTE['hst']))     # France-bound east
    ues.append(RealisticUe(7, 'hst', 50.5, 25.0, 30.0,-78.0, PALETTE['hst']))     # NE→SW
    # 2 maritime (Atlantic)
    ues.append(RealisticUe(8, 'maritime', 48.0, -10.0, 0.5, 11.0, PALETTE['maritime']))
    ues.append(RealisticUe(9, 'maritime', 45.0, -20.0,-0.3,-10.5, PALETTE['maritime']))
    # 2 aviation (FL350)
    ues.append(RealisticUe(10, 'aviation', 50.0, -15.0, 25.0, 235.0,
                            PALETTE['aviation'], alt_m=10668))   # ATL→EUR
    ues.append(RealisticUe(11, 'aviation', 40.0, 5.0, 55.0,-230.0,
                            PALETTE['aviation'], alt_m=10668))   # EUR→ATL
    # 2 IoT (rural sensors)
    ues.append(RealisticUe(12, 'iot', 55.0, 12.0, 0, 0, PALETTE['iot']))
    ues.append(RealisticUe(13, 'iot', 38.5, -7.5, 0, 0, PALETTE['iot']))
    return ues


# =============================================================================
# Constellation geometry helper — simple Walker-Star at 780 km
# =============================================================================

def constellation_positions(t, n_planes=6, sats_per_plane=11, alt_km=780,
                             inc_deg=86.4):
    R = 6371.0
    r = R + alt_km
    mu = 398600.4418
    period = 2*np.pi*np.sqrt(r**3/mu)
    omega = 2*np.pi/period
    omega_e = 2*np.pi/86164.0
    inc = np.radians(inc_deg)
    out = []
    for p in range(n_planes):
        raan = 2*np.pi*p/n_planes
        for s in range(sats_per_plane):
            ma = 2*np.pi*s/sats_per_plane + omega*t
            x = r*np.cos(ma)
            y = r*np.sin(ma)
            x_e = x*np.cos(raan) - y*np.cos(inc)*np.sin(raan)
            y_e = x*np.sin(raan) + y*np.cos(inc)*np.cos(raan)
            z_e = y*np.sin(inc)
            lat = np.degrees(np.arcsin(z_e/r))
            lon = np.degrees(np.arctan2(y_e, x_e)) - np.degrees(omega_e*t)
            lon = ((lon + 180) % 360) - 180
            out.append((lat, lon))
    return out


def best_serving(ue_lat, ue_lon, sats):
    """Return (sat_idx, elev_deg) of the highest-elev visible sat (>10°)."""
    best = (-1, -90.0)
    R = 6371.0
    for i, (lat, lon) in enumerate(sats):
        # great-circle central angle
        c = (np.sin(np.radians(ue_lat))*np.sin(np.radians(lat)) +
             np.cos(np.radians(ue_lat))*np.cos(np.radians(lat))*
             np.cos(np.radians(lon-ue_lon)))
        c = float(np.clip(c, -1, 1))
        gamma = np.arccos(c)  # rad
        d = R*np.sqrt(1 + ((R+780)/R)**2 - 2*((R+780)/R)*np.cos(gamma))
        # elevation
        sin_elev = ((R+780)*np.cos(gamma) - R) / d
        elev = np.degrees(np.arcsin(np.clip(sin_elev, -1, 1)))
        if elev > best[1]:
            best = (i, elev)
    if best[1] < 10.0:
        return (-1, best[1])
    return best


# =============================================================================
# 2. Realistic mobility GIF
# =============================================================================

def generate_realistic_mobility_gif(n_frames=80, sim_step_s=8.0):
    print("[2/3] realistic mobility GIF (this may take ~60 s)...")

    # Two-panel layout: large map on the left, dedicated legend axis on the
    # right.  This keeps the legend completely outside the data area, so it
    # never overlaps UEs or satellite scatter no matter how the map evolves.
    fig = plt.figure(figsize=(18, 10), facecolor='#050a18')
    gs = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.0],
                          left=0.05, right=0.985, top=0.93, bottom=0.07,
                          wspace=0.04)
    ax     = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])
    ax.set_facecolor('#050a18')
    ax_leg.set_facecolor('#0a1628')

    rng = np.random.default_rng(42)
    ues = make_population(rng)

    # map bounds — Europe + North Atlantic
    lon_min, lon_max = -25, 35
    lat_min, lat_max = 30, 65

    def draw_static_map():
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_aspect((lon_max-lon_min) / (lat_max-lat_min) /
                       (np.cos(np.radians(50.0))))
        for ll in range(lon_min, lon_max+1, 10):
            ax.axvline(ll, color='#1a2744', lw=0.6, alpha=0.6, zorder=1)
        for ll in range(lat_min, lat_max+1, 10):
            ax.axhline(ll, color='#1a2744', lw=0.6, alpha=0.6, zorder=1)
        europe = [
            (-10, 36), (-9, 43), (-2, 43), (1, 49), (-5, 51),
            (-1, 60), (10, 65), (30, 65), (35, 50), (28, 41),
            (20, 38), (12, 38), (5, 36), (-6, 36)
        ]
        ax.add_patch(Polygon(europe, fc='#0a1f3a', ec='#1a3556',
                             lw=1.2, alpha=0.85, zorder=2))
        ax.set_xlabel('Longitude (°)', color='#90a4ae', fontsize=13)
        ax.set_ylabel('Latitude (°)',  color='#90a4ae', fontsize=13)
        ax.tick_params(colors='#cfd8dc', labelsize=11)
        for s in ax.spines.values():
            s.set_color('#1a2744')

    legend_classes = [
        ('Static handheld', 'static',     _phone),
        ('Pedestrian',      'pedestrian', _pedestrian),
        ('Vehicular',       'vehicular',  _car),
        ('High-speed train','hst',        _train),
        ('Maritime',        'maritime',   _ship),
        ('Aviation FL350',  'aviation',   _plane),
        ('IoT sensor',      'iot',        _iot),
    ]

    def draw_legend_panel():
        """Static legend in its own axis — drawn once per frame after clear."""
        ax_leg.clear()
        ax_leg.set_facecolor('#0a1628')
        ax_leg.set_xlim(0, 10); ax_leg.set_ylim(0, 24)
        ax_leg.set_xticks([]); ax_leg.set_yticks([])
        for s in ax_leg.spines.values():
            s.set_color('#37474f'); s.set_linewidth(1.2)

        # ── header ──
        ax_leg.text(5.0, 22.6, 'UE Classes',
                    color='#e8eaf6', fontsize=15, fontweight='bold',
                    ha='center', va='center')
        ax_leg.text(5.0, 21.6, '(3GPP TR 38.811 §6.1.1.1)',
                    color='#90a4ae', fontsize=10, fontstyle='italic',
                    ha='center', va='center')
        ax_leg.plot([0.6, 9.4], [20.7, 20.7], color='#37474f', lw=0.8)

        # ── 7 class rows, generously spaced ──
        y0, dy = 19.2, 2.20
        for i, (lab, key, drawer) in enumerate(legend_classes):
            yy = y0 - i*dy
            drawer(ax_leg, 1.6, yy, 1.55, PALETTE[key])
            ax_leg.text(3.4, yy, lab, color='#e8eaf6', fontsize=12.5,
                        va='center', fontweight='bold')

        # ── footer: satellite + link key ──
        ax_leg.plot([0.6, 9.4], [3.5, 3.5], color='#37474f', lw=0.8)
        ax_leg.scatter(1.6, 2.6, marker='^', s=240,
                       c=PALETTE['sat'], edgecolors='white', linewidths=1.0)
        ax_leg.text(3.4, 2.6, 'LEO sat (780 km)',
                    color='#e8eaf6', fontsize=11.5, va='center',
                    fontweight='bold')
        ax_leg.plot([1.0, 2.2], [1.3, 1.3],
                    color=PALETTE['gnd_link'], lw=1.4, alpha=0.85)
        ax_leg.text(3.4, 1.3, 'serving link',
                    color='#cfd8dc', fontsize=11, va='center')
        ax_leg.text(5.0, 0.4,
                    f'14 UEs · 6 × 11 = 66 sats',
                    color='#78909c', fontsize=9.5, fontstyle='italic',
                    ha='center')

    def step(frame):
        ax.clear()
        draw_static_map()
        draw_legend_panel()

        t_sec = frame * sim_step_s
        for ue in ues:
            step_ue(ue, sim_step_s, rng)

        sats = constellation_positions(t_sec)
        sat_lats = np.array([s[0] for s in sats])
        sat_lons = np.array([s[1] for s in sats])
        m = ((sat_lats > lat_min-5) & (sat_lats < lat_max+5) &
             (sat_lons > lon_min-5) & (sat_lons < lon_max+5))
        ax.scatter(sat_lons[m], sat_lats[m], marker='^', s=180,
                    c=PALETTE['sat'], edgecolors='white', linewidths=1.0,
                    zorder=8)

        for ue in ues:
            si, elev = best_serving(ue.lat, ue.lon, sats)
            if si >= 0 and m[si]:
                ax.plot([ue.lon, sats[si][1]], [ue.lat, sats[si][0]],
                         color=PALETTE['gnd_link'], lw=0.9, alpha=0.55,
                         zorder=5)
            # bigger, scenario-aware icon scaling — fills the map cleanly now
            scale = 2.4 if ue.cls in ('hst', 'maritime', 'aviation') else 2.0
            ICONS[ue.cls](ax, ue.lon, ue.lat, scale, ue.color)

        ax.set_title(
            f'ns3-ntn-toolkit — Realistic NTN UE Mobility   ·   t = {t_sec:6.1f} s',
            color='#e8eaf6', fontsize=17, fontweight='bold', pad=14)

    writer = PillowWriter(fps=8)
    anim = FuncAnimation(fig, step, frames=n_frames, interval=125)
    path = os.path.join(OUT_DIR, 'ntn_realistic_mobility.gif')
    anim.save(path, writer=writer, dpi=100, savefig_kwargs={'facecolor': '#050a18'})
    plt.close(fig)
    print(f"  → {path} ({os.path.getsize(path)/1024/1024:.1f} MB, {n_frames} frames)")


# =============================================================================
# 3. Per-class handover behaviour GIF
# =============================================================================

def generate_handover_realistic_gif(n_frames=80, sim_step_s=8.0):
    print("[3/3] handover-with-realistic-UE GIF...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.5),
                              facecolor='#050a18',
                              gridspec_kw={'width_ratios':[3, 2]})
    ax = axes[0]
    ax_kpi = axes[1]
    ax.set_facecolor('#050a18')
    ax_kpi.set_facecolor('#0a1628')

    rng = np.random.default_rng(2026)
    ues = make_population(rng)

    # tracking per-UE: serving sat, HO count
    serving = [-1] * len(ues)
    ho_count = [0] * len(ues)
    ho_history = []  # list of (t, ue_id, class)

    lon_min, lon_max = -25, 35
    lat_min, lat_max = 30, 65

    def draw_map():
        ax.clear()
        ax.set_xlim(lon_min, lon_max); ax.set_ylim(lat_min, lat_max)
        ax.set_aspect((lon_max-lon_min) / (lat_max-lat_min) /
                       np.cos(np.radians(50.0)))
        for ll in range(lon_min, lon_max+1, 10):
            ax.axvline(ll, color='#1a2744', lw=0.4, alpha=0.5, zorder=1)
        for ll in range(lat_min, lat_max+1, 10):
            ax.axhline(ll, color='#1a2744', lw=0.4, alpha=0.5, zorder=1)
        europe = [
            (-10, 36), (-9, 43), (-2, 43), (1, 49), (-5, 51),
            (-1, 60), (10, 65), (30, 65), (35, 50), (28, 41),
            (20, 38), (12, 38), (5, 36), (-6, 36)
        ]
        ax.add_patch(Polygon(europe, fc='#0a1f3a', ec='#1a3556',
                             lw=0.8, alpha=0.85, zorder=2))
        ax.tick_params(colors='#607d8b', labelsize=8)
        for s in ax.spines.values():
            s.set_color('#1a2744')

    def update(frame):
        nonlocal serving, ho_count
        t_sec = frame * sim_step_s
        for ue in ues:
            step_ue(ue, sim_step_s, rng)

        sats = constellation_positions(t_sec)
        sat_lats = np.array([s[0] for s in sats])
        sat_lons = np.array([s[1] for s in sats])

        draw_map()
        m = ((sat_lats > lat_min-5) & (sat_lats < lat_max+5) &
             (sat_lons > lon_min-5) & (sat_lons < lon_max+5))
        ax.scatter(sat_lons[m], sat_lats[m], marker='^', s=85,
                   c=PALETTE['sat'], edgecolors='white', linewidths=0.6,
                   zorder=8)

        for i, ue in enumerate(ues):
            si, elev = best_serving(ue.lat, ue.lon, sats)
            if si >= 0 and si != serving[i]:
                if serving[i] != -1:
                    ho_count[i] += 1
                    ho_history.append((t_sec, ue.id, ue.cls))
                serving[i] = si
            if si >= 0 and m[si]:
                ax.plot([ue.lon, sats[si][1]], [ue.lat, sats[si][0]],
                         color=PALETTE['gnd_link'], lw=0.5, alpha=0.5, zorder=5)
            scale = 2.0 if ue.cls in ('hst', 'maritime', 'aviation') else 1.6
            ICONS[ue.cls](ax, ue.lon, ue.lat, scale, ue.color)

        ax.set_title(
            f'TTE-aware CHO  ·  t = {t_sec:6.1f} s  ·  total HO = {sum(ho_count)}',
            color='#e8eaf6', fontsize=14, fontweight='bold', pad=10)

        # ── KPI panel ──
        ax_kpi.clear()
        ax_kpi.set_facecolor('#0a1628')
        ax_kpi.set_title('Per-class handover counts',
                          color='#e8eaf6', fontsize=12, fontweight='bold',
                          pad=10)
        # aggregate per class
        per_cls = {c: 0 for c in PALETTE if c in ICONS}
        for j, ue in enumerate(ues):
            per_cls.setdefault(ue.cls, 0)
            per_cls[ue.cls] = per_cls.get(ue.cls, 0) + ho_count[j]
        labels = ['static', 'pedestrian', 'vehicular', 'hst',
                  'maritime', 'aviation', 'iot']
        vals = [per_cls.get(l, 0) for l in labels]
        bars = ax_kpi.barh(range(len(labels)), vals,
                            color=[PALETTE[l] for l in labels],
                            edgecolor='white', linewidth=0.5)
        for b, v in zip(bars, vals):
            ax_kpi.text(b.get_width()+0.2, b.get_y()+b.get_height()/2,
                         str(int(v)), color='white', va='center', fontsize=10)
        ax_kpi.set_yticks(range(len(labels)))
        ax_kpi.set_yticklabels(labels, color='#cfd8dc', fontsize=9)
        ax_kpi.tick_params(colors='#607d8b')
        for s in ax_kpi.spines.values():
            s.set_color('#1a2744')
        ax_kpi.set_xlim(0, max(max(vals)+1, 6))
        ax_kpi.set_xlabel('handover count',
                           color='#90a4ae', fontsize=9)
        # caption above the chart so the bottom margin stays clean
        fig.text(0.74, 0.94,
                 'fewer HOs ← static · IoT          aviation · HST → more HOs',
                 color='#78909c', fontsize=8.5, fontstyle='italic',
                 ha='center')

    writer = PillowWriter(fps=8)
    anim = FuncAnimation(fig, update, frames=n_frames, interval=125)
    path = os.path.join(OUT_DIR, 'ntn_handover_realistic.gif')
    anim.save(path, writer=writer, dpi=100,
              savefig_kwargs={'facecolor': '#050a18'})
    plt.close(fig)
    print(f"  → {path} ({os.path.getsize(path)/1024/1024:.1f} MB, {n_frames} frames)")


# =============================================================================
if __name__ == '__main__':
    print('=' * 64)
    print('  ns3-ntn-toolkit — realistic visualisation generator')
    print('=' * 64)
    generate_toolkit_architecture()
    generate_realistic_mobility_gif()
    generate_handover_realistic_gif()
    print('\nAll artefacts written to', OUT_DIR)
