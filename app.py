# ------------------ app.py (complete) ------------------
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Active & Passive Earth Pressure Diagram", layout="wide")
st.markdown("## 🧱 Active & Passive Earth Pressure Diagram (English Units)")
st.caption("Rankine theory with cohesion, dual water tables, excavation depth, and optional seismic/hydrodynamic terms.")

gamma_w = 62.4  # pcf

# ---------------- Sidebar inputs ----------------
st.sidebar.header("Input Parameters")

n_layers = st.sidebar.number_input("Number of soil layers", 1, 10, 3, 1)

# depths are positive downward from the ground surface (we convert to elevation z later)
gw_active_depth = st.sidebar.number_input("Water Table Depth (Active Side, ft)", value=6.0, min_value=0.0)
gw_pass_depth   = st.sidebar.number_input("Water Table Depth (Passive Side, ft)", value=8.0, min_value=0.0)

exc_depth = st.sidebar.number_input("Excavation Depth (Passive Side, ft)", value=10.0, min_value=0.0)
surcharge = st.sidebar.number_input("Surcharge (psf)", value=300.0, min_value=0.0)

st.sidebar.markdown("---")
st.sidebar.subheader("Optional Dynamic Terms")
Kh_soil = st.sidebar.number_input("Seismic horizontal coeff. for soil, Kh", value=0.0, min_value=0.0, step=0.01)
ah      = st.sidebar.number_input("Hydrodynamic input coeff. (a_h)", value=0.0, min_value=0.0, step=0.01)
Cw      = st.sidebar.number_input("Hydrodynamic scale factor Cw", value=0.875, min_value=0.0, step=0.025)

st.sidebar.markdown("---")
st.sidebar.subheader("Soil Layer Properties")

layers = []
for i in range(1, n_layers + 1):
    st.sidebar.markdown(f"**Layer {i}**")
    t   = st.sidebar.number_input(f"Thickness L{i} (ft)",  min_value=0.01, value=5.0,  step=0.5, key=f"t{i}")
    phi = st.sidebar.number_input(f"ϕ L{i} (°)",           min_value=0.0,  value=30.0, step=0.5, key=f"phi{i}")
    c   = st.sidebar.number_input(f"c L{i} (psf)",         min_value=0.0,  value=0.0,  step=5.0, key=f"c{i}")
    g_d = st.sidebar.number_input(f"γ_dry L{i} (pcf)",     min_value=1.0,  value=110.0,step=1.0, key=f"gdy{i}")
    g_s = st.sidebar.number_input(f"γ_sat L{i} (pcf)",     min_value=1.0,  value=120.0,step=1.0, key=f"gsat{i}")
    layers.append(dict(thk=t, phi=phi, c=c, gamma_dry=g_d, gamma_sat=g_s))

st.sidebar.button("Compute")

# ---------------- Derived geometry ----------------
# Elevation z: 0 at ground, negative downward. Build tops/bottoms from thicknesses.
tops = [0.0]
for i in range(n_layers):
    tops.append(tops[-1] - layers[i]["thk"])
tops = tops[:-1]
bottoms = [tops[i] - layers[i]["thk"] for i in range(n_layers)]

df_layers = pd.DataFrame({
    "Layer": [f"L{i+1}" for i in range(n_layers)],
    "Top Elev.": tops,
    "Bottom Elev.": bottoms,
    "ϕ (deg)": [L["phi"] for L in layers],
    "c (psf)": [L["c"] for L in layers],
    "γ_dry (pcf)": [L["gamma_dry"] for L in layers],
    "γ_sat (pcf)": [L["gamma_sat"] for L in layers],
})

col1, col2 = st.columns([1,1])
with col1:
    st.markdown("### Layers (computed elevations)")
    st.dataframe(df_layers, use_container_width=True)

z_bot = bottoms[-1]
gw_active_z = -gw_active_depth
gw_pass_z   = -gw_pass_depth
cut_z       = -exc_depth

# ---------------- Utilities ----------------
def Ka(phi): return math.tan(math.radians(45 - phi/2.0))**2
def Kp(phi): return math.tan(math.radians(45 + phi/2.0))**2

def phi_c_at_z(z):
    """Return (phi, c, gamma_dry, gamma_sat) for the layer containing elevation z."""
    for i in range(n_layers):
        if (z <= tops[i]) and (z >= bottoms[i]):
            L = layers[i]
            return L["phi"], L["c"], L["gamma_dry"], L["gamma_sat"]
    L = layers[-1]
    return L["phi"], L["c"], L["gamma_dry"], L["gamma_sat"]

def vertical_effective_stress(z, gw_z, surcharge):
    """σ'v at elevation z using piecewise γ' and surcharge (z is ≤ 0)."""
    sigma = surcharge
    for i in range(n_layers):
        z_top = tops[i]
        z_bot_i = bottoms[i]
        _, _, g_dry, g_sat = phi_c_at_z(z_top)

        if z_top <= z:
            # we are already deeper than/at this top; continue
            continue

        # thickness that contributes from this layer = from z_top down to max(z, z_bot_i)
        seg_top = z_top
        seg_bot = max(z, z_bot_i)
        dz = max(0.0, seg_top - seg_bot)
        if dz <= 0:
            continue

        # split at GW if intersected
        above_top = seg_top >= gw_z
        above_bot = seg_bot >= gw_z
        if above_top and above_bot:
            sigma += g_dry * dz
        elif (not above_top) and (not above_bot):
            sigma += (g_sat - gamma_w) * dz
        else:
            # crosses
            dz_above = max(0.0, seg_top - gw_z)
            dz_below = max(0.0, gw_z - seg_bot)
            if dz_above > 0: sigma += g_dry * dz_above
            if dz_below > 0: sigma += (g_sat - gamma_w) * dz_below

        if z >= z_bot_i:
            break

    return sigma

# ---------------- Profiles ----------------
npts = 600
z_profile = np.linspace(0.0, z_bot, npts)  # elevations (negative down)
depth = -z_profile                         # positive depth for plotting

sigma_eff_active = np.array([vertical_effective_stress(z, gw_active_z, surcharge) for z in z_profile])
sigma_eff_pass   = np.array([vertical_effective_stress(z, gw_pass_z,   surcharge) for z in z_profile])

phi_prof   = np.zeros_like(z_profile)
c_prof     = np.zeros_like(z_profile)
Ka_prof    = np.zeros_like(z_profile)
Kp_prof    = np.zeros_like(z_profile)
for i, z in enumerate(z_profile):
    phi, c, g_d, g_s = phi_c_at_z(z)
    phi_prof[i] = phi
    c_prof[i]   = c
    Ka_prof[i]  = Ka(phi)
    Kp_prof[i]  = Kp(phi)

# Active soil (effective)
p_active_eff = np.maximum(0.0, Ka_prof * sigma_eff_active - 2.0*np.sqrt(Ka_prof)*c_prof)

# Passive soil (effective), zero above cut, referenced to σ' at cut on passive side
sigma_cut_pass = float(vertical_effective_stress(cut_z, gw_pass_z, surcharge))
p_passive_eff  = Kp_prof * np.maximum(0.0, sigma_eff_pass - sigma_cut_pass) + 2.0*np.sqrt(Kp_prof)*c_prof
p_passive_eff[z_profile > cut_z] = 0.0
p_passive_eff = np.maximum(0.0, p_passive_eff)

# Hydrostatic (water) components
hs_active = gamma_w * np.maximum(0.0, (gw_active_z - z_profile))
hs_pass   = gamma_w * np.maximum(0.0, (gw_pass_z   - z_profile))

# Seismic soil increment (simple): Kh * σ′v
p_seis_active = Kh_soil * sigma_eff_active
p_seis_pass   = Kh_soil * sigma_eff_pass
p_seis_pass[z_profile > cut_z] = 0.0

# Hydrodynamic (simple triangular approx. below passive GW)
water_depth_at = np.maximum(0.0, gw_pass_z - z_profile)
p_hdyn = Cw * gamma_w * ah * water_depth_at
p_hdyn[(z_profile > gw_pass_z) | (z_profile > cut_z)] = 0.0

# ---------------- Interface results table (soil-only) ----------------
rows = []
for i in range(n_layers):
    zt, zb = tops[i], bottoms[i]
    phi, c, gd, gs = phi_c_at_z(zt)
    ka, kp = Ka(phi), Kp(phi)
    st_act = vertical_effective_stress(zt, gw_active_z, surcharge)
    sb_act = vertical_effective_stress(zb, gw_active_z, surcharge)
    at = max(0.0, ka*st_act - 2*math.sqrt(ka)*c)
    ab = max(0.0, ka*sb_act - 2*math.sqrt(ka)*c)
    st_pas = vertical_effective_stress(zt, gw_pass_z, surcharge)
    sb_pas = vertical_effective_stress(zb, gw_pass_z, surcharge)
    pt = 0.0 if zt > cut_z else max(0.0, kp*max(0.0, st_pas - sigma_cut_pass) + 2*math.sqrt(kp)*c)
    pb = 0.0 if zb > cut_z else max(0.0, kp*max(0.0, sb_pas - sigma_cut_pass) + 2*math.sqrt(kp)*c)
    rows.append([f"L{i+1}", zt, zb, phi, c, ka, kp, at, ab, pt, pb])

df = pd.DataFrame(rows, columns=[
    "Layer","Top","Bottom","phi","c","Ka","Kp",
    "Active Top","Active Bottom","Passive Top","Passive Bottom"
])

with col2:
    st.markdown("### Interface Results (soil-only lateral)")
    st.dataframe(df, use_container_width=True)
# ----- Totals for convenience (used by the new plot controls) -----
# Active side total = soil (effective) + hydrostatic + seismic
p_active_total = p_active_eff + hs_active + p_seis_active

# Passive side total = soil (effective) + hydrostatic + seismic + hydrodynamic
p_passive_total = p_passive_eff + hs_pass + p_seis_pass + p_hdyn

# Numerical safety: keep non-negative
p_active_total  = np.maximum(0.0, p_active_total)
p_passive_total = np.maximum(0.0, p_passive_total)

# ---------------- Diagram (configurable by side & component) ----------------
st.markdown("### Earth Pressure Diagram — Select What to Show")

with st.expander("Plot controls", expanded=True):
    colA, colB, colC = st.columns([1.2, 1.6, 1.2])
    with colA:
        sides = st.multiselect(
            "Sides to display",
            ["Active (right)", "Passive (left)"],
            default=["Active (right)", "Passive (left)"]
        )
    with colB:
        components = st.multiselect(
            "Components to display",
            ["Soil (effective)", "Hydrostatic", "Seismic (soil)", "Hydrodynamic water", "Total"],
            default=["Soil (effective)"]
        )
        st.caption("Note: Hydrodynamic water applies to the passive (left) side only.")
    with colC:
        shade = st.checkbox("Shade areas", value=False)
        stacked = st.checkbox("Stack components by side", value=False)

# convenience flags
show_active  = "Active (right)" in sides
show_passive = "Passive (left)" in sides

want_eff   = "Soil (effective)"   in components
want_hs    = "Hydrostatic"        in components
want_seis  = "Seismic (soil)"     in components
want_hdyn  = "Hydrodynamic water" in components  # passive only
want_total = "Total"              in components

# if nothing selected, just show a note
if not (show_active or show_passive):
    st.info("Select at least one **Side to display**.")
else:
    if not any([want_eff, want_hs, want_seis, want_hdyn, want_total]):
        st.info("Select at least one **Component**.")
    else:
        fig, ax = plt.subplots(figsize=(8.5, 8.5))

        # ---------- helper to plot one series -----------
        def plot_series(xvals, lab, color, ls='-', lw=2.0, fill_alpha=0.10, mirror=False):
            x = -xvals if mirror else xvals
            ax.plot(x, depth, color=color, ls=ls, lw=lw, label=lab)
            if shade:
                ax.fill_betweenx(depth, 0, x, color=color, alpha=fill_alpha)

        # ---------- ACTIVE (right, positive) -------------
        if show_active:
            stack_offset = np.zeros_like(depth)
            # for stacked: add components cumulatively on each side
            if stacked and not want_total:
                # order: soil -> hydrostatic -> seismic
                if want_eff:
                    plot_series(p_active_eff, "Active – Soil (effective)", 'royalblue')
                    stack_offset = stack_offset + p_active_eff
                if want_hs:
                    plot_series(stack_offset + hs_active, "Active – Hydrostatic (stacked)", 'steelblue', ls='--', lw=1.6)
                    stack_offset = stack_offset + hs_active
                if want_seis:
                    plot_series(stack_offset + p_seis_active, "Active – Seismic (stacked)", 'navy', ls='-.', lw=1.6)
                    stack_offset = stack_offset + p_seis_active
            else:
                if want_eff:
                    plot_series(p_active_eff, "Active – Soil (effective)", 'royalblue')
                if want_hs:
                    plot_series(hs_active, "Active – Hydrostatic", 'steelblue', ls='--', lw=1.6)
                if want_seis:
                    plot_series(p_seis_active, "Active – Seismic (soil)", 'navy', ls='-.', lw=1.6)
                if want_total:
                    plot_series(p_active_total, "Active – Total", 'midnightblue', ls=':', lw=2.2, fill_alpha=0.06)

        # ---------- PASSIVE (left, negative) -------------
        if show_passive:
            stack_offset = np.zeros_like(depth)
            if stacked and not want_total:
                # order: soil -> hydrostatic -> seismic -> hydrodynamic
                if want_eff:
                    plot_series(p_passive_eff, "Passive – Soil (effective)", 'crimson', mirror=True)
                    stack_offset = stack_offset + p_passive_eff
                if want_hs:
                    plot_series(stack_offset + hs_pass, "Passive – Hydrostatic (stacked)", 'firebrick', ls='--', lw=1.6, mirror=True)
                    stack_offset = stack_offset + hs_pass
                if want_seis:
                    plot_series(stack_offset + p_seis_pass, "Passive – Seismic (stacked)", 'darkred', ls='-.', lw=1.6, mirror=True)
                    stack_offset = stack_offset + p_seis_pass
                if want_hdyn:
                    plot_series(stack_offset + p_hdyn, "Passive – Hydrodynamic (stacked)", 'brown', ls=':', lw=1.6, mirror=True)
                    stack_offset = stack_offset + p_hdyn
            else:
                if want_eff:
                    plot_series(p_passive_eff, "Passive – Soil (effective)", 'crimson', mirror=True)
                if want_hs:
                    plot_series(hs_pass, "Passive – Hydrostatic", 'firebrick', ls='--', lw=1.6, mirror=True)
                if want_seis:
                    plot_series(p_seis_pass, "Passive – Seismic (soil)", 'darkred', ls='-.', lw=1.6, mirror=True)
                if want_hdyn:
                    plot_series(p_hdyn, "Passive – Hydrodynamic (water)", 'brown', ls=':', lw=1.6, mirror=True)
                if want_total:
                    plot_series(p_passive_total, "Passive – Total", 'maroon', ls=':', lw=2.2, mirror=True, fill_alpha=0.06)

        # excavation depth and layer lines
        ax.axhline(exc_depth, color='k', lw=2)
        ax.text(-xmax*0.15, exc_depth-0.2, f'Excavation (z={exc_depth:.1f} ft)', 
        ha='right', va='top', fontsize=9, fontweight='bold')

        for zb in bottoms:
            ax.axhline(-zb, color='k', ls='--', lw=0.7, alpha=0.35)

        # axes, limits, legend
        ax.set_xlabel("Lateral Pressure (psf)")
        ax.set_ylabel("Depth (ft)")
        ax.grid(True, ls="--", alpha=0.3)
        ax.invert_yaxis()

        # Pick a symmetric x-limit based on everything we might plot
        xmax_candidates = [1.0]
        if show_active:
            if want_eff:   xmax_candidates.append(p_active_eff.max())
            if want_hs:    xmax_candidates.append(hs_active.max())
            if want_seis:  xmax_candidates.append(p_seis_active.max())
            if want_total: xmax_candidates.append(p_active_total.max())
        if show_passive:
            if want_eff:   xmax_candidates.append(p_passive_eff.max())
            if want_hs:    xmax_candidates.append(hs_pass.max())
            if want_seis:  xmax_candidates.append(p_seis_pass.max())
            if want_hdyn:  xmax_candidates.append(p_hdyn.max())
            if want_total: xmax_candidates.append(p_passive_total.max())
        xmax = max(xmax_candidates)
        ax.set_xlim(-1.1*xmax, 1.1*xmax)

        title_bits = []
        if show_active:  title_bits.append("Active")
        if show_passive: title_bits.append("Passive")
        ax.set_title(" / ".join(title_bits) + " — Selected Components")
        ax.legend(loc="upper right", fontsize=9)
        st.pyplot(fig)
