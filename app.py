import math
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="Active & Passive Earth Pressure Diagram", layout="wide")
st.markdown("## 🧱 Active & Passive Earth Pressure Diagram (English Units)")
st.caption("Rankine theory with cohesion, dual water tables, excavation depth, and optional seismic/hydrodynamic terms.")

gamma_w = 62.4  # pcf

# ---------------- Sidebar inputs ----------------
st.sidebar.header("Input Parameters")

# Layers
n_layers = st.sidebar.number_input("Number of soil layers", 1, 10, 3, 1)

# Water tables (depth below ground surface; 0 at top)
gw_active_depth = st.sidebar.number_input("Water Table Depth (Active Side, ft)", value=6.0, min_value=0.0)
gw_pass_depth   = st.sidebar.number_input("Water Table Depth (Passive Side, ft)", value=8.0, min_value=0.0)

# Excavation (passive side)
exc_depth = st.sidebar.number_input("Excavation Depth (Passive Side, ft)", value=10.0, min_value=0.0)

# Surcharge
surcharge = st.sidebar.number_input("Surcharge (psf)", value=300.0, min_value=0.0)

# Seismic & hydrodynamic (simple, user-controlled)
st.sidebar.markdown("---")
st.sidebar.subheader("Optional Dynamic Terms")
Kh_soil = st.sidebar.number_input("Seismic horizontal coeff. for soil, Kh (dimensionless)", value=0.0, min_value=0.0, step=0.01)
ah      = st.sidebar.number_input("Hydrodynamic input coeff. (approx. a_h)", value=0.0, min_value=0.0, step=0.01)
Cw      = st.sidebar.number_input("Hydrodynamic scale factor Cw", value=0.875, min_value=0.0, step=0.025, help="Used as multiplier in p≈Cw·γw·ah·(water depth - y).")

st.sidebar.markdown("---")
st.sidebar.subheader("Soil Layer Properties")

# Layer editors (thickness + properties)
layers = []
for i in range(1, n_layers + 1):
    st.sidebar.markdown(f"**Layer {i}**")
    t   = st.sidebar.number_input(f"Thickness L{i} (ft)",  min_value=0.01, value=5.0, step=0.5, key=f"t{i}")
    phi = st.sidebar.number_input(f"ϕ L{i} (°)",            min_value=0.0,  value=30.0, step=0.5, key=f"phi{i}")
    c   = st.sidebar.number_input(f"c L{i} (psf)",          min_value=0.0,  value=0.0,  step=5.0, key=f"c{i}")
    g_d = st.sidebar.number_input(f"γ_dry L{i} (pcf)",      min_value=1.0,  value=110.0,step=1.0, key=f"gdy{i}")
    g_s = st.sidebar.number_input(f"γ_sat L{i} (pcf)",      min_value=1.0,  value=120.0,step=1.0, key=f"gsat{i}")
    layers.append(dict(thk=t, phi=phi, c=c, gamma_dry=g_d, gamma_sat=g_s))

st.sidebar.markdown("")

# Button (optional; app is reactive anyway)
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
    # below last: use last
    L = layers[-1]
    return L["phi"], L["c"], L["gamma_dry"], L["gamma_sat"]

def vertical_effective_stress(z, gw_z, surcharge):
    """σ'v at elevation z using piecewise γ' and surcharge."""
    sigma = surcharge
    # accumulate layer-by-layer from surface to z
    current = 0.0
    for i in range(n_layers):
        z_top = tops[i]
        z_bot_i = bottoms[i]
        phi, c, g_dry, g_sat = layers[i]["phi"], layers[i]["c"], layers[i]["gamma_dry"], layers[i]["gamma_sat"]

        # segment contributing = from current elevation (start) down to either z or bottom of this layer
        seg_top = max(z, z_top)  # shallower (less negative)
        seg_bot = z_top
        if z <= z_bot_i:
            seg_top = z  # we're entirely below this layer
            seg_bot = z_top

        # thickness inside layer from z_top to z (if z is deeper)
        dz = max(0.0, z_top - max(z, z_bot_i))
        if dz > 0:
            # average unit weight with water table
            # split by gw crossing if needed
            parts = []
            # segment 1: above/below GW?
            if (z_top >= gw_z) and (max(z, z_bot_i) >= gw_z):
                # entire segment above GW -> γ = γ_dry
                parts.append((dz, g_dry))
            elif (z_top <= gw_z) and (max(z, z_bot_i) <= gw_z):
                # entire segment below GW -> γ' = γ_sat - γw
                parts.append((dz, g_sat - gamma_w))
            else:
                # crosses GW: split
                dz1 = z_top - gw_z
                dz2 = max(z, z_bot_i) - gw_z
                dz_above = max(0.0, dz1)
                dz_below = max(0.0, -dz2)
                if dz_above > 0: parts.append((dz_above, g_dry))
                if dz_below > 0: parts.append((dz_below, g_sat - gamma_w))

            for thk, g in parts:
                sigma += g * thk

        if z >= z_bot_i:
            # we've reached/ passed z; stop
            break

    return sigma

# ---------------- Profiles ----------------
npts = 600
z_profile = np.linspace(0.0, z_bot, npts)  # elevations (negative down)
depth = -z_profile                         # positive depth for plotting

# σ'v for each side
sigma_eff_active = np.array([vertical_effective_stress(z, gw_active_z, surcharge) for z in z_profile])
sigma_eff_pass   = np.array([vertical_effective_stress(z, gw_pass_z,   surcharge) for z in z_profile])

# strength params per point for Rankine
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

# Active lateral (soil)
p_active_eff = np.maximum(0.0, Ka_prof * sigma_eff_active - 2.0*np.sqrt(Ka_prof)*c_prof)

# Passive lateral (soil), zero above cut, referenced to σ' at cut on passive side
sigma_cut_pass = float(vertical_effective_stress(cut_z, gw_pass_z, surcharge))
p_passive_eff  = Kp_prof * np.maximum(0.0, sigma_eff_pass - sigma_cut_pass) + 2.0*np.sqrt(Kp_prof)*c_prof
p_passive_eff[z_profile > cut_z] = 0.0
p_passive_eff = np.maximum(0.0, p_passive_eff)

# Hydrostatic (water) components
hs_active = gamma_w * np.maximum(0.0, (gw_active_z - z_profile))  # linear below GW
hs_pass   = gamma_w * np.maximum(0.0, (gw_pass_z   - z_profile))

# Seismic soil increment (simple): Kh * σ′v
p_seis_active = Kh_soil * sigma_eff_active
p_seis_pass   = Kh_soil * sigma_eff_pass
p_seis_pass[z_profile > cut_z] = 0.0  # typically no soil at the excavated face above cut

# Hydrodynamic (very simple triangular approx. below passive GW)
# p_hdyn(z) ≈ Cw * γw * ah * (water depth below passive GW), zero above GW or above cut
water_depth_at = np.maximum(0.0, gw_pass_z - z_profile)
p_hdyn = Cw * gamma_w * ah * water_depth_at
p_hdyn[(z_profile > gw_pass_z) | (z_profile > cut_z)] = 0.0

# Totals (for plotting as single curves)
p_active_total  = p_active_eff  + hs_active + p_seis_active
p_passive_total = p_passive_eff + hs_pass   + p_seis_pass + p_hdyn

# ---------------- Results snippet table ----------------
rows = []
for i in range(n_layers):
    zt, zb = tops[i], bottoms[i]
    phi, c, gd, gs = layers[i]["phi"], layers[i]["c"], layers[i]["gamma_dry"], layers[i]["gamma_sat"]
    ka, kp = Ka(phi), Kp(phi)
    # interface values (active side)
    st_act = vertical_effective_stress(zt, gw_active_z, surcharge)
    sb_act = vertical_effective_stress(zb, gw_active_z, surcharge)
    at = max(0.0, ka*st_act - 2*math.sqrt(ka)*c)
    ab = max(0.0, ka*sb_act - 2*math.sqrt(ka)*c)
    # interface values (passive side)
    st_pas = vertical_effective_stress(zt, gw_pass_z, surcharge)
    sb_pas = vertical_effective_stress(zb, gw_pass_z, surcharge)
    pt = 0.0 if zt > cut_z else max(0.0, kp*max(0.0, st_pas - sigma_cut_pass) + 2*math.sqrt(kp)*c)
    pb = 0.0 if zb > cut_z else max(0.0, kp*max(0.0, sb_pas - sigma_cut_pass) + 2*math.sqrt(kp)*c)
    rows.append([f"L{i+1}", zt, zb, phi, c, ka, kp, at, ab, pt, pb])

df = pd.DataFrame(rows, columns=["Layer","Top","Bottom","phi","c","Ka","Kp","Active Top","Active Bottom","Passive Top","Passive Bottom"])

with col2:
    st.markdown("### Interface Results (soil-only lateral)")
    st.dataframe(df, use_container_width=True)

# ---------------- Diagram (Passive left, Active right) ----------------
fig, ax = plt.subplots(figsize=(7.6,7.6))
# soil effective components
ax.plot(-p_passive_eff, depth, 'r-', lw=1.8, label='Passive (soil, eff.)')
ax.plot( p_active_eff,  depth, 'b-', lw=1.8, label='Active (soil, eff.)')
ax.fill_betweenx(depth, 0, -p_passive_eff, color='red',  alpha=0.15)
ax.fill_betweenx(depth, 0,  p_active_eff,  color='blue', alpha=0.15)

# optional adds (hydrostatic, seismic, hydrodynamic) on passive and active as faint lines
ax.plot(-(p_passive_total), depth, color='maroon', lw=1.0, ls='--', label='Passive total (soil+water+dynamic)')
ax.plot( (p_active_total), depth, color='navy',   lw=1.0, ls='--', label='Active total (soil+water+dynamic)')

# excavation line
ax.axhline(exc_depth, color='k', lw=2)
ax.text(0, exc_depth+0.2, f'Excavation (z={exc_depth:.1f} ft)', ha='center', va='bottom')

# layer bottoms
for zb in bottoms:
    ax.axhline(-zb, color='k', ls='--', lw=0.7, alpha=0.4)

ax.set_xlabel("Lateral Pressure (psf)")
ax.set_ylabel("Depth (ft)")
ax.grid(True, ls="--", alpha=0.3)
ax.invert_yaxis()

xmax = max(1.0, p_active_total.max(), p_passive_total.max())
ax.set_xlim(-1.1*xmax, 1.1*xmax)
ax.set_title("Active (Right) and Passive (Left) Earth Pressure Diagram")
ax.legend(loc="upper right", fontsize=9)

st.pyplot(fig)

