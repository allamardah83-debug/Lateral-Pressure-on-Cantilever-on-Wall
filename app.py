# ---------------- Diagram (separate components) ----------------
st.markdown("### Earth Pressure Diagram — Components Shown Separately")
shade = st.checkbox("Shade each component", value=False)

fig, ax = plt.subplots(figsize=(8,8))

# ----- Passive components (left; negative x) -----
lns = []
# Passive soil (effective)
lns += ax.plot(-p_passive_eff, depth, color='crimson', lw=2.0, label='Passive – Soil (effective)')
if shade:
    ax.fill_betweenx(depth, 0, -p_passive_eff, color='crimson', alpha=0.12)

# Passive hydrostatic
lns += ax.plot(-hs_pass, depth, color='firebrick', lw=1.5, ls='--', label='Passive – Hydrostatic')
if shade:
    ax.fill_betweenx(depth, 0, -hs_pass, color='firebrick', alpha=0.08)

# Passive seismic soil increment
lns += ax.plot(-p_seis_pass, depth, color='darkred', lw=1.5, ls='-.', label='Passive – Seismic (soil)')
if shade:
    ax.fill_betweenx(depth, 0, -p_seis_pass, color='darkred', alpha=0.06)

# Passive hydrodynamic water
lns += ax.plot(-p_hdyn, depth, color='brown', lw=1.5, ls=':', label='Passive – Hydrodynamic (water)')
if shade:
    ax.fill_betweenx(depth, 0, -p_hdyn, color='brown', alpha=0.06)

# ----- Active components (right; positive x) -----
# Active soil (effective)
lns += ax.plot(p_active_eff, depth, color='royalblue', lw=2.0, label='Active – Soil (effective)')
if shade:
    ax.fill_betweenx(depth, 0, p_active_eff, color='royalblue', alpha=0.12)

# Active hydrostatic
lns += ax.plot(hs_active, depth, color='steelblue', lw=1.5, ls='--', label='Active – Hydrostatic')
if shade:
    ax.fill_betweenx(depth, 0, hs_active, color='steelblue', alpha=0.08)

# Active seismic soil increment
lns += ax.plot(p_seis_active, depth, color='navy', lw=1.5, ls='-.', label='Active – Seismic (soil)')
if shade:
    ax.fill_betweenx(depth, 0, p_seis_active, color='navy', alpha=0.06)

# Excavation line and layer bottoms
ax.axhline(exc_depth, color='k', lw=2)
ax.text(0, exc_depth+0.25, f'Excavation (z={exc_depth:.1f} ft)', ha='center', va='bottom', fontsize=9)

for zb in bottoms:
    ax.axhline(-zb, color='k', ls='--', lw=0.7, alpha=0.35)

# Axes & legend
ax.set_xlabel("Lateral Pressure (psf)")
ax.set_ylabel("Depth (ft)")
ax.grid(True, ls="--", alpha=0.3)
ax.invert_yaxis()

xmax = max(
    1.0,
    p_active_eff.max(), hs_active.max(), p_seis_active.max(),
    p_passive_eff.max(), hs_pass.max(), p_seis_pass.max(), p_hdyn.max()
)
ax.set_xlim(-1.1*xmax, 1.1*xmax)
ax.set_title("Active (Right) and Passive (Left) — Separate Components")

ax.legend(loc="upper right", fontsize=9, ncol=1)
st.pyplot(fig)


