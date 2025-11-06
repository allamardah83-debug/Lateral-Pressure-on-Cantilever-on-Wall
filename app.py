import math
import pandas as pd
import numpy as np
import plotly.express as px
import streamlit as st

st.set_page_config(page_title="Active & Passive Earth Pressure", layout="wide")

st.title("Active & Passive Earth Pressure (Rankine Theory)")

# Sidebar inputs
st.sidebar.header("Input Parameters")
gw = st.sidebar.number_input("Groundwater Elevation", value=0.0)
cut = st.sidebar.number_input("Cut Elevation (Passive starts from zero)", value=-22.0)
surcharge = st.sidebar.number_input("Surcharge (psf)", value=300.0)
gamma_w = 62.4

# Soil layer input table
st.sidebar.write("### Soil Layers")
default_data = pd.DataFrame({
    "Soil": ["Fill", "Silty GRVL", "Fat CLAY"],
    "Top Elev.": [0, -34, -62],
    "Bottom Elev.": [-34, -62, -70],
    "Density": [120, 122, 123],
    "Friction Angle": [30, 34, 27],
    "Cohesion": [0, 0, 0]
})
data = st.data_editor(default_data, use_container_width=True)

# Helper functions
def Ka(phi): return math.tan(math.radians(45 - phi/2))**2
def Kp(phi): return math.tan(math.radians(45 + phi/2))**2

def effective_stress(z, layers, gw, surcharge):
    sigma = surcharge
    for _, L in layers.iterrows():
        if z <= L["Top Elev."] and z >= L["Bottom Elev."]:
            sigma += (L["Density"] - (gamma_w if z < gw else 0)) * (L["Top Elev."] - z)
            break
        elif z < L["Bottom Elev."]:
            sigma += (L["Density"] - (gamma_w if L["Bottom Elev."] < gw else 0)) * (L["Top Elev."] - L["Bottom Elev."])
    return sigma

# Compute pressures
rows = []
sigma_cut = effective_stress(cut, data, gw, surcharge)
for _, L in data.iterrows():
    ka = Ka(L["Friction Angle"])
    kp = Kp(L["Friction Angle"])
    sig_t = effective_stress(L["Top Elev."], data, gw, surcharge)
    sig_b = effective_stress(L["Bottom Elev."], data, gw, surcharge)
    passive_t = 0 if L["Top Elev."] > cut else kp * max(0, sig_t - sigma_cut) + 2 * math.sqrt(kp) * L["Cohesion"]
    passive_b = 0 if L["Bottom Elev."] > cut else kp * max(0, sig_b - sigma_cut) + 2 * math.sqrt(kp) * L["Cohesion"]
    active_t = ka * sig_t - 2 * math.sqrt(ka) * L["Cohesion"]
    active_b = ka * sig_b - 2 * math.sqrt(ka) * L["Cohesion"]
    rows.append([L["Soil"], L["Top Elev."], L["Bottom Elev."], ka, kp, passive_t, passive_b, active_t, active_b])

df = pd.DataFrame(rows, columns=["Soil","Top","Bottom","Ka","Kp","Passive Top","Passive Bottom","Active Top","Active Bottom"])
st.write("### Results", df)

# Plot
fig = px.line(df.melt(id_vars=["Soil","Top","Bottom"], value_vars=["Passive Top","Active Top"],
                      var_name="Type", value_name="Pressure (psf)"),
              x="Pressure (psf)", y="Top", color="Type", title="Pressure Profile")
fig.update_yaxes(autorange="reversed")_
