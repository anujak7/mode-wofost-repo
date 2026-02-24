"""
 WOFOST Crop Simulation - 
Date: Jan to June 2024
"""

from pcse.models import Wofost72_WLP_FD
from pcse.input import (
    CSVWeatherDataProvider,
    YAMLAgroManagementReader,
    WOFOST72SiteDataProvider
)
from pcse.base import ParameterProvider
import argparse
import csv
from datetime import datetime
import os
import sys
import yaml

print("="*60)
print(" WOFOST CROP SIMULATION")
print("="*60)
print()

# ============================================
# File Paths (Absolute for reliability)
# ============================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RUNTIME_DIR = os.path.join(BASE_DIR, "_runtime")
os.makedirs(RUNTIME_DIR, exist_ok=True)

def _resolve_path(path_value):
    if os.path.isabs(path_value):
        return path_value
    return os.path.join(BASE_DIR, path_value)


def _select_weather_file(weather_arg):
    if not weather_arg:
        raise RuntimeError("Weather file path is required.")
    return _resolve_path(weather_arg)


def _parse_weather_day(raw_value):
    value = raw_value.strip()
    for fmt in ("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", "%d-%m-%Y", "%Y%m%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            continue
    raise ValueError(f"Unsupported weather day format: '{raw_value}'")


def _to_literal_string(key, value):
    value = value.strip().strip(",")
    string_keys = {"Country", "Station", "Description", "Source", "Contact"}
    bool_keys = {"HasSunshine"}
    if key in string_keys:
        if value.startswith("'") and value.endswith("'"):
            return value
        return repr(value.strip("'").strip('"'))
    if key in bool_keys:
        low = value.lower()
        return "True" if low in {"true", "1", "yes"} else "False"
    try:
        return str(float(value))
    except ValueError:
        return value


def _prepare_weather_for_pcse(input_weather_file):
    with open(input_weather_file, "r", encoding="utf-8") as fp:
        lines = fp.read().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("DAY,"):
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError("Could not find DAY header in weather file.")

    metadata = {}
    for raw in lines[:header_idx]:
        line = raw.strip().strip(",")
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        metadata[key.strip()] = value.strip()

    required_metadata = [
        "Country", "Station", "Description", "Source", "Contact",
        "Longitude", "Latitude", "Elevation", "AngstromA", "AngstromB", "HasSunshine"
    ]
    merged_metadata = {}
    missing_metadata = [key for key in required_metadata if key not in metadata]
    if missing_metadata:
        raise RuntimeError(
            "Missing required weather metadata keys in input file: "
            + ", ".join(missing_metadata)
        )
    for key in required_metadata:
        merged_metadata[key] = _to_literal_string(key, metadata[key])

    weather_rows = []
    start_day = None
    end_day = None
    max_vap = None
    for row_idx, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 8:
            raise RuntimeError(
                f"Invalid weather row at line {row_idx}: expected at least 8 columns, got {len(parts)}"
            )
        day = _parse_weather_day(parts[0])
        values = [float(v) for v in parts[1:8]]
        irrad, tmin, tmax, vap, wind, rain, snowdepth = values
        if max_vap is None or vap > max_vap:
            max_vap = vap
        weather_rows.append(
            (
                day.strftime("%Y-%m-%d"),
                irrad,
                tmin,
                tmax,
                vap,
                wind,
                rain,
                snowdepth,
            )
        )
        if start_day is None or day < start_day:
            start_day = day
        if end_day is None or day > end_day:
            end_day = day

    if not weather_rows:
        raise RuntimeError("No weather observations found after parsing weather file.")

    # PCSE CSVWeatherDataProvider expects VAP in kPa and internally multiplies by 10.
    # If input looks like hPa-scale values (~90-110), normalize to kPa using same source data.
    if max_vap is not None and max_vap > 30.0:
        normalized_rows = []
        for day, irrad, tmin, tmax, vap, wind, rain, snowdepth in weather_rows:
            normalized_rows.append((day, irrad, tmin, tmax, vap / 10.0, wind, rain, snowdepth))
        weather_rows = normalized_rows

    prepared_weather_file = os.path.join(RUNTIME_DIR, "weather_prepared_pcse.csv")
    with open(prepared_weather_file, "w", encoding="utf-8", newline="") as fp:
        fp.write("## Station data\n")
        for key in [
            "Country", "Station", "Description", "Source", "Contact",
            "Longitude", "Latitude", "Elevation", "AngstromA", "AngstromB", "HasSunshine"
        ]:
            fp.write(f"{key} = {merged_metadata[key]}\n")
        fp.write("## Daily weather observations\n")
        writer = csv.writer(fp)
        writer.writerow(["DAY", "IRRAD", "TMIN", "TMAX", "VAP", "WIND", "RAIN", "SNOWDEPTH"])
        writer.writerows(weather_rows)

    return prepared_weather_file, start_day, end_day, len(weather_rows)


def _unwrap_pcse_value(raw_value):
    # Official PCSE YAML stores parameters as [value, description, units]
    if isinstance(raw_value, list) and len(raw_value) >= 2 and isinstance(raw_value[1], str):
        return raw_value[0]
    return raw_value


def _load_crop_from_local_file(crop_yaml_file):
    with open(crop_yaml_file, "r", encoding="utf-8") as fp:
        crop_doc = yaml.safe_load(fp) or {}

    crop_parameters = crop_doc.get("CropParameters", {})
    if not isinstance(crop_parameters, dict):
        raise RuntimeError("Missing 'CropParameters' block in local crop file.")

    selected_name = "winter_wheat"
    raw_parameters = None

    ecotypes = crop_parameters.get("EcoTypes")
    if isinstance(ecotypes, dict) and len(ecotypes) > 0:
        raw_parameters = ecotypes.get(selected_name)
        if raw_parameters is None and len(ecotypes) == 1:
            selected_name, raw_parameters = next(iter(ecotypes.items()))

    if raw_parameters is None:
        varieties = crop_parameters.get("Varieties")
        if isinstance(varieties, dict) and len(varieties) > 0:
            raw_parameters = varieties.get(selected_name)
            if raw_parameters is None and len(varieties) == 1:
                selected_name, raw_parameters = next(iter(varieties.items()))

    if raw_parameters is None or not isinstance(raw_parameters, dict):
        raise RuntimeError(
            "Could not find crop parameters for 'winter_wheat' in local crop file "
            "(expected CropParameters -> EcoTypes/Varieties)."
        )

    cropdata = {
        key: _unwrap_pcse_value(value)
        for key, value in raw_parameters.items()
        if key != "Metadata"
    }
    return cropdata, selected_name

parser = argparse.ArgumentParser(description="Run WOFOST simulation with selected weather and agromanagement files.")
parser.add_argument(
    "--weather",
    default=os.path.join("weather", "weather_pcse_2014_2034.csv"),
    help="Weather CSV path (absolute or relative). Only this file is used for weather input.",
)
parser.add_argument("--agro", default=os.path.join("agromanagement", "agro.yaml"), help="Agromanagement YAML path (absolute or relative).")
parser.add_argument(
    "--crop-file",
    default=os.path.join("crop", "winter_wheat.yaml"),
    help="Local crop YAML path. Only this file is used for crop parameters.",
)
parser.add_argument("--soil-file", default=os.path.join("soil", "soil.yaml"), help="Local soil YAML file.")
args = parser.parse_args()

WEATHER_FILE = _select_weather_file(args.weather)
AGRO_FILE = _resolve_path(args.agro)
CROP_FILE = _resolve_path(args.crop_file)
SOIL_FILE = _resolve_path(args.soil_file)

# ============================================
# Validate Files Exist
# ============================================
print(" Checking required files...\n")

missing_files = []
if not os.path.exists(WEATHER_FILE):
    missing_files.append(f"    {WEATHER_FILE}")
if not os.path.exists(AGRO_FILE):
    missing_files.append(f"    {AGRO_FILE}")
if not os.path.exists(CROP_FILE):
    missing_files.append(f"    {CROP_FILE} (crop file)")
if not os.path.exists(SOIL_FILE):
    missing_files.append(f"    {SOIL_FILE} (soil file)")

if missing_files:
    print(" ERROR: Missing files!\n")
    print("\n".join(missing_files))
    print("\n Make sure all files are in the correct folders!")
    sys.exit(1)

print(f"    Weather: {WEATHER_FILE}")
print(f"    Agromanagement: {AGRO_FILE}")
print(f"    Crop file: {CROP_FILE}")
print(f"    Soil file: {SOIL_FILE}")
print()

# ============================================
# Load Data Providers
# ============================================

# 1. Weather
print(" Loading weather data...")
try:
    prepared_weather_file, weather_start_day, weather_end_day, weather_rows = _prepare_weather_for_pcse(WEATHER_FILE)
    weather = CSVWeatherDataProvider(prepared_weather_file, dateformat="%Y-%m-%d", force_reload=True)
    print(f"    Weather rows: {weather_rows} ({weather_start_day} to {weather_end_day})")
    print(f"    Prepared weather file: {prepared_weather_file}")
    print("    Weather data loaded successfully!")
except Exception as e:
    print(f"    Error loading weather: {e}")
    print("\n Make sure weather file has valid PCSE columns: DAY,IRRAD,TMIN,TMAX,VAP,WIND,RAIN,SNOWDEPTH")
    sys.exit(1)

# 2. Agromanagement
print("  Loading agromanagement...")
campaign_start_date = None
crop_start_date_from_agro = None
crop_end_date_from_agro = None
crop_end_type_from_agro = None
try:
    agro = YAMLAgroManagementReader(AGRO_FILE)
    for campaign_item in agro:
        for campaign_date, campaign_data in campaign_item.items():
            if campaign_start_date is None:
                campaign_start_date = campaign_date
            crop_calendar = campaign_data.get("CropCalendar", {}) if campaign_data else {}
            crop_start_date = crop_calendar.get("crop_start_date")
            crop_end_date = crop_calendar.get("crop_end_date")
            if crop_start_date_from_agro is None:
                crop_start_date_from_agro = crop_start_date
                crop_end_date_from_agro = crop_end_date
                crop_end_type_from_agro = crop_calendar.get("crop_end_type")
            if crop_start_date and crop_start_date < weather_start_day:
                print(
                    f"    Warning: crop_start_date {crop_start_date} is before weather start {weather_start_day}. "
                    "Model may fail or truncate behavior."
                )
            if crop_end_date and crop_end_date > weather_end_day:
                print(
                    f"    Warning: crop_end_date {crop_end_date} is after weather end {weather_end_day}. "
                    "Model may fail or truncate behavior."
                )
            if campaign_date < weather_start_day:
                print(
                    f"    Warning: campaign start {campaign_date} is before weather start {weather_start_day}."
                )
    print("    Agromanagement loaded successfully!")
except Exception as e:
    print(f"    Error loading agromanagement: {e}")
    sys.exit(1)

# 3. Crop
print(" Loading crop parameters...")
try:
    crop, selected_profile = _load_crop_from_local_file(CROP_FILE)
    print(f"    Crop source: local file ({CROP_FILE})")
    print(f"    Selected profile: {selected_profile}")
    print(f"    Loaded crop parameter keys: {len(crop)}")
    print("    Crop parameters loaded successfully!")
except Exception as e:
    print(f"    Error loading crop: {e}")
    print("    Tip: update crop/winter_wheat.yaml with all required WOFOST parameters.")
    sys.exit(1)

REQUIRED_SITE_KEYS = ("SSMAX", "WAV", "NOTINF", "SSI", "SMLIM")
REQUIRED_SOIL_KEYS = ("CRAIRC", "K0", "KSUB", "RDMSOL", "SM0", "SMFCF", "SMW", "SOPE")


def _read_required_float(name, source_nodes, source_label):
    for node in source_nodes:
        if isinstance(node, dict) and name in node:
            try:
                return float(node[name])
            except (TypeError, ValueError):
                raise RuntimeError(f"'{name}' in {source_label} is not a valid number: {node[name]}")
    return None


def _read_weighted_layer_value(layers, name, source_label):
    if not isinstance(layers, list) or not layers:
        return None

    weighted_sum = 0.0
    total_thickness = 0.0
    previous_depth = 0.0
    for layer in layers:
        if not isinstance(layer, dict):
            continue
        if "Depth" not in layer:
            continue
        try:
            current_depth = float(layer["Depth"])
        except (TypeError, ValueError):
            raise RuntimeError(f"'Depth' in {source_label} is not a valid number: {layer['Depth']}")
        thickness = current_depth - previous_depth
        previous_depth = current_depth
        if thickness <= 0:
            continue
        if name not in layer:
            continue
        try:
            value = float(layer[name])
        except (TypeError, ValueError):
            raise RuntimeError(f"'{name}' in {source_label} is not a valid number: {layer[name]}")
        weighted_sum += value * thickness
        total_thickness += thickness

    if total_thickness > 0:
        return weighted_sum / total_thickness

    simple_values = []
    for layer in layers:
        if isinstance(layer, dict) and name in layer:
            try:
                simple_values.append(float(layer[name]))
            except (TypeError, ValueError):
                raise RuntimeError(f"'{name}' in {source_label} is not a valid number: {layer[name]}")
    if simple_values:
        return sum(simple_values) / len(simple_values)

    return None


def _load_site_and_soil_from_yaml(soil_file):
    with open(soil_file, "r", encoding="utf-8") as fp:
        soil_doc = yaml.safe_load(fp) or {}

    soil_params = soil_doc.get("SoilParameters")
    if not isinstance(soil_params, dict):
        raise RuntimeError("Missing 'SoilParameters' mapping in soil YAML.")

    site_node = soil_params.get("Site")
    if not isinstance(site_node, dict):
        raise RuntimeError("Missing 'SoilParameters -> Site' mapping in soil YAML.")

    soil_node = soil_params.get("SoilPhysical", soil_params.get("Soil"))
    if not isinstance(soil_node, dict):
        raise RuntimeError("Missing 'SoilParameters -> SoilPhysical' mapping in soil YAML.")
    layers = soil_params.get("Layers", [])

    site_values = {
        key: _read_required_float(key, [site_node], "SoilParameters -> Site")
        for key in REQUIRED_SITE_KEYS
    }

    if any(value is None for value in site_values.values()):
        missing_site = [k for k, v in site_values.items() if v is None]
        raise RuntimeError("Missing required Site keys in soil YAML: " + ", ".join(missing_site))

    soil_values = {}
    for key in REQUIRED_SOIL_KEYS:
        value = _read_required_float(
            key,
            [soil_node, soil_params],
            "SoilParameters -> SoilPhysical/SoilParameters",
        )
        if value is None and key in {"K0", "SM0", "SMFCF", "SMW"}:
            value = _read_weighted_layer_value(layers, key, "SoilParameters -> Layers")
        if value is None:
            raise RuntimeError(
                f"Missing required key '{key}' in SoilParameters -> SoilPhysical/SoilParameters/Layers"
            )
        soil_values[key] = value

    return site_values, soil_values


# 4. Site + Soil data
print(" Creating site and soil parameters...")
try:
    site_values, soil = _load_site_and_soil_from_yaml(SOIL_FILE)
    site = WOFOST72SiteDataProvider(**site_values)
    print(f"    Soil source: local soil YAML ({SOIL_FILE})")
    print("    Loaded site keys: " + ", ".join(sorted(site_values.keys())))
    print("    Loaded soil keys: " + ", ".join(sorted(soil.keys())))

    print("    Site parameters created successfully!")
    print("    Soil parameters created successfully!")
except Exception as e:
    print(f"    Error creating site/soil parameters: {e}")
    sys.exit(1)

# ============================================
# Run WOFOST Model
# ============================================
print("\n" + "="*60)
print(" STARTING WOFOST SIMULATION")
print("="*60)
print()

try:
    # Initialize model
    print(" Initializing model...")
    parameters = ParameterProvider(sitedata=site, soildata=soil, cropdata=crop)
    model = Wofost72_WLP_FD(
        parameterprovider=parameters,
        weatherdataprovider=weather,
        agromanagement=agro
    )
    print("    Model initialized!")
    
    # Run simulation
    print(" Running simulation (this may take a few seconds)...")
    model.run_till_terminate()
    print("    Simulation completed!\n")
    
except Exception as e:
    print(f"\n Simulation failed: {e}\n")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================
# Display Results
# ============================================

# Get output
output = model.get_output()
final = output[-1]

def _nz(value, default=0.0):
    return default if value is None else value

max_dvs = max((_nz(d.get("DVS"), 0.0) for d in output), default=0.0)
total_biomass = _nz(final.get("TAGP"), 0.0) + _nz(final.get("TWRT"), 0.0)

print("="*60)
print(" SIMULATION RESULTS")
print("="*60)
print(f"\n FINAL OUTPUTS (Simulation Day {final['day']}):\n")

# Main metrcs
print("    YIELD & BIOMASS:")
print(f"      • Crop Yield (TWSO):           {_nz(final.get('TWSO'), 0.0):.2f} kg/ha")
print(f"      • Total Above-Ground (TAGP):   {_nz(final.get('TAGP'), 0.0):.2f} kg/ha")
print(f"      • Total Biomass:               {total_biomass:.2f} kg/ha")

print(f"\n    GROWTH INDICATORS:")
print(f"      • Leaf Area Index (LAI):       {_nz(final.get('LAI'), 0.0):.3f}")
print(f"      • Development Stage (DVS):     {_nz(final.get('DVS'), 0.0):.2f}")
print(f"      • Leaf Weight (TWLV):          {_nz(final.get('TWLV'), 0.0):.2f} kg/ha")
print(f"      • Stem Weight (TWST):          {_nz(final.get('TWST'), 0.0):.2f} kg/ha")
print(f"      • Root Weight (TWRT):          {_nz(final.get('TWRT'), 0.0):.2f} kg/ha")

print(f"\n    WATER STATUS:")
print(f"      • Soil Moisture (SM):          {_nz(final.get('SM'), 0.0):.2f} cm")
print(f"      • Transpiration (TRA):         {_nz(final.get('TRA'), 0.0):.3f} cm/day")

print(f"\n    SIMULATION INFO:")
print(f"      • Total Days Simulated:        {len(output)}")
if campaign_start_date is not None:
    print(f"      • Campaign Start (Agro):       {campaign_start_date}")
if crop_start_date_from_agro is not None:
    print(f"      • Crop Start (Agro):           {crop_start_date_from_agro}")
if crop_end_date_from_agro is not None:
    print(f"      • Crop End Date (Agro):        {crop_end_date_from_agro}")
if crop_end_type_from_agro is not None:
    print(f"      • Crop End Type (Agro):        {crop_end_type_from_agro}")
print(f"      • Final Date:                  {final['day']}")
print(f"      • Max DVS reached:             {max_dvs:.2f}")

if max_dvs < 1.0:
    print("\n    NOTE:")
    print("      • Crop did not reach reproductive stage (DVS < 1.0).")
    print("      • Low/zero grain yield can occur due to weather-crop mismatch.")

print("\n" + "="*60)

# ============================================
# Save to CSV (Optional)
# ============================================
try:
    import pandas as pd
    df = pd.DataFrame(output)
    output_csv = os.path.join(BASE_DIR, "simulation_output.csv")
    df.to_csv(output_csv, index=False)
    print(f" Results saved to: simulation_output.csv")
    print(f"   (You can open this in Excel for analysis)")

    # Optional crop-period-only output for easier date validation
    if crop_start_date_from_agro is not None and "day" in df.columns:
        crop_period_df = df[pd.to_datetime(df["day"]).dt.date >= crop_start_date_from_agro]
        crop_period_csv = os.path.join(BASE_DIR, "simulation_output_crop_period.csv")
        crop_period_df.to_csv(crop_period_csv, index=False)
        print(" Crop-period output saved to: simulation_output_crop_period.csv")

    # Save only the final simulation day in separate files
    final_df = pd.DataFrame([final])
    final_df.to_csv(os.path.join(BASE_DIR, "simulation_final_output.csv"), index=False)
    final_df.to_json(
        os.path.join(BASE_DIR, "simulation_final_output.json"),
        orient="records",
        indent=2,
        date_format="iso"
    )
    print(" Final day summary saved to:")
    print("   • simulation_final_output.csv")
    print("   • simulation_final_output.json")

    # Create matplotlib dashboard
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch

        df["day"] = pd.to_datetime(df["day"])
        fig, axes = plt.subplots(3, 2, figsize=(15, 11))
        fig.suptitle("WOFOST Simulation Dashboard", fontsize=16, fontweight="bold")

        # 1) Development and canopy
        axes[0, 0].plot(df["day"], df["DVS"], label="DVS", color="#1f77b4")
        axes[0, 0].plot(df["day"], df["LAI"], label="LAI", color="#2ca02c")
        axes[0, 0].set_title("Development and Canopy")
        axes[0, 0].legend()
        axes[0, 0].grid(alpha=0.3)

        # 2) Production
        axes[0, 1].plot(df["day"], df["TAGP"], label="TAGP", color="#ff7f0e")
        axes[0, 1].plot(df["day"], df["TWSO"], label="TWSO", color="#d62728")
        axes[0, 1].set_title("Biomass and Yield (kg/ha)")
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3)

        # 3) Organ weights
        axes[1, 0].plot(df["day"], df["TWLV"], label="TWLV", color="#9467bd")
        axes[1, 0].plot(df["day"], df["TWST"], label="TWST", color="#8c564b")
        axes[1, 0].plot(df["day"], df["TWRT"], label="TWRT", color="#e377c2")
        axes[1, 0].set_title("Organ Weights (kg/ha)")
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)

        # 4) Water status
        axes[1, 1].plot(df["day"], df["SM"], label="SM", color="#17becf")
        axes[1, 1].plot(df["day"], df["TRA"], label="TRA", color="#bcbd22")
        axes[1, 1].set_title("Water Status")
        axes[1, 1].legend()
        axes[1, 1].grid(alpha=0.3)

        # 5) Rooting depth
        axes[2, 0].plot(df["day"], df["RD"], color="#1f77b4")
        axes[2, 0].set_title("Rooting Depth (cm)")
        axes[2, 0].grid(alpha=0.3)

        # 6) Water stress reduction factor
        axes[2, 1].plot(df["day"], df["RFTRA"], color="#d62728")
        axes[2, 1].set_ylim(0, 1.05)
        axes[2, 1].set_title("Transpiration Reduction Factor (RFTRA)")
        axes[2, 1].grid(alpha=0.3)

        for ax in axes.flat:
            ax.tick_params(axis="x", rotation=20)

        fig.text(
            0.01,
            0.01,
            f"Final Day: {final['day']} | Final Yield TWSO: {final['TWSO']:.2f} kg/ha | Final TAGP: {final.get('TAGP', 0):.2f} kg/ha",
            fontsize=10
        )

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        dashboard_png = os.path.join(BASE_DIR, "simulation_dashboard.png")
        fig.savefig(dashboard_png, dpi=150)
        print(" Dashboard saved to: simulation_dashboard.png")

        # Additional key metrics trend chart
        metrics_fig, metrics_ax = plt.subplots(figsize=(14, 5))
        key_series = [
            ("TAGP", "#ff7f0e"),
            ("TWSO", "#d62728"),
            ("LAI", "#2ca02c"),
            ("DVS", "#1f77b4"),
        ]
        for series_name, color in key_series:
            if series_name in df.columns:
                metrics_ax.plot(df["day"], df[series_name], label=series_name, color=color, linewidth=2)

        metrics_ax.set_title("Key Simulation Metrics Over Time")
        metrics_ax.set_xlabel("Date")
        metrics_ax.grid(alpha=0.3)
        metrics_ax.legend()
        metrics_ax.tick_params(axis="x", rotation=20)
        metrics_fig.tight_layout()
        metrics_png = os.path.join(BASE_DIR, "simulation_key_metrics.png")
        metrics_fig.savefig(metrics_png, dpi=150)
        print(" Key metrics plot saved to: simulation_key_metrics.png")

        # Input -> Processing -> Output flowchart image
        flow_fig, flow_ax = plt.subplots(figsize=(12, 3.8))
        flow_ax.set_axis_off()
        flow_ax.set_xlim(0, 12)
        flow_ax.set_ylim(0, 4)

        boxes = [
            (0.5, 1.2, 3.2, 1.6, "INPUT\nweather + agro + crop + soil"),
            (4.4, 1.2, 3.2, 1.6, "PROCESSING\nnormalize + load + run WOFOST"),
            (8.3, 1.2, 3.2, 1.6, "OUTPUT\nCSV/JSON + plots + dashboard"),
        ]

        for x, y, w, h, label in boxes:
            patch = FancyBboxPatch(
                (x, y), w, h,
                boxstyle="round,pad=0.2",
                linewidth=1.5,
                edgecolor="#1f1f1f",
                facecolor="#f2f6ff",
            )
            flow_ax.add_patch(patch)
            flow_ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=10)

        arrow_props = dict(arrowstyle="->", lw=2.0, color="#1f1f1f")
        flow_ax.annotate("", xy=(4.2, 2.0), xytext=(3.7, 2.0), arrowprops=arrow_props)
        flow_ax.annotate("", xy=(8.1, 2.0), xytext=(7.6, 2.0), arrowprops=arrow_props)
        flow_ax.set_title("WOFOST Pipeline: Input -> Processing -> Output", fontsize=12, fontweight="bold")
        flow_fig.tight_layout()
        flow_png = os.path.join(BASE_DIR, "simulation_pipeline_flowchart.png")
        flow_fig.savefig(flow_png, dpi=150)
        print(" Flowchart saved to: simulation_pipeline_flowchart.png")

        # Optional live dashboard window (off by default to avoid blocking runs)
        if os.environ.get("WOFOST_SHOW_DASHBOARD", "0") == "1":
            plt.show()
        else:
            plt.close(fig)
            plt.close(metrics_fig)
            plt.close(flow_fig)
    except ImportError:
        print("  matplotlib not installed - skipping dashboard plot")
        print("   To enable dashboard: pip install matplotlib")
    except Exception as e:
        print(f"  Could not generate dashboard plot: {e}")
except ImportError:
    print("  pandas not installed - skipping CSV export")
    print("   To enable CSV export: pip install pandas")
except Exception as e:
    print(f"  Could not save CSV: {e}")

print("\n WOFOST simulation complete! \n")
print("="*60)
