#pragma once
#include "FuelModel.h"

// ---------------------------------------------------------------------------
// RothermelInput.h
// Everything needed to compute ROS for one cell.
// Weather and terrain are per-cell (interpolated from your GP rasters).
// ---------------------------------------------------------------------------

struct MoistureInputs {
    float M_1;      // 1-hr dead fuel moisture   [fraction, e.g. 0.06 = 6%]
    float M_10;     // 10-hr dead fuel moisture  [fraction]
    float M_100;    // 100-hr dead fuel moisture [fraction]
    float M_lh;     // live herbaceous moisture  [fraction]
    float M_lw;     // live woody moisture       [fraction]
};

struct WindInputs {
    float midflame_speed;   // midflame wind speed [mi/hr]  — Rothermel uses mi/hr internally
    float direction_deg;    // met convention: direction wind comes FROM [0-360°, 0=N]
};

struct TerrainInputs {
    float slope_tan;        // tan(slope angle) [dimensionless, e.g. 0.30 for 30% slope]
    float aspect_deg;       // aspect: direction the slope faces [0-360°, 0=N]
};

struct RothermelInput {
    const FuelModel*  fuel;      // pointer into lookup table (never null when used)
    MoistureInputs    moisture;
    WindInputs        wind;
    TerrainInputs     terrain;
};

// ---------------------------------------------------------------------------
// Output from the ROS calculation — one per cell per solver step.
// ---------------------------------------------------------------------------
struct RothermelOutput {
    float ROS_head;     // rate of spread in head direction [ft/min]
    float ROS_back;     // rate of spread backing into wind [ft/min]
    float ROS_flank;    // rate of spread at flanks         [ft/min]
    float head_dir_deg; // direction fire head moves toward [0-360°, 0=N]
                        // = (wind.direction_deg + 180) % 360, slope-corrected

    float phi_w;        // wind coefficient   (diagnostic)
    float phi_s;        // slope coefficient  (diagnostic)
    float I_R;          // reaction intensity [BTU/ft²/min] (diagnostic)
    float xi;           // propagating flux ratio           (diagnostic)
};
