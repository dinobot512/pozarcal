#pragma once
#include <array>
#include <string>

// ---------------------------------------------------------------------------
// FuelModel.h
// Fuel model parameters per Anderson 13 (and extensible to Scott-Burgan 40).
// Units are Rothermel-native: lb, ft, BTU — NOT SI.
// All loads are lb/ft² (converted from tons/acre at load time: 1 ton/acre = 0.04591 lb/ft²).
// SAV ratios are ft²/ft³ (= 1/ft).
// ---------------------------------------------------------------------------

struct FuelModel {
    int    id;               // Anderson 1-13, or Scott-Burgan 101-189
    std::string name;

    // Fuel loads [lb/ft²]
    float  w_1;              // 1-hr dead
    float  w_10;             // 10-hr dead
    float  w_100;            // 100-hr dead
    float  w_lh;             // live herbaceous
    float  w_lw;             // live woody

    // Surface-area-to-volume ratios [ft²/ft³ = 1/ft]
    float  sigma_1;          // 1-hr dead SAV
    float  sigma_10;         // 10-hr dead SAV  (fixed at 109 for all Anderson models)
    float  sigma_100;        // 100-hr dead SAV (fixed at 30  for all Anderson models)
    float  sigma_lh;         // live herbaceous SAV
    float  sigma_lw;         // live woody SAV

    float  delta;            // fuel bed depth [ft]
    float  M_x;              // dead fuel moisture of extinction [fraction, e.g. 0.25]
    float  h;                // heat content [BTU/lb] — typically 8000 dead, 8000 live
    float  rho_p;            // particle density [lb/ft³] — typically 32
};

// Conversion: tons/acre → lb/ft²
static constexpr float TONS_ACRE_TO_LB_FT2 = 0.04591f;

// ---------------------------------------------------------------------------
// Anderson 13 standard fuel models.
// Source: Albini 1976 / Andrews 1986 tabulated values.
// Loads stored pre-converted to lb/ft².
// ---------------------------------------------------------------------------
static const std::array<FuelModel, 13> ANDERSON_13 = {{
    // id  name            w_1      w_10     w_100    w_lh     w_lw    σ_1    σ_10  σ_100  σ_lh   σ_lw   δ(ft)  Mx     h      ρp
    {  1, "Short grass",  0.0459,  0.0000,  0.0000,  0.0000,  0.0000, 3500,  109,   30,    1500,  1500,  1.0f,  0.12f, 8000,  32 },
    {  2, "Timber grass", 0.0918,  0.0459,  0.0459,  0.0459,  0.0000, 3000,  109,   30,    1500,  1500,  1.0f,  0.15f, 8000,  32 },
    {  3, "Tall grass",   0.1377,  0.0000,  0.0000,  0.0000,  0.0000, 1500,  109,   30,    1500,  1500,  2.5f,  0.25f, 8000,  32 },
    {  4, "Chaparral",    0.2296,  0.1838,  0.1607,  0.0000,  0.2525, 2000,  109,   30,    1500,  1500,  6.0f,  0.20f, 8000,  32 },
    {  5, "Brush",        0.0459,  0.0230,  0.0000,  0.0000,  0.0918, 2000,  109,   30,    1500,  1500,  2.0f,  0.20f, 8000,  32 },
    {  6, "Dormant brush",0.0688,  0.0459,  0.0459,  0.0000,  0.0000, 1750,  109,   30,    1500,  1500,  2.5f,  0.25f, 8000,  32 },
    {  7, "Southern rough",0.0528, 0.0183,  0.0183,  0.0000,  0.0459, 1750,  109,   30,    1500,  1500,  2.5f,  0.40f, 8000,  32 },
    {  8, "Closed timber",0.0688,  0.0459,  0.1147,  0.0000,  0.0000, 2000,  109,   30,    1500,  1500,  0.2f,  0.30f, 8000,  32 },
    {  9, "Hardwood litter",0.1341,0.0194,  0.0077,  0.0000,  0.0000, 2500,  109,   30,    1500,  1500,  0.2f,  0.25f, 8000,  32 },
    { 10, "Timber litter",0.1377,  0.0918,  0.1377,  0.0000,  0.0459, 2000,  109,   30,    1500,  1500,  1.0f,  0.25f, 8000,  32 },
    { 11, "Light slash",  0.0688,  0.2066,  0.2525,  0.0000,  0.0000, 1500,  109,   30,    1500,  1500,  1.0f,  0.15f, 8000,  32 },
    { 12, "Medium slash", 0.1837,  0.5510,  0.5969,  0.0000,  0.0000, 1500,  109,   30,    1500,  1500,  2.3f,  0.20f, 8000,  32 },
    { 13, "Heavy slash",  0.3214,  0.9603,  1.3776,  0.0000,  0.0000, 1500,  109,   30,    1500,  1500,  3.0f,  0.25f, 8000,  32 },
}};

// Lookup by Anderson ID (1-13). Returns nullptr if not found.
inline const FuelModel* getFuelModel(int id) {
    if (id < 1 || id > 13) return nullptr;
    return &ANDERSON_13[id - 1];
}
