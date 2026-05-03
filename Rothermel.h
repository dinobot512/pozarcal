#pragma once
#include <cmath>
#include "RothermelInput.h"

// ---------------------------------------------------------------------------
// Rothermel.h
// Implementation of Rothermel (1972) surface fire spread model.
// All intermediate functions are exposed for unit testing.
// Units throughout: lb, ft, BTU, min — Rothermel-native.
//
// Key reference:
//   Rothermel, R.C. 1972. A mathematical model for predicting fire spread
//   in wildland fuels. USDA Forest Service Research Paper INT-115.
// ---------------------------------------------------------------------------

namespace Rothermel {

// ---------------------------------------------------------------------------
// STEP 1 — Characteristic SAV ratio (σ) and bulk density (ρ_b)
// ---------------------------------------------------------------------------

struct FuelArrays {
    static constexpr int N = 5;
    float w[N];       // load        [lb/ft²]
    float sigma[N];   // SAV ratio   [1/ft]
    float M[N];       // moisture    [fraction]
    float h[N];       // heat content [BTU/lb]
    float rho_p[N];   // particle density [lb/ft³]
    bool  is_live[N];
};

inline FuelArrays buildArrays(const RothermelInput& in) {
    const FuelModel& f  = *in.fuel;
    const MoistureInputs& m = in.moisture;
    FuelArrays a;

    a.w[0] = f.w_1;   a.sigma[0] = f.sigma_1;   a.M[0] = m.M_1;   a.is_live[0] = false;
    a.w[1] = f.w_10;  a.sigma[1] = f.sigma_10;  a.M[1] = m.M_10;  a.is_live[1] = false;
    a.w[2] = f.w_100; a.sigma[2] = f.sigma_100; a.M[2] = m.M_100; a.is_live[2] = false;
    a.w[3] = f.w_lh;  a.sigma[3] = f.sigma_lh;  a.M[3] = m.M_lh;  a.is_live[3] = true;
    a.w[4] = f.w_lw;  a.sigma[4] = f.sigma_lw;  a.M[4] = m.M_lw;  a.is_live[4] = true;

    for (int i = 0; i < FuelArrays::N; ++i) {
        a.h[i]     = f.h;
        a.rho_p[i] = f.rho_p;
    }
    return a;
}

// σ [ft²/ft³] — Rothermel eq. 71
inline float characteristicSAV(const FuelArrays& a) {
    float num = 0.f, den = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) {
        num += a.w[i] * a.sigma[i] * a.sigma[i];
        den += a.w[i] * a.sigma[i];
    }
    return (den > 0.f) ? num / den : 0.f;
}

// ρ_b [lb/ft³] — Rothermel eq. 40
inline float bulkDensity(const FuelArrays& a, float delta) {
    float total_w = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) total_w += a.w[i];
    return (delta > 0.f) ? total_w / delta : 0.f;
}

// β [-] — Rothermel eq. 31
inline float packingRatio(float rho_b, float rho_p) {
    return rho_b / rho_p;
}

// β_op [-] — Rothermel eq. 37
inline float optimumPackingRatio(float sigma) {
    return 3.348f * std::pow(sigma, -0.8189f);
}

// ---------------------------------------------------------------------------
// STEP 2 — Moisture damping η_M
// ---------------------------------------------------------------------------

inline float moistureDampingDead(float M, float M_x) {
    float r = (M_x > 0.f) ? std::min(M / M_x, 1.f) : 1.f;
    return 1.f - 2.59f*r + 5.11f*r*r - 3.52f*r*r*r;
}

// Live extinction moisture — Rothermel eq. 88
inline float liveExtinctionMoisture(const FuelArrays& a, float M_x_dead) {
    float dead_factor = 0.f, live_factor = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) {
        if (!a.is_live[i]) dead_factor += a.w[i] * std::exp(-138.f / a.sigma[i]);
        else               live_factor += a.w[i] * std::exp(-500.f / a.sigma[i]);
    }
    float W = (live_factor > 0.f) ? dead_factor / live_factor : 0.f;
    float M_x_live = 2.9f * W * (1.f - M_x_dead / 0.3f) - 0.226f;
    return std::max(M_x_live, M_x_dead);
}

inline float moistureDampingLive(float M_live_weighted, float M_x_live) {
    float r = (M_x_live > 0.f) ? std::min(M_live_weighted / M_x_live, 1.f) : 1.f;
    return 1.f - 2.59f*r + 5.11f*r*r - 3.52f*r*r*r;
}

// ---------------------------------------------------------------------------
// STEP 3 — Reaction intensity I_R [BTU/ft²/min]
// ---------------------------------------------------------------------------

// Γ' [1/min] — Rothermel eq. 36 and 38
inline float optimumReactionVelocity(float sigma, float beta, float beta_op) {
    float sigma15  = std::pow(sigma, 1.5f);
    float Gamma_max = sigma15 / (495.f + 0.0594f * sigma15);
    float A        = 1.f / (4.774f * std::pow(sigma, 0.1f) - 7.27f);
    float ratio    = beta / beta_op;
    return Gamma_max * std::pow(ratio, A) * std::exp(A * (1.f - ratio));
}

// η_s [-] — Rothermel eq. 30
inline float mineralDamping(float S_e = 0.01f) {
    return std::min(0.174f * std::pow(S_e, -0.19f), 1.f);
}

// w_n [lb/ft²] — Rothermel eq. 24
inline float netFuelLoad(float w, float S_T = 0.055f) {
    return w * (1.f - S_T);
}

// I_R [BTU/ft²/min] — Rothermel eq. 27
inline float reactionIntensity(const FuelArrays& a, const FuelModel& f) {
    float sigma   = characteristicSAV(a);
    float rho_b   = bulkDensity(a, f.delta);
    float beta    = packingRatio(rho_b, f.rho_p);
    float beta_op = optimumPackingRatio(sigma);
    float Gamma   = optimumReactionVelocity(sigma, beta, beta_op);
    float eta_s   = mineralDamping();

    float M_x_dead = f.M_x;
    float M_x_live = liveExtinctionMoisture(a, M_x_dead);

    float dead_wt = 0.f, dead_Mwt = 0.f;
    float live_wt = 0.f, live_Mwt = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) {
        if (a.w[i] <= 0.f) continue;
        if (!a.is_live[i]) { dead_wt += a.w[i]; dead_Mwt += a.w[i] * a.M[i]; }
        else               { live_wt += a.w[i]; live_Mwt += a.w[i] * a.M[i]; }
    }

    float eta_M_dead = (dead_wt > 0.f) ? moistureDampingDead(dead_Mwt / dead_wt, M_x_dead) : 0.f;
    float eta_M_live = (live_wt > 0.f) ? moistureDampingLive(live_Mwt / live_wt, M_x_live) : 0.f;

    float I_R = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) {
        if (a.w[i] <= 0.f) continue;
        float w_n   = netFuelLoad(a.w[i]);
        float eta_M = a.is_live[i] ? eta_M_live : eta_M_dead;
        I_R += Gamma * w_n * a.h[i] * eta_M * eta_s;
    }
    return std::max(I_R, 0.f);
}

// ---------------------------------------------------------------------------
// STEP 4 — Propagating flux ratio ξ [-]
//
// ξ = exp((0.792 + 0.681·√σ)·(β + 0.1)) / (192 + 0.2595·σ)
// Rothermel eq. 42 — fraction of I_R that propagates as a heat flux
// ahead of the flame front to preheat unburned fuel.
// ---------------------------------------------------------------------------
inline float propagatingFluxRatio(float sigma, float beta) {
    float num = std::exp((0.792f + 0.681f * std::sqrt(sigma)) * (beta + 0.1f));
    float den = 192.f + 0.2595f * sigma;
    return num / den;
}

// ---------------------------------------------------------------------------
// STEP 5 — Wind coefficient φ_w [-]
//
// φ_w = C · U^B · (β / β_op)^(-E)
//   C = 7.47 · exp(-0.133 · σ^0.55)
//   B = 0.02526 · σ^0.54
//   E = 0.715 · exp(-3.59e-4 · σ)
// U is midflame wind speed in ft/min.
// Rothermel eq. 47.
// ---------------------------------------------------------------------------
inline float windCoefficient(float sigma, float beta, float beta_op, float U_ft_min) {
    float C = 7.47f  * std::exp(-0.133f  * std::pow(sigma, 0.55f));
    float B = 0.02526f * std::pow(sigma, 0.54f);
    float E = 0.715f * std::exp(-3.59e-4f * sigma);
    return C * std::pow(U_ft_min, B) * std::pow(beta / beta_op, -E);
}

// ---------------------------------------------------------------------------
// STEP 6 — Slope coefficient φ_s [-]
//
// φ_s = 5.275 · β^(-0.3) · tan²(slope)
// tan(slope) is the tangent of the slope angle (rise/run).
// Rothermel eq. 51.
// ---------------------------------------------------------------------------
inline float slopeCoefficient(float beta, float tan_slope) {
    return 5.275f * std::pow(beta, -0.3f) * tan_slope * tan_slope;
}

// ---------------------------------------------------------------------------
// STEP 7 — No-wind / no-slope ROS  R_0 [ft/min]
//
// R_0 = (I_R · ξ) / (ρ_b · ε · Q_ig)
//
//   ε     = effective heating number = exp(-138 / σ)      (eq. 14 — dead fine fuel approx)
//   Q_ig  = heat of pre-ignition [BTU/lb]
//         = 250 + 1116·M   (dead fuel weighted moisture)  (eq. 12)
//
// Rothermel eq. 1 (without φ terms).
// ---------------------------------------------------------------------------

// Effective heating number ε [-] — Rothermel eq. 14
inline float effectiveHeatingNumber(float sigma) {
    return std::exp(-138.f / sigma);
}

// Heat of pre-ignition Q_ig [BTU/lb] — Rothermel eq. 12
inline float heatOfPreignition(float M_dead_weighted) {
    return 250.f + 1116.f * M_dead_weighted;
}

// R_0 [ft/min]
inline float noWindNoSlopeROS(float I_R, float xi, float rho_b, float epsilon, float Q_ig) {
    float denom = rho_b * epsilon * Q_ig;
    return (denom > 0.f) ? (I_R * xi) / denom : 0.f;
}

// ---------------------------------------------------------------------------
// STEP 8 — Vector ROS combination and ellipse geometry
//
// Wind and slope vectors are combined before applying to R_0.
// The combined effect vector magnitude gives the total φ (wind+slope).
//
// Head ROS:  R     = R_0 · (1 + φ_combined)
// Back ROS:  R_b   = R_0 / (1 + φ_combined)         (Anderson 1983 approx)
// Flank ROS: R_f   = (R + R_b) / (2 · LW_ratio)     (ellipse minor semi-axis speed)
//   where LW_ratio = 0.936·exp(0.2566·U_mph) + 0.461·exp(-0.1548·U_mph) − 0.397
//   (Andrews 1986 eq. for length-to-width ratio of the fire ellipse)
// ---------------------------------------------------------------------------

// Fire ellipse length-to-width ratio — Andrews 1986
// U_mph = midflame wind speed in mi/hr
inline float ellipseLWRatio(float U_mph) {
    float lw = 0.936f * std::exp(0.2566f * U_mph)
             + 0.461f * std::exp(-0.1548f * U_mph)
             - 0.397f;
    return std::max(lw, 1.0f);
}

// ---------------------------------------------------------------------------
// computeROS — master function returning RothermelOutput
//
// Wind direction convention (met): direction wind comes FROM, 0=N.
// Aspect convention: direction slope faces, 0=N.
// Head direction: fire moves INTO wind + upslope, resolved as a vector sum.
// ---------------------------------------------------------------------------
inline RothermelOutput computeROS(const RothermelInput& in) {
    const FuelModel& f = *in.fuel;
    FuelArrays a = buildArrays(in);

    // --- Derived fuel bed properties ---
    float sigma   = characteristicSAV(a);
    float rho_b   = bulkDensity(a, f.delta);
    float beta    = packingRatio(rho_b, f.rho_p);
    float beta_op = optimumPackingRatio(sigma);

    // Weighted dead moisture (for Q_ig and η_M)
    float dead_wt = 0.f, dead_Mwt = 0.f;
    for (int i = 0; i < FuelArrays::N; ++i) {
        if (!a.is_live[i] && a.w[i] > 0.f) {
            dead_wt  += a.w[i];
            dead_Mwt += a.w[i] * a.M[i];
        }
    }
    float M_dead_bar = (dead_wt > 0.f) ? dead_Mwt / dead_wt : 0.f;

    // --- Reaction intensity ---
    float I_R = reactionIntensity(a, f);

    // --- Propagating flux ---
    float xi = propagatingFluxRatio(sigma, beta);

    // --- No-wind / no-slope ROS ---
    float epsilon = effectiveHeatingNumber(sigma);
    float Q_ig    = heatOfPreignition(M_dead_bar);
    float R_0     = noWindNoSlopeROS(I_R, xi, rho_b, epsilon, Q_ig);

    // --- Wind coefficient ---
    // Convert midflame speed: mi/hr → ft/min  (1 mi/hr = 88 ft/min)
    float U_ft_min = in.wind.midflame_speed * 88.f;
    float phi_w    = windCoefficient(sigma, beta, beta_op, U_ft_min);

    // --- Slope coefficient ---
    float phi_s = slopeCoefficient(beta, in.terrain.slope_tan);

    // ---------------------------------------------------------------------------
    // Vector combination of wind and slope effects.
    //
    // Both wind and slope push fire in a direction. We resolve them as 2D vectors
    // in (East, North) space, then find the magnitude and direction of the resultant.
    //
    // Wind pushes fire DOWNWIND → fire head moves toward (wind_from + 180°).
    // Slope pushes fire UPSLOPE → fire head moves toward (aspect).
    //
    // wind_from_deg: direction wind comes FROM (met convention, 0=N, CW)
    // aspect_deg:    direction slope faces (0=N, CW)
    //
    // Unit vectors (trig uses math convention: 0=E, CCW — convert from CW/N=0):
    //   heading_rad = (90 - deg) * π/180
    // ---------------------------------------------------------------------------
    static constexpr float DEG2RAD = 3.14159265f / 180.f;

    // Wind: fire moves toward (wind_from + 180) % 360
    float wind_head_deg = std::fmod(in.wind.direction_deg + 180.f, 360.f);
    float wind_head_rad = (90.f - wind_head_deg) * DEG2RAD;
    float wx = phi_w * std::cos(wind_head_rad);
    float wy = phi_w * std::sin(wind_head_rad);

    // Slope: fire moves upslope, toward aspect direction
    float slope_head_rad = (90.f - in.terrain.aspect_deg) * DEG2RAD;
    float sx = phi_s * std::cos(slope_head_rad);
    float sy = phi_s * std::sin(slope_head_rad);

    // Combined vector
    float cx = wx + sx;
    float cy = wy + sy;
    float phi_combined = std::sqrt(cx*cx + cy*cy);

    // Direction of combined fire spread (math→met: deg = 90 - atan2(cy,cx)*180/π)
    float head_dir_deg = 0.f;
    if (phi_combined > 0.f) {
        float math_deg = std::atan2(cy, cx) * 180.f / 3.14159265f;
        head_dir_deg = std::fmod(90.f - math_deg + 360.f, 360.f);
    }

    // --- Head / back / flank ROS ---
    float R_head  = R_0 * (1.f + phi_combined);
    float R_back  = R_0 / (1.f + phi_combined);
    float lw      = ellipseLWRatio(in.wind.midflame_speed);
    float R_flank = (R_head + R_back) / (2.f * lw);

    RothermelOutput out{};
    out.ROS_head    = std::max(R_head,  0.f);
    out.ROS_back    = std::max(R_back,  0.f);
    out.ROS_flank   = std::max(R_flank, 0.f);
    out.head_dir_deg = head_dir_deg;
    out.phi_w       = phi_w;
    out.phi_s       = phi_s;
    out.I_R         = I_R;
    out.xi          = xi;

    return out;
}

} // namespace Rothermel
