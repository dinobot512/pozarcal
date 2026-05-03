#pragma once
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <functional>
#include "Rothermel.h"

// ---------------------------------------------------------------------------
// FireGrid.h
// Cellular automaton fire propagation on a uniform 2D raster.
//
// Propagation model:
//   Each burning cell ignites its 8 neighbors. The time for a burning cell
//   (r_src, c_src) to ignite neighbor (r_dst, c_dst) is:
//
//       t_ignite = t_src + dist / R(θ)
//
//   where dist = cell_size_ft * geometric distance (1.0 or √2 for diagonal),
//   and R(θ) is the ROS of the Huygens fire ellipse evaluated in the
//   direction from src → dst.
//
// Huygens ellipse ROS at angle θ from head direction:
//   Uses the standard ellipse formula with semi-axes derived from
//   R_head, R_back, R_flank (Anderson 1983).
//
//       a = (R_head + R_back) / 2   — semi-major axis (along head)
//       b = R_flank                 — semi-minor axis
//       c = a - R_back              — focus offset from center
//
//   ROS(θ) = b² / (a - c·cos(θ))    (focal form of ellipse, θ from head dir)
//
// All cells share the same fuel model, wind, and terrain (uniform grid).
// ---------------------------------------------------------------------------

struct GridCell {
    float burn_time = -1.f;  // sim time [min] when cell ignited; -1 = unburned
};

class FireGrid {
public:
    int   rows, cols;
    float cell_size_ft;   // physical size of each cell [ft]

    std::vector<GridCell> cells;

    // Uniform environment across grid
    RothermelInput env;
    RothermelOutput ros_cache;  // computed once in constructor

    // Simulation time of last step
    float sim_time = 0.f;

    // ---------------------------------------------------------------------------
    // Constructor
    // env_in:       uniform Rothermel input (fuel, moisture, wind, terrain)
    // cell_size_ft: resolution of grid [ft]
    // ---------------------------------------------------------------------------
    FireGrid(int rows_, int cols_, float cell_size_ft_, const RothermelInput& env_in)
        : rows(rows_), cols(cols_), cell_size_ft(cell_size_ft_), env(env_in)
    {
        cells.resize(rows * cols);
        ros_cache = Rothermel::computeROS(env);
    }

    GridCell& at(int r, int c)             { return cells[r * cols + c]; }
    const GridCell& at(int r, int c) const { return cells[r * cols + c]; }

    bool inBounds(int r, int c) const {
        return r >= 0 && r < rows && c >= 0 && c < cols;
    }

    // ---------------------------------------------------------------------------
    // Huygens ellipse ROS [ft/min] in direction θ from head direction.
    //
    //   a = (R_head + R_back) / 2
    //   b = R_flank
    //   c = a - R_back
    //   R(θ) = b² / (a - c·cos(θ))
    //
    // θ is the angle at the source cell between the head direction and the
    // direction toward the target cell (neighbor).
    // ---------------------------------------------------------------------------
    // Anderson (1983) focal form — ignition at REAR focus.
    // r(θ) = p / (1 - e·cos(θ))
    //   e = (R_head - R_back) / (R_head + R_back)
    //   a = R_head / (1 + e)
    //   p = a·(1 - e²)
    // Guarantees r(0)=R_head, r(π)=R_back exactly.
    float ellipseROS(float theta_rad) const {
        float R_h = ros_cache.ROS_head;
        float R_b = ros_cache.ROS_back;

        if (R_h <= 0.f) return 0.f;
        if (R_b <= 0.f) R_b = 1e-6f;

        float e = (R_h - R_b) / (R_h + R_b);
        float a = R_h / (1.f + e);
        float p = a * (1.f - e * e);

        float denom = 1.f - e * std::cos(theta_rad);
        if (denom <= 0.f) return 0.f;
        return p / denom;
    }

    // ---------------------------------------------------------------------------
    // Ignite a single cell at time t_ignite.
    // Returns false if already burned.
    // ---------------------------------------------------------------------------
    bool ignite(int r, int c, float t_ignite) {
        if (!inBounds(r, c)) return false;
        GridCell& cell = at(r, c);
        if (cell.burn_time >= 0.f) return false;
        cell.burn_time = t_ignite;
        return true;
    }

    // ---------------------------------------------------------------------------
    // run() — propagate fire using a minimum-time priority queue (Dijkstra-like).
    //
    // Starts from all already-ignited cells. Expands until no cells remain
    // or max_time_min is reached.
    //
    // This is equivalent to solving the eikonal equation on the grid with
    // anisotropic speed (the Huygens ellipse).
    // ---------------------------------------------------------------------------
    void run(float max_time_min = std::numeric_limits<float>::infinity()) {
        // Priority queue: (burn_time, row, col) — min-heap
        using Entry = std::tuple<float, int, int>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

        // Seed with all currently ignited cells
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
                if (cells[r * cols + c].burn_time >= 0.f)
                    pq.push({cells[r * cols + c].burn_time, r, c});

        // Head direction vector (unit, in grid East/North space)
        // Grid convention: col increases East, row increases South.
        // head_dir_deg: 0=N, 90=E, 180=S, 270=W (met, fire travels toward)
        static constexpr float DEG2RAD = 3.14159265f / 180.f;
        float head_rad = (90.f - ros_cache.head_dir_deg) * DEG2RAD;
        float hx = std::cos(head_rad);   // East component
        float hy = std::sin(head_rad);   // North component (grid row decreases northward)

        // 8-neighbor offsets (dr, dc)
        static const int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
        static const int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

        while (!pq.empty()) {
            auto [t_src, r_src, c_src] = pq.top(); pq.pop();

            // Skip stale entries (cell was already ignited earlier)
            if (at(r_src, c_src).burn_time < t_src - 1e-6f) continue;

            if (t_src > max_time_min) break;
            sim_time = std::max(sim_time, t_src);

            for (int k = 0; k < 8; ++k) {
                int r_dst = r_src + DR[k];
                int c_dst = c_src + DC[k];
                if (!inBounds(r_dst, c_dst)) continue;
                if (at(r_dst, c_dst).burn_time >= 0.f) continue;

                // Direction vector from src → dst in grid space
                // dc = East offset, -dr = North offset (row down = South)
                float dx =  (float)DC[k];   // East
                float dy = -(float)DR[k];   // North

                // Physical distance
                float dist_cells = std::sqrt(dx*dx + dy*dy);  // 1 or √2
                float dist_ft    = dist_cells * cell_size_ft;

                // Angle between head direction and spread direction
                // cos(θ) = dot(head_unit, spread_unit)
                float spread_len = dist_cells;   // already computed
                float dot = (hx * dx + hy * dy) / spread_len;
                dot = std::max(-1.f, std::min(1.f, dot));  // clamp for acos
                float theta = std::acos(dot);

                // ROS in this direction [ft/min]
                float ros = ellipseROS(theta);
                if (ros <= 0.f) continue;

                float t_arrive = t_src + dist_ft / ros;
                if (t_arrive > max_time_min) continue;

                // Relax: update if this path is faster
                GridCell& dst = at(r_dst, c_dst);
                if (dst.burn_time < 0.f || t_arrive < dst.burn_time) {
                    dst.burn_time = t_arrive;
                    pq.push({t_arrive, r_dst, c_dst});
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    // ASCII render at a given simulation time.
    // '.' = unburned, '*' = burning (ignited by t), '#' = burned (ignited before t-1)
    // ---------------------------------------------------------------------------
    void printASCII(float t) const {
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                float bt = at(r, c).burn_time;
                if (bt < 0.f)       std::putchar('.');   // unburned
                else if (bt > t)    std::putchar('.');   // not yet reached
                else if (bt >= t - 5.f) std::putchar('*');  // active flame front
                else                std::putchar('#');   // burned out
            }
            std::putchar('\n');
        }
    }
};
