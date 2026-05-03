#include <cstdio>
#include <cmath>
#include <string>
#include "FuelModel.h"
#include "RothermelInput.h"
#include "FireGrid.h"

// ---------------------------------------------------------------------------
// main.cpp — Fire propagation simulation driver
//
// Configures a uniform grid with:
//   - Anderson fuel model (selectable)
//   - Uniform wind and terrain
//   - Single point ignition at grid center
//
// Outputs:
//   - ASCII snapshots printed to stdout at regular intervals
//   - burn_times.csv  — one row per cell: row,col,burn_time_min
// ---------------------------------------------------------------------------

int main() {

    // -----------------------------------------------------------------------
    // 1. Simulation parameters
    // -----------------------------------------------------------------------
    constexpr int   GRID_ROWS     = 40;        // cells
    constexpr int   GRID_COLS     = 40;        // cells
    constexpr float CELL_SIZE_FT  = 164.f;     // ft per cell (164 ft ≈ 50 m)
    constexpr float MAX_TIME_MIN  = 60.f;      // simulate 60 minutes
    constexpr float SNAPSHOT_INTERVAL = 10.f; // print ASCII every N minutes
    constexpr int   FUEL_MODEL_ID = 2;         // Anderson 1-13

    // -----------------------------------------------------------------------
    // 2. Environment — moisture, wind, terrain
    // -----------------------------------------------------------------------
    MoistureInputs moisture {
        .M_1   = 0.06f,   // 6%  — dry 1-hr dead fuel
        .M_10  = 0.07f,   // 7%  — 10-hr dead
        .M_100 = 0.08f,   // 8%  — 100-hr dead
        .M_lh  = 0.60f,   // 60% — live herbaceous
        .M_lw  = 0.90f,   // 90% — live woody
    };

    WindInputs wind {
        .midflame_speed = 0.f,   // mi/hr
        .direction_deg  = 270.f,  // wind FROM west → fire spreads east
    };

    TerrainInputs terrain {
        .slope_tan  = 0.20f,   // ~11° slope (tan 11.3° ≈ 0.20)
        .aspect_deg = 180.f,   // slope faces south → upslope is northward
    };

    // -----------------------------------------------------------------------
    // 3. Build RothermelInput
    // -----------------------------------------------------------------------
    const FuelModel* fuel = getFuelModel(FUEL_MODEL_ID);
    if (!fuel) {
        std::fprintf(stderr, "Invalid fuel model ID %d\n", FUEL_MODEL_ID);
        return 1;
    }

    RothermelInput env { fuel, moisture, wind, terrain };

    // -----------------------------------------------------------------------
    // 4. Print computed ROS diagnostics
    // -----------------------------------------------------------------------
    RothermelOutput ros = Rothermel::computeROS(env);
    std::printf("=== Rothermel ROS diagnostics ===\n");
    std::printf("  Fuel model     : %d — %s\n", fuel->id, fuel->name.c_str());
    std::printf("  I_R            : %.2f BTU/ft²/min\n", ros.I_R);
    std::printf("  φ_w (wind)     : %.4f\n", ros.phi_w);
    std::printf("  φ_s (slope)    : %.4f\n", ros.phi_s);
    std::printf("  ξ (flux ratio) : %.4f\n", ros.xi);
    std::printf("  ROS head       : %.3f ft/min\n", ros.ROS_head);
    std::printf("  ROS back       : %.3f ft/min\n", ros.ROS_back);
    std::printf("  ROS flank      : %.3f ft/min\n", ros.ROS_flank);
    std::printf("  Head direction : %.1f° (0=N, CW)\n", ros.head_dir_deg);
    std::printf("  Grid cell size : %.0f ft\n", CELL_SIZE_FT);
    std::printf("  Grid size      : %d × %d = %.2f × %.2f ac\n",
                GRID_ROWS, GRID_COLS,
                (GRID_ROWS * CELL_SIZE_FT * CELL_SIZE_FT) / 43560.f,
                (GRID_COLS * CELL_SIZE_FT * CELL_SIZE_FT) / 43560.f);
    std::printf("=================================\n\n");

    // -----------------------------------------------------------------------
    // 5. Initialize grid and ignite center cell
    // -----------------------------------------------------------------------
    FireGrid grid(GRID_ROWS, GRID_COLS, CELL_SIZE_FT, env);
    int ign_r = GRID_ROWS / 2;
    int ign_c = GRID_COLS / 2;
    grid.ignite(ign_r, ign_c, 0.f);

    // -----------------------------------------------------------------------
    // 6. Run propagation
    // -----------------------------------------------------------------------
    grid.run(MAX_TIME_MIN);

    // -----------------------------------------------------------------------
    // 7. ASCII snapshots
    // -----------------------------------------------------------------------
    std::printf("Legend: '.' unburned   '*' active front   '#' burned\n\n");
    for (float t = SNAPSHOT_INTERVAL; t <= MAX_TIME_MIN + 0.5f; t += SNAPSHOT_INTERVAL) {
        std::printf("--- t = %.0f min ---\n", t);
        grid.printASCII(t);
        std::printf("\n");
    }

    // -----------------------------------------------------------------------
    // 8. Write burn_times.csv
    //    Columns: row, col, burn_time_min
    //    -1 = never burned within simulation window
    // -----------------------------------------------------------------------
    const char* csv_path = "burn_times.csv";
    FILE* fp = std::fopen(csv_path, "w");
    if (!fp) {
        std::fprintf(stderr, "Cannot open %s for writing\n", csv_path);
        return 1;
    }
    std::fprintf(fp, "row,col,burn_time_min\n");
    for (int r = 0; r < GRID_ROWS; ++r)
        for (int c = 0; c < GRID_COLS; ++c)
            std::fprintf(fp, "%d,%d,%.4f\n", r, c, grid.at(r, c).burn_time);
    std::fclose(fp);

    std::printf("Burn times written to %s\n", csv_path);

    // -----------------------------------------------------------------------
    // 9. Summary stats
    // -----------------------------------------------------------------------
    int burned = 0;
    for (int r = 0; r < GRID_ROWS; ++r)
        for (int c = 0; c < GRID_COLS; ++c)
            if (grid.at(r, c).burn_time >= 0.f) ++burned;

    float burned_ac = (float)burned * CELL_SIZE_FT * CELL_SIZE_FT / 43560.f;
    std::printf("Cells burned: %d / %d  (%.2f acres)\n", burned, GRID_ROWS * GRID_COLS, burned_ac);

    return 0;
}
