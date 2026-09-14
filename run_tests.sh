#!/bin/bash
#
# run_tests.sh
#
# Runs the three blacklight test cases, then generates a spectrum plot for
# each using plot_spectra.py, with the MC overlay applied to all three
# cases. The MC directory for each case is taken from that input file's own
# "simulation_file" line (simulation_file points at a data file inside the
# run directory, so the script uses its containing directory).
#
# The .npz image path for each run is read directly from that input file's
# own "output_file" line, so there's no guessing about naming conventions.
# If an input file has no output_file line for some reason, this falls back
# to ./tests/<input_basename>.npz.
#
# Usage:
#   Place this script (and plot_spectra.py) in the top-level "blacklight"
#   directory (the same directory that contains ./bin/blacklight and
#   ./input/), then run:
#
#       ./run_tests.sh

#Tegan to do: modify the script so that ymin and ymax are treated similar to temperature
# also come up with a test for the scattering error bars 
BLACKLIGHT_BIN="./bin/blacklight"
INPUT_DIR="./input/"
LOG_DIR="./test_logs"
PLOT_SCRIPT="./scripts/plot_spectra.py"

# Fallback output location if an input file has no output_file line.
FALLBACK_OUTPUT_DIR="./tests"

INPUTS=(
    "example_isoth.input"
    "example_isoth_gr.input"
    "example_mc_isoth_gr.input"
)

# Blackbody temperature (K) and radius (cm) for each case, passed through to
# plot_spectra.py as --temperature/--radius. Indices line up with INPUTS
# above - i.e. TEMPERATURES[0]/RADII[0] go with INPUTS[0], and so on.
TEMPERATURES=(
    1e6
    1e6
    1e7
)
RADII=(
    1e11
    8859750.2283
    8859750.2283
)

# Redshift for each case, passed through to plot_spectra.py as
# --redshift_re. Indices line up with INPUTS above.
REDSHIFT=(
    0.0
    6.0
    6.0
)

mkdir -p "$LOG_DIR" "$FALLBACK_OUTPUT_DIR"

if [[ ${#TEMPERATURES[@]} -ne ${#INPUTS[@]} || ${#RADII[@]} -ne ${#INPUTS[@]} || ${#REDSHIFT[@]} -ne ${#INPUTS[@]} ]]; then
    echo "Error: TEMPERATURES, RADII, and REDSHIFT must each have exactly one entry per INPUTS entry."
    exit 1
fi

# ----------------------------------------------------------------------------
# Sanity checks
# ----------------------------------------------------------------------------
if [[ ! -x "$BLACKLIGHT_BIN" ]]; then
    echo "Error: $BLACKLIGHT_BIN not found or not executable."
    echo "Make sure you're running this script from the top-level blacklight directory."
    exit 1
fi

if [[ ! -f "$PLOT_SCRIPT" ]]; then
    echo "Error: $PLOT_SCRIPT not found next to this script."
    exit 1
fi

# ----------------------------------------------------------------------------
# Extract a "key = value  # comment" style setting from a blacklight input
# file. Matches the format confirmed in example_mc_isoth_gr.input, e.g.:
#   output_file   = /path/to/file.npz  # file to be (over)written with output data
#   mc_file       = /path/to/dir/isothgr.out3.00000.athdf
# ----------------------------------------------------------------------------
extract_input_value() {
    local input_file="$1"
    local key="$2"
    grep -m1 -iE "^[[:space:]]*${key}[[:space:]]*=" "$input_file" \
        | sed -E "s/^[[:space:]]*${key}[[:space:]]*=[[:space:]]*//I" \
        | sed -E 's/#.*$//' \
        | sed -E 's/[[:space:]]+$//' \
        | xargs
}

declare -A RUN_RESULTS
declare -A PLOT_RESULTS

for i in "${!INPUTS[@]}"; do
    input="${INPUTS[$i]}"
    temperature="${TEMPERATURES[$i]}"
    radius="${RADII[$i]}"
    redshift="${REDSHIFT[$i]}"
    input_path="${INPUT_DIR}${input}"
    base="${input%.input}"
    log_file="${LOG_DIR}/${base}.log"
    plot_log_file="${LOG_DIR}/${base}_plot.log"

    if [[ ! -f "$input_path" ]]; then
        echo "Warning: $input_path not found, skipping."
        RUN_RESULTS["$input"]="MISSING"
        PLOT_RESULTS["$input"]="SKIPPED"
        continue
    fi

    echo "=============================================="
    echo "Running: $BLACKLIGHT_BIN $input_path"
    echo "Log:     $log_file"
    echo "=============================================="

    "$BLACKLIGHT_BIN" "$input_path" > "$log_file" 2>&1
    status=$?

    if [[ $status -eq 0 ]]; then
        echo "  -> SUCCESS"
        RUN_RESULTS["$input"]="SUCCESS"
    else
        echo "  -> FAILED (exit code $status). See $log_file for details."
        RUN_RESULTS["$input"]="FAILED (exit $status)"
        PLOT_RESULTS["$input"]="SKIPPED (run failed)"
        echo
        continue
    fi

    # Determine the .npz image path from this input file's own output_file line.
    image_path="$(extract_input_value "$input_path" "output_file")"
    if [[ -z "$image_path" ]]; then
        image_path="${FALLBACK_OUTPUT_DIR}/${base}.npz"
        echo "  -> No output_file line found in $input_path; falling back to $image_path"
    fi

    if [[ ! -f "$image_path" ]]; then
        echo "  -> Warning: expected image $image_path not found; skipping plot."
        PLOT_RESULTS["$input"]="SKIPPED (no image found: $image_path)"
        echo
        continue
    fi

    plot_output="$(dirname "$image_path")/${base}_spectrum.png"
    plot_args=(--image "$image_path" --output "$plot_output" --label "$base" \
                --temperature "$temperature" --radius "$radius" --redshift_re "$redshift")

    # Use the MC overlay for all cases, with the directory taken from this
    # input file's simulation_file line (simulation_file points at a data
    # file inside the run directory, so we use its containing directory).
    simulation_file_value="$(extract_input_value "$input_path" "simulation_file")"
    if [[ -z "$simulation_file_value" ]]; then
        echo "  -> Warning: no simulation_file line was found in $input_path."
        echo "     Plotting without MC overlay."
    else
        mc_dir="$(dirname "$simulation_file_value")/"
        echo "  -> Using simulation_file directory: $mc_dir"
        plot_args+=(--mc-dir "$mc_dir")
    fi

    echo "Plotting: python3 $PLOT_SCRIPT ${plot_args[*]}"
    python3 "$PLOT_SCRIPT" "${plot_args[@]}" > "$plot_log_file" 2>&1
    plot_status=$?

    if [[ $plot_status -eq 0 ]]; then
        echo "  -> Plot SUCCESS: $plot_output"
        PLOT_RESULTS["$input"]="SUCCESS ($plot_output)"
    else
        echo "  -> Plot FAILED (exit code $plot_status). See $plot_log_file for details."
        PLOT_RESULTS["$input"]="FAILED (exit $plot_status)"
    fi
    echo
done

echo "=============================================="
echo "Summary"
echo "=============================================="
printf "%-28s %-22s %s\n" "INPUT" "BLACKLIGHT RUN" "PLOT"
for input in "${INPUTS[@]}"; do
    printf "%-28s %-22s %s\n" "$input" "${RUN_RESULTS[$input]:-NOT RUN}" "${PLOT_RESULTS[$input]:-NOT RUN}"
done