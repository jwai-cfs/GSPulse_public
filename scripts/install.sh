#!/usr/bin/env bash

# --- Error Handling ---
# -e: Exit immediately if a command exits with a non-zero status.
# -u: Treat unset variables as an error when substituting.
# -o pipefail: The return value of a pipeline is the status of the last command
#              to exit with a non-zero status, or zero if all commands exit successfully.
set -euo pipefail

# --- Python interpreter selection ---
DEFAULT_PYTHON_VERSION="3.12"
REQUESTED_PYTHON=""
USE_LATEST=false
SKIP_OCTAVE_BUILD=false
SKIP_MATLAB_BUILD=false

# --- Get Git Repository Root Directory ---
GIT_ROOT_DIR=$(git rev-parse --show-toplevel)

# --- Get the operating system name ---
OS=$(uname -s)

# --- Detect Octave path ---
if command -v octave &> /dev/null; then
    OCTAVE_BIN=$(dirname "$(command -v octave)")
    echo "Detected Octave at: $OCTAVE_BIN"
else
    echo "Warning: Octave not found in PATH. MEQ compilation may fail."
    OCTAVE_BIN=""
fi

# --- Detect Matlab path ---
if [ -z "${MATPATH:-}" ]; then
    echo "Warning: MATPATH environment variable is not set. Matlab MEQ build may fail."
    echo "To set MATPATH, use a command similar to:"
    echo '  export MATPATH="/Applications/MATLAB_R2024a.app"'
else
    echo "Detected MATPATH: $MATPATH"
fi

# --- Default variables ---
CLEAN=false
UPDATE_LOCKS=false


# --- Function to display help message ---
show_help() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo ""
    echo "This script runs a command based on the operating system."
    echo ""
    echo "Options:"
    echo "  --make-clean        Run 'make clean' before building MEQ."
    echo "  --update-locks      Resolve to the latest packages permitted by pyproject.toml for all projects (new lock files written)."
    echo "  --python <version>  Pin uv to a Python interpreter/version (default: 3.13, use 'latest' for highest allowed for each package independently)."
    echo "  --skip-octave-build  Skip building MEQ for Octave."
    echo "  --skip-matlab-build  Skip building MEQ for Matlab."
    echo "  -h, --help          Display this help message and exit."
    echo ""
}


apply_python_selection() {
    local context="$1"
    if [ "$USE_LATEST" = true ]; then
        if [ -d ".venv" ]; then
            echo "Removing existing virtual environment in $context to test latest interpreter..."
            rm -rf ".venv"
        fi

        # Try Python versions in descending order to find the highest compatible one.
        # Uses uv sync --dry-run to validate against the project's requires-python constraint.
        for version in 3.14 3.13 3.12 3.11; do
            if (cd "$context" && UV_PYTHON=$version uv sync --dry-run >/dev/null 2>&1); then
                export UV_PYTHON="$version"
                echo "Using Python interpreter for uv ($context): $UV_PYTHON (highest compatible)"
                return
            fi
        done
        # No compatible version found—fail loudly
        echo "Error: No compatible Python version found in [$context] for versions 3.14–3.11" >&2
        return 1
    else
        export UV_PYTHON="$PINNED_PYTHON"
        echo "Using Python interpreter for uv ($context): $UV_PYTHON (pinned)"
    fi
}

# --- Parse command-line arguments ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --make-clean)
            CLEAN=true
            shift
            ;;
        --python)
            if [[ $# -lt 2 ]]; then
                echo "Error: --python requires an argument" >&2
                show_help
                exit 1
            fi
            REQUESTED_PYTHON="$2"
            # Validate the requested Python version
            if [[ "$REQUESTED_PYTHON" != "latest" ]] && ! [[ "$REQUESTED_PYTHON" =~ ^3\.[0-9]{1,2}$ ]]; then
                echo "Error: Invalid Python version '$REQUESTED_PYTHON'. Must be 'latest' or match pattern '3.xx'" >&2
                show_help
                exit 1
            fi
            shift 2
            ;;
        --skip-octave-build)
            SKIP_OCTAVE_BUILD=true
            shift
            ;;
        --skip-matlab-build)
            SKIP_MATLAB_BUILD=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'" >&2
            show_help
            exit 1
            ;;
    esac
done

# --- Determine Python strategy ---
PINNED_PYTHON=""
# Case-insensitive comparison for "latest"
if [[ -n "$REQUESTED_PYTHON" ]] && [[ "$REQUESTED_PYTHON" == "latest" ]]; then
    USE_LATEST=true
elif [[ -n "$REQUESTED_PYTHON" ]]; then
    PINNED_PYTHON="$REQUESTED_PYTHON"
elif [[ -n "${UV_PYTHON:-}" ]]; then
    PINNED_PYTHON="$UV_PYTHON"
else
    PINNED_PYTHON="$DEFAULT_PYTHON_VERSION"
fi

if [[ "$USE_LATEST" = true ]]; then
    echo "Python selection: using highest compatible Python version for each project."
else
    export UV_PYTHON="$PINNED_PYTHON"
    echo "Python selection: pinned to $UV_PYTHON."
fi

# --- Change to Git Root Directory ---
# All subsequent relative paths will be based from here.
echo "Navigating to Git repository root: $GIT_ROOT_DIR"
cd "$GIT_ROOT_DIR"

# --- Update Git Submodules ---
echo "Updating all Git submodules (this may take a while)..."
git submodule update --init --recursive submodules/scs-matlab
git submodule update --init --recursive submodules/meq
echo "All Git submodules updated."

# --- sync uv ---
cd "$GIT_ROOT_DIR"
(
    cd src/python
    echo "Syncing uv in $(pwd)"
    uv sync
)

build_meq() {
    local use_octave="$1"  # "true" for Octave, "false" for Matlab
    local DEST_DIR="$2"    # Destination directory for the build

    # Build in a specified directory
    mkdir -p "${DEST_DIR}"
    # copy directory, use tar because cp -r fails
    tar -C "$GIT_ROOT_DIR/submodules/meq" --exclude='.git' -cf - . | tar -C "${DEST_DIR}" -xf - 
    cd "${DEST_DIR}"
   
    # conditional 'make clean'
    if [ "$CLEAN" = true ]; then
        echo "Running 'make clean'..."
        make clean
    fi

    # Compile MEQ for Octave
    if [ "$use_octave" = "true" ]; then
        if [ "$OS" == "Linux" ]; then
            echo "Compiling MEQ for Octave on Linux..."
            make -j USE_OCTAVE=1 USE_OPENMP=no BLASLIB=CBLAS MATPATH="$OCTAVE_BIN" tbx
            echo "... Finished compiling MEQ for Octave on Linux in ${DEST_DIR}"
        elif [ "$OS" == "Darwin" ]; then
            echo "Compiling MEQ for Octave on OSX..."
            make -j tbx USE_OCTAVE=1 USE_OPENMP=no BLASLIB=CBLAS CC=gcc \
            CBLASPATH=/opt/homebrew/opt/openblas/ LAPACKPATH=/opt/homebrew/opt/lapack/ MATPATH="$OCTAVE_BIN"
            echo "... Finished compiling MEQ for Octave on OSX in ${DEST_DIR}"
        else
            echo "Unsupported operating system: $OS, cannot compile MEQ for Octave."
        fi
    else
    # Compile MEQ for Matlab   
        if [ "$OS" == "Linux" ]; then
            make -j tbx MATPATH="$MATPATH" CC=gcc
            echo "... Finished compiling MEQ for Matlab on Linux in ${DEST_DIR}"
        elif [ "$OS" == "Darwin" ]; then
            echo "Compiling MEQ for Matlab on OSX..."
            make -j tbx CC=gcc USE_OPENMP=no MATPATH="$MATPATH"
            echo "... Finished compiling MEQ for Matlab on OSX in ${DEST_DIR}"
        else
            echo "Unsupported operating system: $OS, cannot compile MEQ for Matlab."
        fi
    fi
    
}

# --- compile MEQ -------
echo -e "\n--- compiling MEQ ---"
(
    # build for octave
    if [ "$SKIP_OCTAVE_BUILD" = false ]; then
        build_meq "true" "$GIT_ROOT_DIR/build/meq_octave_build" 
    else
        echo "Skipping MEQ build for Octave (--skip-octave-build set)"
    fi

    # build for matlab
    if [ "$SKIP_MATLAB_BUILD" = false ]; then
       build_meq "false" "$GIT_ROOT_DIR/build/meq_matlab_build" 
    else
        echo "Skipping MEQ build for Matlab (--skip-matlab-build set)"
    fi  
)

echo -e "\n--- All setup complete ---"
