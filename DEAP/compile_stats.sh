#!/bin/bash
set -euo pipefail

# compile_stats.sh - Generate build information JSON file
# This script collects build metadata and creates a JSON file with compilation details
#
# Usage: compile_stats.sh [target_directory]
#   target_directory: Directory where the compile_stats.json should be placed
#                    If not provided, defaults to current directory

# Determine output directory
if [ $# -ge 1 ]; then
    TARGET_DIR="$1"
    # Create target directory if it doesn't exist
    mkdir -p "$TARGET_DIR"
    OUTPUT_FILE="$TARGET_DIR/compile_stats.json"
else
    OUTPUT_FILE="./compile_stats.json"
fi



# Function to get git information
get_git_info() {
    local git_branch=""
    local git_hash=""
    local git_status=""
    local build_type=""
    
    # Check if we're in a git repository
    if git rev-parse --git-dir > /dev/null 2>&1; then
        # Get current branch name
        git_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
        
        # Get current commit hash
        git_hash=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
        
        # Check if there are uncommitted changes
        if git diff-index --quiet HEAD -- 2>/dev/null; then
            git_status="clean"
            build_type="prod"
        else
            git_status="dirty"
            build_type="dev"
        fi
    else
        git_branch="unknown"
        git_hash="unknown"
        git_status="unknown"
        build_type="unknown"
    fi
    
    echo "$git_branch" "$git_hash" "$git_status" "$build_type"
}

# Function to create JSON
create_compile_stats_json() {
    local build_host=$(hostname)
    local build_user=$(whoami)
    local build_timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    
    # Get git information
    read -r git_branch git_hash git_status build_type <<< "$(get_git_info)"
    
    # Create directory if it doesn't exist
    mkdir -p "$(dirname "$OUTPUT_FILE")"
    
    # Write JSON structure directly to output file
    cat > "$OUTPUT_FILE" << EOF
{
    "build_host": "$build_host",
    "build_user": "$build_user", 
    "git_branch": "$git_branch",
    "git_hash": "$git_hash",
    "git_status": "$git_status",
    "build_type": "$build_type",
    "build_timestamp": "$build_timestamp"
}
EOF
    
    #echo "Build statistics written to: $OUTPUT_FILE"
}

# Main execution
main() {
    echo "Generating compile statistics..."
    create_compile_stats_json
}

# Run main function
main "$@"
