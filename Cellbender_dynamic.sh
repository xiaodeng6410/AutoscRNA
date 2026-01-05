#!/bin/bash
set -e
set -u
set -o pipefail

# ===============================================
# CellBender Simplified Pipeline
# Uses CellBender's built-in learning rate tuning
# ===============================================

CALL_DIR=$1
GPU_ID="${2:-}"  # Optional GPU ID parameter
GETCELL="/data_result/dengys/scRNA/AutoScRNA/Script/Getcellnumber.py"

# Fixed CellBender path (always use cellbender2)
cellbenderpath="/home/dengys/anaconda3/envs/cellbender2/bin/cellbender"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_gpu() {
    echo -e "${CYAN}[GPU]${NC} $1"
}

# ===============================================
# GPU Selection and Status Check
# ===============================================
check_gpu_available() {
    if ! command -v nvidia-smi &> /dev/null; then
        log_error "nvidia-smi not found. Cannot detect GPU status."
        return 1
    fi
    return 0
}

get_gpu_memory_usage() {
    local gpu_id=$1
    # Get memory usage in MB
    nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$gpu_id" 2>/dev/null || echo "N/A"
}

get_gpu_memory_total() {
    local gpu_id=$1
    # Get total memory in MB
    nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits -i "$gpu_id" 2>/dev/null || echo "N/A"
}

get_gpu_utilization() {
    local gpu_id=$1
    # Get GPU utilization percentage
    nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -i "$gpu_id" 2>/dev/null || echo "N/A"
}

show_gpu_status() {
    echo ""
    log_gpu "Current GPU Status:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    for i in 0 1; do
        local mem_used=$(get_gpu_memory_usage $i)
        local mem_total=$(get_gpu_memory_total $i)
        local gpu_util=$(get_gpu_utilization $i)
        local gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader -i $i 2>/dev/null || echo "Unknown")
        
        if [ "$mem_used" != "N/A" ] && [ "$mem_total" != "N/A" ]; then
            local mem_percent=$((mem_used * 100 / mem_total))
            local status_color=$GREEN
            
            if [ $mem_percent -gt 80 ]; then
                status_color=$RED
            elif [ $mem_percent -gt 50 ]; then
                status_color=$YELLOW
            fi
            
            echo -e "GPU $i: ${CYAN}${gpu_name}${NC}"
            echo -e "  Memory: ${status_color}${mem_used}MB / ${mem_total}MB (${mem_percent}%)${NC}"
            echo -e "  Utilization: ${gpu_util}%"
        else
            echo -e "GPU $i: ${RED}Not available${NC}"
        fi
        echo ""
    done
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}

select_best_gpu() {
    log_info "Auto-selecting GPU with most free memory..."
    
    local gpu0_used=$(get_gpu_memory_usage 0)
    local gpu1_used=$(get_gpu_memory_usage 1)
    
    if [ "$gpu0_used" = "N/A" ] && [ "$gpu1_used" = "N/A" ]; then
        log_error "Cannot detect GPU memory usage"
        return 1
    fi
    
    if [ "$gpu0_used" = "N/A" ]; then
        echo "1"
        return 0
    fi
    
    if [ "$gpu1_used" = "N/A" ]; then
        echo "0"
        return 0
    fi
    
    # Select GPU with less memory usage
    if [ "$gpu0_used" -lt "$gpu1_used" ]; then
        echo "0"
    else
        echo "1"
    fi
}

# Check GPU availability
if ! check_gpu_available; then
    log_error "GPU not available. Exiting."
    exit 1
fi

# Show current GPU status
show_gpu_status

# Determine which GPU to use
if [ -z "$GPU_ID" ]; then
    GPU_ID=$(select_best_gpu)
    log_success "Auto-selected GPU: $GPU_ID"
else
    # Validate user-specified GPU ID
    if ! [[ "$GPU_ID" =~ ^[0-1]$ ]]; then
        log_error "GPU_ID must be 0 or 1 (got: '$GPU_ID')"
        exit 1
    fi
    log_info "Using user-specified GPU: $GPU_ID"
fi

# Set CUDA device
export CUDA_VISIBLE_DEVICES=$GPU_ID

# Get GPU name for logging
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader -i $GPU_ID 2>/dev/null || echo "Unknown GPU")

# ===============================================
# Verify CellBender Installation
# ===============================================
if [ ! -f "$cellbenderpath" ]; then
    log_error "CellBender executable not found at $cellbenderpath"
    exit 1
fi

# ===============================================
# Path Configuration
# ===============================================
filename="${CALL_DIR}/filename.txt"
INPUT_DIR="${CALL_DIR}/cellranger_out"
OUTPUT_DIR="${CALL_DIR}/cellbender_out"
LOG_DIR="${OUTPUT_DIR}/logs"

# ===============================================
# Input Validation
# ===============================================
if [ ! -d "$CALL_DIR" ]; then
    log_error "Directory $CALL_DIR does not exist"
    exit 1
fi

if [ ! -f "$filename" ]; then
    log_error "filename.txt not found at $filename"
    exit 1
fi

if [ ! -f "$GETCELL" ]; then
    log_error "Getcellnumber.py not found at $GETCELL"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

echo ""
log_success "Validation passed"
log_info "Input directory: $INPUT_DIR"
log_info "Output directory: $OUTPUT_DIR"
log_info "Using CellBender: $cellbenderpath"
log_gpu "Using GPU $GPU_ID: $GPU_NAME"
echo ""

# ===============================================
# Process Each Sample
# ===============================================
FAILED_SAMPLES=()
SUCCESSFUL_SAMPLES=()
START_TIME=$(date +%s)

while IFS= read -r id; do
    # Skip empty lines and comments
    [[ -z "$id" || "$id" =~ ^[[:space:]]*# ]] && continue
    
    # Trim whitespace
    id=$(echo "$id" | xargs)
    
    echo "========================================="
    log_info "Processing sample: $id"
    echo "========================================="
    
    SAMPLE_START=$(date +%s)
    
    # Check input file
    INPUT_H5="$INPUT_DIR/${id}/outs/raw_feature_bc_matrix.h5"
    if [ ! -f "$INPUT_H5" ]; then
        log_warning "Input file not found: $INPUT_H5"
        log_warning "Skipping sample $id"
        FAILED_SAMPLES+=("$id (input not found)")
        echo ""
        continue
    fi
    
    # ===============================================
    # Extract Cell Count from CellRanger Output
    # ===============================================
    SAMPLE_DIR="$INPUT_DIR/${id}"
    log_info "Extracting cell count from CellRanger output..."
    
    set +e
    N_CELLS_RAW=$(python "$GETCELL" "$SAMPLE_DIR" 2>/dev/null)
    GETCELL_EXIT=$?
    set -e
    
    # Validate N_CELLS is a positive integer
    if [ $GETCELL_EXIT -eq 0 ] && [[ "$N_CELLS_RAW" =~ ^[0-9]+$ ]] && [ "$N_CELLS_RAW" -gt 0 ]; then
        N_CELLS=$N_CELLS_RAW
        log_success "Estimated cells from CellRanger: $N_CELLS"
    else
        log_warning "Could not extract valid cell count (got: '$N_CELLS_RAW'), using default value"
        N_CELLS=10000
    fi
    
    # ===============================================
    # Calculate CellBender Parameters
    # ===============================================
    EXPECTED_CELLS=$((N_CELLS + 2000))
    TOTAL_DROPLETS=$((N_CELLS * 3))
    
    log_info "CellBender parameters:"
    echo "  --expected-cells: $EXPECTED_CELLS"
    echo "  --total-droplets-included: $TOTAL_DROPLETS"
    echo ""
    
    # ===============================================
    # Run CellBender
    # ===============================================
    SAMPLE_OUTPUT="$OUTPUT_DIR/${id}"
    mkdir -p "$SAMPLE_OUTPUT"
    
    OUTPUT_FILE="$SAMPLE_OUTPUT/output_cellbender.h5"
    LOG_FILE="$LOG_DIR/${id}_cellbender.log"
    
    log_info "Running CellBender with automatic learning rate tuning..."
    log_info "Log file: $LOG_FILE"
    
    # Show GPU memory before running
    MEM_BEFORE=$(get_gpu_memory_usage $GPU_ID)
    log_gpu "GPU $GPU_ID memory before: ${MEM_BEFORE}MB"
    
    set +e
    "$cellbenderpath" remove-background \
        --cuda \
        --input "$INPUT_H5" \
        --output "$OUTPUT_FILE" \
        --expected-cells "$EXPECTED_CELLS" \
        --total-droplets-included "$TOTAL_DROPLETS" \
        --learning-rate 0.0001 \
        --num-training-tries 20 \
        --final-elbo-fail-fraction 0.3 \
        --learning-rate-retry-mult 0.5 \
        --fpr 0.01 \
        --epochs 150 \
        > "$LOG_FILE" 2>&1
    
    CB_EXIT=$?
    set -e
    
    # Show GPU memory after running
    MEM_AFTER=$(get_gpu_memory_usage $GPU_ID)
    log_gpu "GPU $GPU_ID memory after: ${MEM_AFTER}MB"
    
    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))
    
    if [ $CB_EXIT -eq 0 ]; then
        log_success "Sample $id completed successfully in ${SAMPLE_TIME}s"
        log_info "Output: $OUTPUT_FILE"
        
        # Check if output files exist
        if [ -f "$OUTPUT_FILE" ]; then
            FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
            log_success "Output H5 file created (size: $FILE_SIZE)"
        else
            log_warning "Output H5 file not found"
        fi
        
        REPORT_FILE="${OUTPUT_FILE%.h5}_report.html"
        if [ -f "$REPORT_FILE" ]; then
            log_success "Report HTML created: $REPORT_FILE"
        else
            log_warning "Report HTML not found"
        fi
        
        SUCCESSFUL_SAMPLES+=("$id")
    else
        log_error "CellBender execution failed for sample: $id (exit code: $CB_EXIT)"
        log_error "Check log file for details: $LOG_FILE"
        FAILED_SAMPLES+=("$id (exit code: $CB_EXIT)")
    fi
    
    echo ""
done < "$filename"

# ===============================================
# Final Summary
# ===============================================
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
HOURS=$((TOTAL_TIME / 3600))
MINUTES=$(((TOTAL_TIME % 3600) / 60))
SECONDS=$((TOTAL_TIME % 60))

echo "========================================="
echo "PROCESSING SUMMARY"
echo "========================================="

TOTAL_SAMPLES=$((${#SUCCESSFUL_SAMPLES[@]} + ${#FAILED_SAMPLES[@]}))
echo "Total samples: $TOTAL_SAMPLES"
echo "Successful: ${#SUCCESSFUL_SAMPLES[@]}"
echo "Failed: ${#FAILED_SAMPLES[@]}"
echo "Total time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""

if [ ${#SUCCESSFUL_SAMPLES[@]} -gt 0 ]; then
    log_success "Successful samples:"
    for sample in "${SUCCESSFUL_SAMPLES[@]}"; do
        echo "  ✓ $sample"
    done
    echo ""
fi

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    log_warning "Failed samples:"
    for sample in "${FAILED_SAMPLES[@]}"; do
        echo "  ✗ $sample"
    done
    echo ""
    
    # Show final GPU status
    show_gpu_status
    
    exit 1
fi

log_success "All samples processed successfully!"

# Show final GPU status
show_gpu_status

echo "========================================="