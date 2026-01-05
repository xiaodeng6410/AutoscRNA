#!/bin/bash
set -e
set -u
set -o pipefail

# ===============================================
# CellBender Adaptive Learning Rate Pipeline
# Auto-adjusts learning rate based on convergence
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
MAGENTA='\033[0;35m'
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

log_lr() {
    echo -e "${MAGENTA}[LEARNING RATE]${NC} $1"
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

# ===============================================
# Learning Rate Convergence Check
# ===============================================
check_learning_curve_normal() {
    local report_file=$1
    
    if [ ! -f "$report_file" ]; then
        log_warning "Report file not found: $report_file"
        return 1
    fi
    
    # Check if the report contains "This learning curve looks normal"
    if grep -q "This learning curve looks normal" "$report_file"; then
        log_success "Learning curve is normal - convergence achieved!"
        return 0
    else
        log_warning "Learning curve not normal - needs adjustment"
        return 1
    fi
}

# ===============================================
# Run CellBender with Adaptive Learning Rate
# ===============================================
run_cellbender_adaptive() {
    local input_h5=$1
    local output_file=$2
    local expected_cells=$3
    local total_droplets=$4
    local log_file=$5
    local sample_id=$6
    
    local initial_lr=0.0001
    local current_lr=$initial_lr
    local max_attempts=20
    local lr_reduction_factor=0.5
    
    # Report file is in the same directory as output file
    # Pattern: output_cellbender.h5 -> output_cellbender_report.html
    local output_dir=$(dirname "$output_file")
    local output_basename=$(basename "$output_file" .h5)
    local report_file="${output_dir}/${output_basename}_report.html"
    
    for attempt in $(seq 1 $max_attempts); do
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log_lr "Attempt $attempt/$max_attempts - Learning rate: $current_lr"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        # Show GPU memory before running
        MEM_BEFORE=$(get_gpu_memory_usage $GPU_ID)
        log_gpu "GPU $GPU_ID memory before: ${MEM_BEFORE}MB"
        
        # Run CellBender
        log_info "Running CellBender..."
        
        set +e
        "$cellbenderpath" remove-background \
            --cuda \
            --input "$input_h5" \
            --output "$output_file" \
            --expected-cells "$expected_cells" \
            --total-droplets-included "$total_droplets" \
            --learning-rate "$current_lr" \
            --num-training-tries 20 \
            --learning-rate-retry-mult 0.5 \
            --fpr 0.01 \
            --epochs 150 \
            > "$log_file" 2>&1
        
        CB_EXIT=$?
        set -e
        
        # Show GPU memory after running
        MEM_AFTER=$(get_gpu_memory_usage $GPU_ID)
        log_gpu "GPU $GPU_ID memory after: ${MEM_AFTER}MB"
        
        # Check if CellBender execution failed
        if [ $CB_EXIT -ne 0 ]; then
            log_error "CellBender execution failed (exit code: $CB_EXIT)"
            
            if [ $attempt -lt $max_attempts ]; then
                current_lr=$(echo "scale=10; $current_lr * $lr_reduction_factor" | bc)
                log_warning "Reducing learning rate to $current_lr and retrying..."
                continue
            else
                log_error "Maximum attempts reached. Giving up on sample: $sample_id"
                return 1
            fi
        fi
        
        # Check if output file was created
        if [ ! -f "$output_file" ]; then
            log_error "Output file not created"
            
            if [ $attempt -lt $max_attempts ]; then
                current_lr=$(echo "scale=10; $current_lr * $lr_reduction_factor" | bc)
                log_warning "Reducing learning rate to $current_lr and retrying..."
                continue
            else
                log_error "Maximum attempts reached. Giving up on sample: $sample_id"
                return 1
            fi
        fi
        
        # Check learning curve in report
        log_info "Checking learning curve convergence..."
        
        if check_learning_curve_normal "$report_file"; then
            log_success "Sample converged successfully with learning rate: $current_lr"
            log_success "Total attempts needed: $attempt"
            return 0
        else
            if [ $attempt -lt $max_attempts ]; then
                log_warning "Learning curve not normal, adjusting parameters..."
                
                # Backup previous attempt
                local backup_dir="$(dirname "$output_file")/attempts"
                mkdir -p "$backup_dir"
                
                if [ -f "$output_file" ]; then
                    mv "$output_file" "$backup_dir/output_cellbender_attempt${attempt}.h5"
                    log_info "Backed up attempt $attempt output"
                fi
                
                if [ -f "$report_file" ]; then
                    mv "$report_file" "$backup_dir/output_cellbender_report_attempt${attempt}.html"
                    log_info "Backed up attempt $attempt report"
                fi
                
                # Reduce learning rate
                current_lr=$(echo "scale=10; $current_lr * $lr_reduction_factor" | bc)
                log_lr "Reducing learning rate to $current_lr for next attempt"
            else
                log_error "Maximum attempts ($max_attempts) reached without achieving normal learning curve"
                log_warning "Using last generated output, but results may not be optimal"
                return 2  # Return 2 to indicate partial success
            fi
        fi
    done
    
    return 1
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

# Check if bc is available (needed for floating point arithmetic)
if ! command -v bc &> /dev/null; then
    log_error "bc command not found. Please install bc for learning rate calculations."
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
PARTIAL_SUCCESS_SAMPLES=()
START_TIME=$(date +%s)

while IFS= read -r id; do
    # Skip empty lines and comments
    [[ -z "$id" || "$id" =~ ^[[:space:]]*# ]] && continue
    
    # Trim whitespace
    id=$(echo "$id" | xargs)
    
    echo ""
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
    # Run CellBender with Adaptive Learning Rate
    # ===============================================
    SAMPLE_OUTPUT="$OUTPUT_DIR/${id}"
    mkdir -p "$SAMPLE_OUTPUT"
    
    OUTPUT_FILE="$SAMPLE_OUTPUT/output_cellbender.h5"
    LOG_FILE="$LOG_DIR/${id}_cellbender.log"
    
    log_info "Running CellBender with adaptive learning rate..."
    log_info "Log file: $LOG_FILE"
    
    set +e
    run_cellbender_adaptive "$INPUT_H5" "$OUTPUT_FILE" "$EXPECTED_CELLS" "$TOTAL_DROPLETS" "$LOG_FILE" "$id"
    ADAPTIVE_EXIT=$?
    set -e
    
    SAMPLE_END=$(date +%s)
    SAMPLE_TIME=$((SAMPLE_END - SAMPLE_START))
    
    if [ $ADAPTIVE_EXIT -eq 0 ]; then
        log_success "Sample $id completed successfully in ${SAMPLE_TIME}s"
        log_info "Output: $OUTPUT_FILE"
        
        # Check if output files exist
        if [ -f "$OUTPUT_FILE" ]; then
            FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
            log_success "Output H5 file created (size: $FILE_SIZE)"
        fi
        
        REPORT_FILE="${OUTPUT_FILE%.h5}_report.html"
        if [ -f "$REPORT_FILE" ]; then
            log_success "Report HTML created: $REPORT_FILE"
        fi
        
        SUCCESSFUL_SAMPLES+=("$id")
    elif [ $ADAPTIVE_EXIT -eq 2 ]; then
        log_warning "Sample $id completed with suboptimal convergence in ${SAMPLE_TIME}s"
        PARTIAL_SUCCESS_SAMPLES+=("$id (suboptimal convergence)")
    else
        log_error "Sample $id failed after all attempts in ${SAMPLE_TIME}s"
        log_error "Check log file for details: $LOG_FILE"
        FAILED_SAMPLES+=("$id (all attempts failed)")
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

echo ""
echo "========================================="
echo "PROCESSING SUMMARY"
echo "========================================="

TOTAL_SAMPLES=$((${#SUCCESSFUL_SAMPLES[@]} + ${#PARTIAL_SUCCESS_SAMPLES[@]} + ${#FAILED_SAMPLES[@]}))
echo "Total samples: $TOTAL_SAMPLES"
echo "Fully successful: ${#SUCCESSFUL_SAMPLES[@]}"
echo "Partial success: ${#PARTIAL_SUCCESS_SAMPLES[@]}"
echo "Failed: ${#FAILED_SAMPLES[@]}"
echo "Total time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""

if [ ${#SUCCESSFUL_SAMPLES[@]} -gt 0 ]; then
    log_success "Fully successful samples:"
    for sample in "${SUCCESSFUL_SAMPLES[@]}"; do
        echo "  ✓ $sample"
    done
    echo ""
fi

if [ ${#PARTIAL_SUCCESS_SAMPLES[@]} -gt 0 ]; then
    log_warning "Partial success samples (check results carefully):"
    for sample in "${PARTIAL_SUCCESS_SAMPLES[@]}"; do
        echo "  ⚠ $sample"
    done
    echo ""
fi

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    log_error "Failed samples:"
    for sample in "${FAILED_SAMPLES[@]}"; do
        echo "  ✗ $sample"
    done
    echo ""
    
    # Show final GPU status
    show_gpu_status
    
    exit 1
fi

if [ ${#PARTIAL_SUCCESS_SAMPLES[@]} -gt 0 ]; then
    log_warning "All samples processed, but some with suboptimal convergence"
    log_warning "Please review the samples marked with ⚠"
else
    log_success "All samples processed successfully with optimal convergence!"
fi

# Show final GPU status
show_gpu_status

echo "========================================="