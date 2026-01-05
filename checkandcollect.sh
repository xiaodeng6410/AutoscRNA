#!/bin/bash

# ===============================================
# FastQ Processing Script v3.0
# Purpose: MD5 check → mapping → rename → organize
# Author: Generated for scRNA data processing
# ===============================================

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ===============================================
# Utility Functions
# ===============================================

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

# ===============================================
# Step 1: MD5 Validation
# ===============================================

check_md5() {
    local input_dir="$1"
    local log_file="${input_dir}/md5check.log"
    
    log_info "Starting MD5 validation..."
    
    # Find all .md5 files
    local md5_files=$(find "$input_dir" -maxdepth 1 -name "*.md5" -type f)
    
    if [ -z "$md5_files" ]; then
        log_warning "No .md5 files found in $input_dir"
        echo "NO MD5 FILES FOUND" > "$log_file"
        return 1
    fi
    
    local all_ok=true
    
    # Process each .md5 file
    while IFS= read -r md5_file; do
        log_info "Checking: $(basename "$md5_file")"
        
        # Change to the directory containing the md5 file for relative path checking
        pushd "$(dirname "$md5_file")" > /dev/null
        
        # Run md5sum check and append to log
        if md5sum -c "$(basename "$md5_file")" >> "$log_file" 2>&1; then
            log_success "$(basename "$md5_file") - All files OK"
        else
            log_error "$(basename "$md5_file") - Some files FAILED"
            all_ok=false
        fi
        
        popd > /dev/null
    done <<< "$md5_files"
    
    # Write final status
    if [ "$all_ok" = true ]; then
        echo "" >> "$log_file"
        echo "ALL FILES OK" >> "$log_file"
        log_success "MD5 validation completed - ALL FILES OK"
    else
        echo "" >> "$log_file"
        echo "ERROR: md5 mismatch detected" >> "$log_file"
        log_error "MD5 validation completed - ERRORS DETECTED (continuing anyway)"
    fi
}

# ===============================================
# Step 2: Parse Sample Information from .md5
# ===============================================

parse_md5_mapping() {
    local input_dir="$1"
    local projection_file="${input_dir}/projection.txt"
    
    log_info "Parsing sample information from .md5 files..."
    
    # Clear projection file
    > "$projection_file"
    
    # Find all .md5 files
    local md5_files=$(find "$input_dir" -maxdepth 1 -name "*.md5" -type f)
    
    declare -A folder_oes_map
    
    while IFS= read -r md5_file; do
        while IFS= read -r line; do
            # Skip empty lines
            [ -z "$line" ] && continue
            
            # Parse line: md5hash  path/to/file
            local filepath=$(echo "$line" | awk '{print $2}')
            
            # Extract folder name and OES number
            if [[ "$filepath" =~ ^([^/]+)/OES([0-9]+) ]]; then
                local folder_name="${BASH_REMATCH[1]}"
                local oes_number="OES${BASH_REMATCH[2]}"
                
                folder_oes_map["$folder_name"]="$oes_number"
            fi
        done < "$md5_file"
    done <<< "$md5_files"
    
    # Write to projection.txt
    for folder in "${!folder_oes_map[@]}"; do
        echo -e "${folder}\t${folder_oes_map[$folder]}" >> "$projection_file"
    done
    
    log_success "Mapping saved to projection.txt ($(wc -l < "$projection_file") entries)"
}

# ===============================================
# Step 3: Rename FastQ Files (OE Company Rules)
# ===============================================

rename_fastq_OE() {
    local input_dir="$1"
    
    log_info "Renaming fastq files (OE format)..."
    
    # Read projection mapping
    declare -A folder_to_oes
    while IFS=$'\t' read -r folder oes; do
        folder_to_oes["$folder"]="$oes"
    done < "${input_dir}/projection.txt"
    
    # Process each folder
    for folder_name in "${!folder_to_oes[@]}"; do
        local folder_path="${input_dir}/${folder_name}"
        local oes_number="${folder_to_oes[$folder_name]}"
        
        if [ ! -d "$folder_path" ]; then
            log_warning "Folder not found: $folder_path"
            continue
        fi
        
        log_info "Processing folder: $folder_name → $oes_number"
        
        # Find all R1 files (sorted for consistent pairing)
        local r1_files=($(find "$folder_path" -name "*R1.fastq.gz" -type f | sort))
        local lane_counter=1
        
        for r1_file in "${r1_files[@]}"; do
            # Derive R2 file name
            local r2_file="${r1_file/R1.fastq.gz/R2.fastq.gz}"
            
            if [ ! -f "$r2_file" ]; then
                log_warning "R2 file not found for: $(basename "$r1_file")"
                continue
            fi
            
            # Generate new names
            local lane_num=$(printf "L%03d" $lane_counter)
            local new_r1="${folder_path}/${oes_number}_S1_${lane_num}_R1_001.fastq.gz"
            local new_r2="${folder_path}/${oes_number}_S1_${lane_num}_R2_001.fastq.gz"
            
            # Rename files
            mv "$r1_file" "$new_r1"
            mv "$r2_file" "$new_r2"
            
            log_success "  Renamed: $(basename "$r1_file") → $(basename "$new_r1")"
            log_success "  Renamed: $(basename "$r2_file") → $(basename "$new_r2")"
            
            ((lane_counter++))
        done
    done
}

# ===============================================
# Step 4: Rename Sample Folders
# ===============================================

rename_folders() {
    local input_dir="$1"
    local filename_file="${input_dir}/filename.txt"
    
    log_info "Renaming sample folders..."
    
    # Clear filename.txt
    > "$filename_file"
    
    # Read projection mapping
    declare -A folder_to_oes
    while IFS=$'\t' read -r folder oes; do
        folder_to_oes["$folder"]="$oes"
    done < "${input_dir}/projection.txt"
    
    # Rename each folder
    for folder_name in "${!folder_to_oes[@]}"; do
        local old_path="${input_dir}/${folder_name}"
        local oes_number="${folder_to_oes[$folder_name]}"
        local new_path="${input_dir}/${oes_number}"
        
        if [ -d "$old_path" ]; then
            mv "$old_path" "$new_path"
            log_success "Renamed folder: $folder_name → $oes_number"
            echo "$oes_number" >> "$filename_file"
        else
            log_warning "Folder not found: $old_path"
        fi
    done
    
    log_success "Folder names saved to filename.txt"
}

# ===============================================
# Step 5: Organize into rawdata/
# ===============================================

organize_rawdata() {
    local input_dir="$1"
    local rawdata_dir="${input_dir}/rawdata"
    
    log_info "Organizing files into rawdata/..."
    
    # Create rawdata directory
    mkdir -p "$rawdata_dir"
    
    # Read all OES folder names from filename.txt
    while IFS= read -r oes_folder; do
        local source_path="${input_dir}/${oes_folder}"
        local dest_path="${rawdata_dir}/${oes_folder}"
        
        if [ -d "$source_path" ]; then
            mv "$source_path" "$dest_path"
            log_success "Moved: $oes_folder → rawdata/"
        fi
    done < "${input_dir}/filename.txt"
    
    log_success "All samples organized in rawdata/"
}

# ===============================================
# Company-Specific Rename Functions (Extensible)
# ===============================================

rename_fastq_by_company() {
    local input_dir="$1"
    local company="$2"
    
    case "$company" in
        OE)
            rename_fastq_OE "$input_dir"
            ;;
        *)
            log_error "Unsupported company: $company"
            exit 1
            ;;
    esac
}

# ===============================================
# Main Execution
# ===============================================

main() {
    # Check arguments
    if [ $# -ne 2 ]; then
        echo "Usage: $0 <input_dir> <company>"
        echo "Example: $0 /data_result/dengys/scRNA/2511_ZTT_lipid OE"
        exit 1
    fi
    
    local input_dir="$1"
    local company="$2"
    
    # Validate input directory
    if [ ! -d "$input_dir" ]; then
        log_error "Input directory does not exist: $input_dir"
        exit 1
    fi
    
    # Validate company
    if [ "$company" != "OE" ]; then
        log_error "Currently only 'OE' company is supported"
        exit 1
    fi
    
    log_info "=========================================="
    log_info "FastQ Processing Pipeline Started"
    log_info "Input Directory: $input_dir"
    log_info "Company: $company"
    log_info "=========================================="
    echo ""
    
    # Execute pipeline
    check_md5 "$input_dir"
    echo ""
    
    parse_md5_mapping "$input_dir"
    echo ""
    
    rename_fastq_by_company "$input_dir" "$company"
    echo ""
    
    rename_folders "$input_dir"
    echo ""
    
    organize_rawdata "$input_dir"
    echo ""
    
    log_success "=========================================="
    log_success "Pipeline completed successfully!"
    log_success "=========================================="
    log_info "Generated files:"
    log_info "  - md5check.log"
    log_info "  - projection.txt"
    log_info "  - filename.txt"
    log_info "  - rawdata/ (containing all samples)"
}

# Run main function
main "$@"
