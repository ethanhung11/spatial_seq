#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"/../..

# ================== PROCESS ARGS ==================

# arg defaults
transcriptome="./references/refdata-gex-GRCm39-2024-A"
probe="./references/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv"
cores=10
memusage=500

show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i, --inputcsv CSV      CSV with inputs to process (required)"
    echo "  -o, --output FILE       Output file (optional)"
    echo "  --transcriptome FILE    Transcriptome (optional), default is ./references/refdata-gex-GRCm39-2024-A"
    echo "  --probe FILE            Probe (optional), default is ./references/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv"
    echo "  --cores CORES           Cores (optional), default is 10"
    echo "  --mem MBS               Memory Usage (optional), default is 500"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --inputcsv experiment.csv --output output.out"
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--inputcsv)
            input_csv="$2"
            shift 2
            ;;
        -o|--output)
            outs="$2"
            shift 2
            ;;
        --transcriptome)
            transcriptome="$2"
            shift 2
            ;;
        --probe)
            probe="$2"
            shift 2
            ;;
        --cores)
            cores="$2"
            shift 2
            ;;
        --mem)
            memusage="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$input_csv" ]; then
    echo "Error: --inputcsv is required"
    show_usage
    exit 1
fi

if [ ! -f "$input_csv" ]; then
    echo "Error: Input CSV file '$input_csv' not found"
    exit 1
fi

if [ -n "$outs" ]; then
    output_file=./outs/"$outs"
    mkdir -p "$(dirname "$output_file")"
    exec > "$output_file"

    echo "========== Script Output =========="
    echo "Date: $(date)"
    echo "Input CSV: $input_csv"
    echo "==================================="
    echo
fi

# ================== BEGIN SCRIPT ==================

echo $(pwd)

while IFS=',' read -r sample directory slide area cytassist image; do
    # Skip header
    [ "$sample" = "sample" ] && continue
    
    # Print variables
    echo "Processing $sample"
    echo "* FastQ files:" $(find "$directory/$sample" -name "*.fastq.gz" -type f)
    echo "* Slide: $slide"
    echo "* Area: $area"
    echo "* Cytassist: $cytassist"
    echo "* Full Image: $image"

    # Create output directory
    resultdir="./data/spaceranger/$(basename $directory)/$sample"
    echo "$resultdir"
    mkdir -p "$resultdir"

    # Run spaceranger
    # Cmdline usage: https://www.10xgenomics.com/support/software/space-ranger/latest/tutorials/count-ffpe-tutorial
    # Inputs info: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/input-overview
    time spaceranger count --id "$sample" \
        --fastqs "$directory/$sample" \
        --image "$directory/$sample/$image" \
        --cytaimage "$directory/$sample/$cytassist" \
        --slide "$slide" \
        --area "$area" \
        --transcriptome $transcriptome \
        --probe-set $probe \
        --reorient-images true \
        --create-bam false \
        --output-dir "$resultdir" \
        --localcores $cores \
        --localmem $memusage

    echo
done < "$input_csv"

echo "COMPLETE!"