#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"/../..

# ================== PROCESS ARGS ==================

# arg defaults
transcriptome="./references/refdata-gex-GRCm39-2024-A"
cores=20
memusage=50

show_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -i, --inputcsv CSV      CSV with inputs to process (required)"
    echo "  -i, --inputdir DIR      DIR with inputs to process (required)"
    echo "  -o, --outputdir DIR     Output directory (optional)"
    echo "  --transcriptome FILE    Transcriptome (optional), default is ./references/refdata-gex-GRCm39-2024-A"
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
        -i|--inputcsv|--inputdir)
            input="$2"
            shift 2
            ;;
        -o|--outputdir)
            outs="$2"
            shift 2
            ;;
        --transcriptome)
            transcriptome="$2"
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
if [[ -z "$input" ]]; then
    echo "Error: --inputcsv is required"
    show_usage
    exit 1
fi

if [[ ! -f "$input" && ! -d "$input" ]]; then
    echo "Error: Input file/directory '$input_csv' not found"
    exit 1
fi

if [[ -n "$outs" ]]; then
    output_file=./outs/cellcounting_"$outs".out
    mkdir -p "$(dirname "$output_file")"
    exec > "$output_file"

    echo "========== Script Output =========="
    echo "Date: $(date)"
    echo "Input: $input"
    echo "==================================="
fi

# ================== BEGIN SCRIPT ==================


echo $(pwd)

if [[ -f "$input" && "$filename" == *.csv ]]; then
    while IFS=',' read -r sample directory; do
        # Skip header
        [ "$sample" = "sample" ] && continue


        # Create output directory
        resultdir="./data/cellranger/$(basename $directory)/$sample"
        echo "$resultdir"
        
        # Print variables
        echo "Processing $sample"
        echo "* FastQ files:" $(find "$directory/$sample" -name "*.fastq.gz" -type f)

        # Run cellranger
        time cellranger count --id "$sample" \
            --fastqs "$directory/$sample" \
            --transcriptome $transcriptome \
            --create-bam false \
            --output-dir "$resultdir" \
            --localcores $cores \
            --localmem $memusage

        echo
    done < "$input"

elif [[ -d "$input" ]]; then
    
    # Create output directory
    resultdir="./data/cellranger/$(basename $input)/$sample"
    echo "$resultdir"

    find "$input" -maxdepth 1 -type d ! -name "$(basename "$input")" | while read -r sample; do
        echo $sample
        time cellranger count --id "$(basename $sample)" \
            --fastqs "$sample" \
            --transcriptome $transcriptome \
            --create-bam false \
            --output-dir "$resultdir" \
            --localcores $cores \
            --localmem $memusage
    done

fi