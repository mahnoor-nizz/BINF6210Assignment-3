module load sra-toolkit
module load fastqc/0.12.1
module load StdEnv/2023
module load gcc/12.3
module load kraken2/2.1.6
module load bracken/3.0

module load python/3.11

pip install --user kraken-biom

VEGAN=("SRR8146944" "SRR8146963" "SRR8146968")
OMNI=("SRR8146935" "SRR8146936" "SRR8146938")

ALL_SAMPLES=("${VEGAN[@]}" "${OMNI[@]}")

mkdir -p raw_data
cd raw_data


#Download data and convert to FASTQ
for SAMPLE in "${ALL_SAMPLES[@]}"; do
    echo "Downloading $SAMPLE from NCBI..."
    wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/${SAMPLE}/${SAMPLE} \
        -O ${SAMPLE}.sra

    echo "Converting $SAMPLE to FASTQ..."
    fasterq-dump ${SAMPLE}.sra \
        --outdir . \
        --temp . \
        --threads 32 \
        --progress 


    echo "Cleaning up $SAMPLE .sra file..."
    rm -f ${SAMPLE}.sra

    echo "Done $SAMPLE"
done

#QC

fastp \
      --in1 "raw_data/${SAMPLE}_1.fastq" \
      --in2 "raw_data/${SAMPLE}_2.fastq" \
      --out1 "qc_reads/${SAMPLE}_1.fastq.gz" \
      --out2 "qc_reads/${SAMPLE}_2.fastq.gz" \
      --thread 16 \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --html "qc_reports/${SAMPLE}_fastp.html" \
      --json "qc_reports/${SAMPLE}_fastp.json"

    echo "Done $SAMPLE"
done


#Download the Kraken2 Database
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz --progress=bar:force 2>&1

#unzip to use
tar -xvzf kraken2_db/k2_standard_08_GB_20251015.tar.gz -C kraken2_db


SAMPLES=("SRR8146944" "SRR8146963" "SRR8146968" "SRR8146935" "SRR8146936" "SRR8146938")

mkdir -p kraken2_output

#Kraken2 Classification

kraken2 --db ~/scratch/BINF6110/kraken2_db --paired --threads 16 --confidence 0.15 raw_data/SRR8146938_1.fastq raw_data/SRR8146938_2.fastq --output kraken2_output/SRR8146938.kraken --report kraken2_output/SRR8146938.report

# ran this for each sample
#didnt use --memory-mapping because Compute Canada cluster has lots of RAM, so it should be fine without it.
#reports show read length of 149, so using that for Bracken.

SAMPLES=("SRR8146944" "SRR8146963" "SRR8146968" "SRR8146935" "SRR8146936" "SRR8146938")

mkdir -p bracken_output

#Bracken Abundance Re-estimation
for SRR in "${SAMPLES[@]}"; do
    echo "Bracken: $SRR"
    bracken -d "$HOME/scratch/BINF6110/kraken2_db" -i "kraken2_output/${SRR}.report" -o "bracken_output/${SRR}.bracken" -w "bracken_output/${SRR}_bracken.report" -r 150 -l S
    echo "Done Bracken: $SRR"
done


pip install --user h5py
pip install kraken-biom

#Convert to BIOM Format 
REPORTS=$(ls bracken_output/*_bracken.report | tr '\n' ' ')

kraken-biom $REPORTS

kraken-biom $REPORTS \
  --fmt json \
  -o combined_table.biom
