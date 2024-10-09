#!/usr/bin/env bash

#SBATCH -J 20241008_calculating_n_reads_per_promoter
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 3500
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Set variables and directories
longLine="--------------------"
THREADS=13

mamba init
mamba activate py3
module load samtools
module load bedtools2

cd /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/H3K4me3
touch 'snCUTRUN_90_cells_reads_per_promoter.txt'

# Define the path to your BED file

bed_file_path='promoters_hg38.bed'

# Loop over each line in the BED file
while IFS= read -r line; do
    
# Process the line (for now, we just print it)
    echo "$line"  > line.bed

for i in `cat 'filenames_HN120Pri.txt'`

do
    #Intersect bed files with bedtools
    bedtools intersect -u -a "$i"'.bed' -b line.bed  > temp.bed
    line_count=$(wc -l < temp.bed)
    echo "$line_count" > "$i"'_reads_in_promoters.txt'
    paste line.bed "$i"'_reads_in_promoters.txt' > combined.bed
    cat combined.bed >> 'snCUTRUN_90_cells_reads_per_promoter.txt'

done


done < "$bed_file_path"





