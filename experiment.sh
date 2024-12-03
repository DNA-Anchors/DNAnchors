#!/bin/sh
set -x

# Set up directories

BASE_DIR='./data'  # data directory, update this line correspondingly

CODE_DIR='.'  # code directory, update this line correspondingly

# The DNA based data storage pipeline consists of the following steps:
# 1) Encoding data into DNA strands
# 2) Synthesizing, storing and sequencing DNA strands - wetlab activities which introduce errors
# 3) Clustering sequenced reads
# 4) Reconstructing original DNA strands from clusters of reads 
# 5) Decoding data from reconstructed strands
# 
# We will be executing the entire pipeline here:

# Pick configuration file for Input Data

file_name1='baba_smol_woAnchors' # Directory in which we will encode and decode without Anchors
file_name2='baba_smol_wAnchors' # Directory in which we will encode and decode with Anchors
config_file='baba_smol'

skipRS='0' # '0' if you want to use Reed Solomon Code for redundancy, '1' if not

# Probability of error of each type (insertions, deletions, substitutions)
P='0.02'

# Number of reads per DNA strand
coverage='10'

# 1) Encoding data into DNA strands

# First we encode file 'baba_smol' into strands with baseline encoding

cd $CODE_DIR/1-encoding-decoding
python3 codec.py ../$BASE_DIR configs/${config_file}_woAnchors.cfg 0 $file_name1 $skipRS 
cd ..

# Next we encode file 'bab_smol' into strands with encoding for anchors

cd $CODE_DIR/1-encoding-decoding
python3 codec.py ../$BASE_DIR configs/${config_file}_wAnchors.cfg 0 $file_name2 $skipRS 

# Next we add anchors to the encoded strands

python3 generate_anchors.py
python3 generate_anchored_strands.py
cd ..


# 2) Simulating wetlab activities which introduce errors - Synthesizing, storing and sequencing DNA strands - 


python3 $CODE_DIR/2-simulating_wetlab/naive/noise.py --N ${coverage}  --subs $P --dels $P --inss $P  --i $BASE_DIR/output/${file_name1}/EncodedStrands.txt  --o $BASE_DIR/output/${file_name1}_P${P}_N${coverage}/UnderlyingClusters.txt
python3 $CODE_DIR/2-simulating_wetlab/naive/shuffle.py $BASE_DIR/output/${file_name1}_P${P}_N${coverage}/UnderlyingClusters.txt  $BASE_DIR/output/${file_name1}_P${P}_N${coverage}/NoisyStrands.txt



python3 $CODE_DIR/2-simulating_wetlab/naive/noise.py --N ${coverage}  --subs $P --dels $P --inss $P  --i $BASE_DIR/output/${file_name2}/EncodedStrands_wAnchors.txt  --o $BASE_DIR/output/${file_name2}_P${P}_N${coverage}/UnderlyingClusters.txt
python3 $CODE_DIR/2-simulating_wetlab/naive/shuffle.py $BASE_DIR/output/${file_name2}_P${P}_N${coverage}/UnderlyingClusters.txt  $BASE_DIR/output/${file_name2}_P${P}_N${coverage}/NoisyStrands.txt


# 3) Clustering sequenced reads
 

cd $CODE_DIR/3-clustering
cp ../$BASE_DIR/output/${file_name1}_P${P}_N${coverage}/NoisyStrands.txt input/.
cp ../$BASE_DIR/output/${file_name1}_P${P}_N${coverage}/UnderlyingClusters.txt input/.
make run

cp ./output/ClusteredStrands.txt  ../$BASE_DIR/output/${file_name1}_P${P}_N${coverage}/.
cd ..

cd $CODE_DIR/3-clustering
cp ../$BASE_DIR/output/${file_name2}_P${P}_N${coverage}/NoisyStrands.txt input/.
cp ../$BASE_DIR/output/${file_name2}_P${P}_N${coverage}/UnderlyingClusters.txt input/.
make run

cp ./output/ClusteredStrands.txt  ../$BASE_DIR/output/${file_name2}_P${P}_N${coverage}/.
cd ..

# 4) Reconstructing original DNA strands from clusters of reads 

# Recon without anchors

python3 $CODE_DIR/4-reconstruction/recon.py --L 176 --i $BASE_DIR/output/${file_name1}_P${P}_N${coverage}/ClusteredStrands.txt  --o  $BASE_DIR/output/${file_name1}_P${P}_N${coverage}/ReconstructedStrands.txt --path /dev/shm --coverage $coverage


# Recon with anchors

python3 $CODE_DIR/4-reconstruction/reconstruct_with_anchor.py --anchor GTCTGAAC --L 176 --i $BASE_DIR/output/${file_name2}_P${P}_N${coverage}/ClusteredStrands.txt  --o  $BASE_DIR/output/${file_name2}_P${P}_N${coverage}/ReconstructedStrands.txt --path /dev/shm --coverage $coverage

# 5) Decoding data from reconstructed strands

# Decoding files without anchors

cd $CODE_DIR/1-encoding-decoding
python3 codec.py ../$BASE_DIR configs/${config_file}_woAnchors.cfg 1 ${file_name1}_P${P}_N${coverage} $skipRS 
cd ..

# Decoding file with anchors

cd $CODE_DIR/1-encoding-decoding
python3 codec.py ../$BASE_DIR configs/${config_file}_wAnchors.cfg 1 ${file_name2}_P${P}_N${coverage} $skipRS 
cd ..
