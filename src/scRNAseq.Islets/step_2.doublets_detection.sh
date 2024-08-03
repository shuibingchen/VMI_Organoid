# step_2.doublets_detection.sh
# --- process scRNAseq data from islets samples ---
# step 2: Identify doublets
# Author: Tuo Zhang
# Date: 8/1/2024
# 

# folders
workdir="."
sourcedir=${workdir}/source
infodir=${workdir}/info

# solo model
model=${sourcedir}/solo/model.json

# activate conda solo environment
source softwares/anaconda3/etc/profile.d/conda.sh
conda activate solo

# run solo
for sid in D1S1 D1S2 D2S1 D2S2 D2S3 D3S1 D3S2 D4S1 D4S2 D4S3
do
        # create output directory
        outdir=${infodir}/solo/${sid}
        logdir=${outdir}/logs
        if [ ! -d ${outdir} ]
        then
                mkdir ${outdir}
                mkdir ${logdir}
        fi

        solo --set-reproducible-seed 98 -o ${outdir} ${model} ${infodir}/${sid}.h5ad >${logdir}/solo.${sid}.log 2>&1
done

echo "Complete!"
