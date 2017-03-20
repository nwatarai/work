#$1=nmr_of_chm $2=nmr_of_spc $3=nmr_of_reactions $4=nmr_of_environment $5=nmr_of_steps $6=stable(true ofr false) $7=suffix
FIX=${1}_${2}_${3}_${4}_${7}
PARAMFIX=${1}_${2}_${3}

mkdir ${FIX}
python make_parameters.py $1 $2 $3 ${PARAMFIX}
java starving $1 $2 $3 $4 $5 param_${PARAMFIX}.csv param_reaction_${PARAMFIX}.csv param_metabolism_${PARAMFIX}.csv param_chm_${PARAMFIX}.csv ${FIX} $6 $7

echo "analyzing results..."
python total_abundance.py result_${FIX}.csv 
python alpha_diversity.py result_${FIX}.csv
python extract_result.py result_${FIX}.csv 20
python extract_experiment_like.py result_${FIX}.csv
python delete_allZero_column.py result_${FIX}.csv.extract.csv
python delete_allZero_column.py result_${FIX}.csv.experiment_like.csv
python rank_abundance.py result_${FIX}.csv.extract.csv.del0.csv
python rank_abundance.py result_${FIX}.csv.experiment_like.csv.del0.csv
python rank_abundance.py result_chm_${FIX}.csv
python heatmap.py result_${FIX}.csv logging=True figx=100 figy=65 branch_color_threshold=0.5 clusterout=True
python heatmap.py result_chm_${FIX}.csv logging=True

python heatmap.py result_${FIX}.csv.extract.csv.del0.csv logging=True figx=20 figy=50 branch_color_threshold=0.5 clusterout=True
python heatmap.py result_${FIX}.csv.experiment_like.csv.del0.csv logging=True figx=20 figy=50 branch_color_threshold=0.5 clusterout=True
python reaction_to_sif.py param_reaction_${PARAMFIX}.csv
python reaction_abundance.py result_${FIX}.csv param_metabolism_${PARAMFIX}.csv

mv *${FIX}.* ${FIX}
cp *${PARAMFIX}.* ${FIX}