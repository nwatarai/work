#$1=nmr_of_chm $2=nmr_of_spc $3=nmr_of_reactions $4=nmr_of_environment $5=nmr_of_steps $6=stable(true ofr false) $7=suffix

PARAMFIX=${1}_${2}_${3}

cp header outfile.csv
cp header outfile_add.csv

a=0
while [ $a -ne 500 ]
do
FIX=$a
java starving $1 $2 $3 $4 $5 param_${PARAMFIX}.csv param_reaction_${PARAMFIX}.csv param_metabolism_${PARAMFIX}.csv param_chm_${PARAMFIX}.csv ${FIX} $6 ${a}
python extract_result.py result_${FIX}.csv 20
tail -n 1 result_${FIX}.csv > temp_tail
tail -n 11 result_${FIX}.csv.extract.csv | head -n 1 > temp_tail_add
cp outfile.csv temp_out
cp outfile_add.csv temp_out_add
cat temp_out temp_tail > outfile.csv
cat temp_out_add temp_tail_add > outfile_add.csv
a=`expr $a + 1`
rm result_${FIX}.csv
rm result_chm_${FIX}.csv
rm result_${FIX}.csv.extract.csv
done