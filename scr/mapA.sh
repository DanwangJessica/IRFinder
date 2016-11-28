#Map IR to all the genomes which is shown in a txt file
#$1 is the IR Seq to be maped
#$2 is the faNamesPre.txt file
#$3 is the path for the genomes
#$4 is the path for blast output

declare -i totGenome="$(grep -n '^' $2|tail -n 1|cut -d ':' -f1)"

export PATH="$PATH:$HOME/blast-2.2.26/bin"

for (( i=1; i<$totGenome+1; i=i+1 ))
do
    namePre=$(sed -n ''$i'p' $2) 
    
    faName=$(find $3 -name "*$namePre*"|sed -n '1p')
    
    bl2seq -i $1 -j $faName -p blastn -e 1e-50 -o $4/"$i"_$1map.txt
done
