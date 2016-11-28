#!/bin/bash
#$1 is the scaffold_map file,like mgr_macro.txt
#$2 is the block.txt file
#$3 is the maf file
#$4 is the absolute genome path
#$5 is the absolute path for reportA
#$6 is the Maximum distance cutoff

export PATH="$PATH:$HOME/GRIMM/GRIMM-2.01"
export PATH="$PATH:$HOME/blast-2.2.26/bin"

declare -i totGenomes="$(grep '^a' $3|sed '2,$d'|cut -d '=' -f4)"
declare -i RANGE=5000 #the search range of A
declare -i revCut=5 #the threshold for reversd no. of blocks
declare -i ALCut=500 #the length of A threshold
declare -i ASimCut=90 #the similarity threshold of A and -A
declare -i curDis=$6

declare -i distance=0
declare -i LBlk=0
declare -i RBlk=0

declare -i count=0
declare -i totRev=0
declare -i countBoth=0

#grimm -f $1 -C >$5/reportA.txt
echo "A/-A findings between genomes with distance less than $6" >$5/reportA.txt


declare -i endLine=2+$totGenomes
grep -A$totGenomes 'mult='$totGenomes'' $3|sed ''$endLine',$d'|sed '1d'|cut -d ' ' -f2|cut -d '.' -f1 >faNamesPre.txt


for (( i=1; i<$totGenomes+1; i=i+1 ))
do
    for (( j=$i+1; j<$totGenomes+1; j=j+1 ))
    do
        sh deleTransp.sh $1 $i $j #delete the transposition and block interchange cases in genome i and genome j
        grimm -f mgr_macro2.txt -C -g 1,2 >$5/sn"$i"To"$j".sn #reversal senario from genome i to j
        
        distance="$(grep 'Reversal Distance' $5/sn"$i"To"$j".sn|cut -d$'\t' -f2)"
    
        if [ $distance -gt $curDis ]; then
            continue
            #echo ">Genome $i to $j: $distance steps : Distance is too large" >>$5/reportA.txt
        else
            echo ">Genome $i to $j: $distance steps" >>$5/reportA.txt
            for (( k=1; k<$distance+1; k=k+1 )) #check reversal in each steps
            do
                LBlk="$(grep 'Step '$k': ' $5/sn"$i"To"$j".sn|cut -d '[' -f2|sed 's/\']'.*$//g')"
                RBlk="$(grep 'Step '$k': ' $5/sn"$i"To"$j".sn|cut -d '[' -f3|sed 's/\']'.*$//g')"
                
                if [ $LBlk -gt 0 ]; then
                    LBlk=$(sh getleftBlk.sh $LBlk)
                else
                    LBlk=$LBlk*-1
                    LBlk=$(sh getrightBlk.sh $LBlk)
                    LBlk=$LBlk*-1
                fi

                if [ $RBlk -gt 0 ]; then
                    RBlk=$(sh getrightBlk.sh $RBlk)
                else
                    RBlk=$RBlk*-1
                    RBlk=$(sh getleftBlk.sh $RBlk)
                    RBlk=$RBlk*-1
                fi
                                 
                
                totRev=$totRev+1
                sh getleftSq.sh $2 $i $LBlk $RANGE faNamesPre.txt $4 >$5/sLeft$LBlk.fa 
                sh getRightSq.sh $2 $i $RBlk $RANGE faNamesPre.txt $4 >$5/sRight$RBlk.fa
                bl2seq -i $5/sLeft$LBlk.fa -j $5/sRight$RBlk.fa -p blastn -e 1e-50 -o sBlast.txt
                searchSA=$(sh findA.sh $LBlk $RBlk $2 $i sBlast.txt $ALCut $ASimCut)
                    
                               
                sh getleftSq.sh $2 $j $LBlk $RANGE faNamesPre.txt $4 >$5/dLeft$LBlk.fa
                sh getRightSq.sh $2 $j $RBlk $RANGE faNamesPre.txt $4 >$5/dRight$RBlk.fa
                bl2seq -i $5/dLeft$LBlk.fa -j $5/dRight$RBlk.fa -p blastn -e 1e-50 -o dBlast.txt
                searchDA=$(sh findA.sh $LBlk $RBlk $2 $j dBlast.txt $ALCut $ASimCut)
                

                echo "Step $k: $LBlk through $RBlk Reversal: Source:$searchSA Destination:$searchDA" >>$5/reportA.txt

                if [[ $searchSA == *"Found"* ]] || [[ $searchDA == *"Found"* ]]; then
                    count=$count+1                   
                fi

                if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
                    countBoth=$countBoth+1
                fi

                
                #else
                    #echo "Step $k: $LBlk through $RBlk Reversal ($revGeneNum blocks): Reversed Segment Too Short" >>$5/reportA.txt
                #fi

            done
        fi

    done
done

echo "Total $totRev large reversals, and $count of them have A/-A. $countBoth of A/-A occurs in both genomes." >>$5/reportA.txt
rm *.fa
rm *.sn