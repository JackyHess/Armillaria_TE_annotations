#!/bin/bash

### Step 1: Get validated TE consensus seqs

TEannot.py -P Abisporus -C TEannot.cfg -S 1 &> step1.txt

TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a CEN  &> step2_CEN.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a RM  &> step2_RM.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a BLR  &> step2_BLR.txt

TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a BLR -r  &> step2_BLRr.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a CEN -r &> step2_CENr.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a RM -r &> step2_RMr.txt

TEannot.py -P Abisporus -C TEannot.cfg -S 3 -c BLR+RM+CEN &> step3.txt

TEannot.py -P Abisporus -C TEannot.cfg -S 7 &> step7.txt

# Extract validated TEs

export REPET_HOST="localhost"
export REPET_USER=""
export REPET_PW=""
export REPET_DB="abisporus"

# update genome size for each species
PostAnalyzeTELib.py -a 3 -p Abisporus_chr_allTEs_nr_join_path -s Abisporus_refTEs_seq -g 33100000
#
GetSpecificTELibAccordingToAnnotation.py -i Abisporus_chr_allTEs_nr_join_path.annotStatsPerTE.tab -t Abisporus_refTEs_seq -v 1
#
rm -r Abisporus_db/
rm -r Abisporus_TEdetect*
mv Abisporus_chr_allTEs_nr_join_path.annotStatsPerTE_FullLengthFrag.fa Abisporus_refTEs.fa
rm *.txt
rm Abisporus_chr_allTEs_nr_*
#
#### Step 2: Run full pipeline
#
TEannot.py -P Abisporus -C TEannot.cfg -S 1 &> step1.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a CEN  &> step2_CEN.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a RM  &> step2_RM.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a BLR  &> step2_BLR.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a BLR -r  &> step2_BLRr.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a CEN -r &> step2_CENr.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 2 -a RM -r &> step2_RMr.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 3 -c BLR+RM+CEN &> step3.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 4 -s TRF &> step4_TRF.txt
TEannot.py -P Abisporus -C TEannot.cfg -S 4 -s Mreps &> step4_Mreps.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 5 &> step5.txt
#
ln -s /home/jacky/Software/REPET_linux-x64-2.5/repeat_libraries/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa
#
TEannot.py -P Abisporus -C TEannot.cfg -S 6 -b blastx &> step6_blastx.txt
#
TEannot.py -P Abisporus -C TEannot.cfg -S 7 &> step7.txt &
#
TEannot.py -P Abisporus -C TEannot.cfg -S 8 -o GFF3 &> step8.txt &
#
