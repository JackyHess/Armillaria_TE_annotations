#!/bin/bash

TEdenovo.py -P Abisporus -S 1 -C TEdenovo.cfg >& step1.txt

TEdenovo.py -P Abisporus -S 2 -C TEdenovo.cfg >& step2.txt

TEdenovo.py -P Abisporus -S 2 -C TEdenovo.cfg --struct >& step2_struct.txt

# The next three can be run in parallel although they are sufficiently fast to not worry too much about
# submitting them at the same time

TEdenovo.py -P Abisporus -S 3 -s Blaster -c Grouper -C TEdenovo.cfg >& step3_grouper.txt

TEdenovo.py -P Abisporus -S 3 -s Blaster -c Recon -C TEdenovo.cfg >& step3_recon.txt

TEdenovo.py -P Abisporus -S 3 -s Blaster -c Piler -C TEdenovo.cfg >& step3_piler.txt

TEdenovo.py -P Abisporus -S 3 --struct -C TEdenovo.cfg >& step3_struct.txt

TEdenovo.py -P Abisporus -S 4 -s Blaster -c Grouper -m Map -C TEdenovo.cfg >& step4_grouper.txt

TEdenovo.py -P Abisporus -S 4 -s Blaster -c Piler -m Map -C TEdenovo.cfg >& step4_piler.txt

TEdenovo.py -P Abisporus -S 4 -s Blaster -c Recon -m Map -C TEdenovo.cfg >& step4_recon.txt

TEdenovo.py -P Abisporus -S 4 --struct -m Map -C TEdenovo.cfg >& step4_struct.txt

#LaunchRepeatScout.py -i Abisporus.fa -c -v 2
# For some reason RepeatScout did not work. Abandoned.

ln -s /home/jacky/Software/REPET_linux-x64-2.5/repeat_libraries/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa .

ln -s /home/jacky/Software/REPET_linux-x64-2.5/repeat_libraries/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa .

ln -s /home/jacky/Software/REPET_linux-x64-2.5/repeat_libraries/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm .

TEdenovo.py -P Abisporus  -C TEdenovo.cfg -S 5 -s Blaster -c GrpRecPil -m Map --struct >&step5.txt

TEdenovo.py -P Abisporus  -C TEdenovo.cfg -S 6 -s Blaster -c GrpRecPil -m Map --struct >&step6.txt 

#Final filtering step - remove SSRs, NoCat elements with less than 10 repeats (may be putative gene families, etc) and chimeric elements. This is quite stringent, but for this analysis we'd rather have a high quality, conservative set of TEs for downstream analysis

TEdenovo.py -P Abisporus  -C TEdenovo.cfg -S 7 -s Blaster -c GrpRecPil -m Map --struct >&step7.txt 

## Step 8 - clustering by TE family will be performed on the set of TEs consensus sequences from all species in the analysis

# cat ../REPET_pipeline/*/*Blaster_GrpRecPil_Struct_Map_TEclassif_Filtered/*denovoLibTEs_filtered.fa > Pan_species_TElibrary.fa
# TEdenovo.py -P PanSpecies -C TEdenovo.cfg -S 8 --struct -s Blaster -c GrpRecPil -m Map -f MCL

