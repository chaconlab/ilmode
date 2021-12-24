########################################################
# MORPHING IN ilmode
########################################################
prog="ilmode"
ntoff=0  # N-terminal offset from "start" residue in input file
ctoff=0  # C-terminal offset from "end" residue in input file

suff="_aF0"  # general output suffix
opts="-m 2 --skip_missingatoms -a 1 --drmsd 0.25"

suff="_aF0x"  # general output suffix
opts="-m 2 -x --skip_missingatoms -a 1 --drmsd 0.25"

#suff="_af0"  # general output suffix
#opts="-m 3 --skip_missingatoms -a 1 --drmsd 0.25"

stat="stat$suff.txt"  # general output statistics file
echo "# Morphing          RmsdI RmsdF  RmsdINCAC  RmsdFNCAC  N_res N_dh  RmsdD delta0 Motion" > $stat  # dump statistics file header
echo "# ---------------------------------------------------------------------------------------------------------" >> $stat  # dump statistics file header
input=DeaneHV.txt  # RCD input file with loop indices and chain ids, e.g. DeaneHV.txt
base=DHV    # case basename
ali=_ali  # aligned files suffix





time for((i=1;i<=30;i++))
do
 # Get case number formated
 num=`echo $i | awk '{printf "%02d", $1}'`

 # Get pdbs list for current case
 pdbs=`ls $base$num\_??????.pdb`

 # Morphing all PDBs against each other (all-against-all)
 for a in $pdbs
 do
  for b in $pdbs
  do
   
   if [[ "$a" != "$b" ]]  # Non-self combinations
   then

    # Get loop indices and chain ids from RCD input file
    name=`basename $a $ali.pdb`
    row=`grep $name $input`
    start=`echo $row | awk -v ntoff=$ntoff '{print $2+ntoff}'`
    end=`echo $row | awk -v ctoff=$ctoff  '{print $3+ctoff}'`
    chain=`echo $row | awk '{print $4}'`

    # Get target ID
    bb=`basename $b .pdb`
    bid=${bb:6:6}
    a1=${a:0:12}$ali.pdb # aligned (mobile) id
    b1=${b:0:12}\_${a:6:6}$ali.pdb    
    
    
    # Perform morphing
    echo -e "MORPHING $prog $a1  $start $end --chain $chain -t $b1 $opts -o _$bid$suff"
    $prog $a1 $start $end --chain $chain -t $b1 $opts -o _$bid$suff 

    # Collect STATISTICS
    name2=`basename $a1 .pdb`
    log=$name2\_$bid$suff\_traj.log
    log2=$name2\_$bid$suff
    line=`grep "ilmode> Morphing:  Initial_RMSD" $log`
    rmsd0=`echo $line | awk '{print $4}'`
    rmsd0NCAC=`echo $line | awk '{print $6}'`
    rmsdF=`echo $line | awk '{print $8}'`
    rmsdNCAC=`echo $line | awk '{print $10}'`
    rmsdD=`echo $line | awk '{print $12}'`
    delta0=`echo $line | awk '{print $14}'`
    motion=`echo $line | awk '{print $16}'`
    line=`grep "ilmode> Residues:" $log`
    nres=`echo $line | awk '{print $3}'`
    ndh=`echo $line | awk '{print $5}'`
    echo "$log2  $rmsd0 $rmsdF $rmsd0NCAC $rmsdNCAC $nres $ndh $delta0 $motion" >> $stat  # dump statistics (append)
   fi
  done
 done
done


# Computing averages from statistics
stat0=`basename $stat .txt`
stat2=$stat0'A0S.txt'
head -n 1 $stat > $stat2
statext.pl $stat 1 10 1 all 1 >> $stat2
echo -e "Results full set --> $stat2"

# Averages for RmsdI > 1
stat3=$stat0'A1'.txt
awk '{if($4>1) print $0}' $stat > $stat3
stat3b=$stat0'A1S.txt'
head -n 1 $stat > $stat3b
statext.pl $stat3 1 10 1 all 1 >> $stat3b
echo -e "Results RmsdI>1 --> $stat3b"

# Averages for RmsdI > 2
stat4=$stat0'A2'.txt
awk '{if($4>2) print $0}' $stat > $stat4
stat4b=$stat0'A2S'.txt
head -n 1 $stat > $stat4b
statext.pl $stat4 1 10 1 all 1 >> $stat4b
echo -e "Results RmsdI>2 --> $stat4b"



