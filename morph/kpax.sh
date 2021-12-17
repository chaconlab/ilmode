################################
# STRUCTURAL ALIGNMENT IN KPAX #
################################

base=DHV    # case basename
ali=_ali  # aligned files suffix
for((i=1;i<=30;i++))
do
 # Get case number
 num=`echo $i | awk '{printf "%02d", $1}'`

 # Get pdbs list for current case
 pdbs=`ls $base$num\_??????.pdb`

 echo "PDBs $pdbs"
  
 # Get the first one for reference
 ref=`echo $pdbs | awk '{print $1}'`

 # Get the remaining as for targets
 tar=`echo $pdbs | cut -d " " -f 2-`

 echo -e "PROCESSING $base$num\n\tReference: $ref  Targets: $tar"

 # Reference must have "suffix" for simplicity, but it is not moved at all
 cp $ref $base$num\_${ref:6:6}$ali.pdb

 # Aligning each target to the reference (common reference)
 for j in $tar
 do
  echo -e "\taligning $j (ref) to $ref (mobile)  --> ${j:0:12}$ali.pdb"
  /home/pablo/PROGS/kpax/bin/kpax $ref $j -nowrite -nosubdirs -pdb
  # echo "MV kpax_results/${j:0:12}_${ref:0:12}.pdb ${j:0:12}\$ali.pdb"
  mv kpax_results\/${j:0:12}\_${ref:0:12}.pdb ${j:0:12}$ali.pdb
  
 done

 # Aligning all pdbs aganist each other (non-redundant)
 for a in $pdbs
 do
  for b in $pdbs
  do
   if [[ "${a:0:12}" != "${b:0:12}" ]]  # Non-redundant pair-wise combinations
   then
   
    # Move KPAX output into our naming scheme
    a1=${a:0:12}$ali.pdb # aligned (mobile) id
    b1=${b:0:12}$ali.pdb # reference (fixed) id
    a2=kpax_results\/${b:0:12}$ali\_${a:0:12}$ali.pdb   
    b2=${b:0:12}\_${a:6:6}$ali.pdb       
    
    # Align in KPAX 
    echo -e "\taligning $a1 (ref) to $b1 (mobile) --> $b2\n\n"
    /home/pablo/PROGS/kpax/bin/kpax $a1 $b1 -nowrite -nosubdirs -pdb

    # Move KPAX output into our naming scheme
    echo -e "mv $a2 $b2"
    mv $a2 $b2
   fi
  done
 done
done
# Delete Kpax's shit
rm -r kpax_results


