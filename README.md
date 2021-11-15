# ilmode
Normal Mode Analysis with constraints in internal coordinates


 
### sampling #### 
<pre>
../sbg/bin/ilmode 3hsz.pdb --start 81 --end 93 --chain A -m 2 -s 2 --skip_missingatoms -a 1 --ns 100 --rmsd 6 -o _sampling


-a Maximum angular increment 
-s   Sampling strategy
      0 "i"-th mode-following sampling, 
      1  Single-NMA Monte-Carlo, 
      2  Radial-Mode-Following sampling 
-i  Index of the selected mode
       
-ns Number of samples for selected strategy

-m 0=CA, 1=C5, 2=Heavy-Atom (default=2). 

ilmode 3hsz.pdb --start 81 --end 93 --chain A -m 3 -s 2  -a 3  --rmsd 2 -o _sampling --verb 0 --ns 1000

ilmode 3hsz.pdb --start 81 --end 93 --chain A -m 3 -s 1  -a 3  --rmsd 2 -o _sampling --verb 0 --ns 1000


move modes 

ilmode temp.pdb  81 93 --chain A -i 1 -o _mod1F -m 2 -s 0 -a  1  --rmsd 3 --drmsd 0.25 --verb 0
ilmode temp.pdb  81 93 --chain A -i 1 -o _mod1B -m 2 -s 0 -a -1  --rmsd 3 --drmsd 0.25 --verb 0    
../scripts/renum_tr.pl 3hsz_mod1F_traj.pdb 3hsz_mod1B_traj.pdb > 3hsz_mod1.pdb



moving a helix...

ilmode 3hsz.pdb --start 63 --end 81 --chain A -m 3 -s 2  -a 1  --rmsd 2 -o _sampling --verb 0 --ns 100

</pre>

../sbg/bin/ilmode 3hsz.pdb --start 81 --end 93 --chain A -m 2 -s 0 --skip_missingatoms -a 1 --ns 100 --rmsd 6 -o _sampling

<pre>
for ((i=1;i<=17;i++)); 
do
     echo "processing mode $i";
     ../sbg/bin/ilmode 3irs.pdb --start 66 --end 76 --chain C -m 2 -i $i -a  1 -r 2  --drmsd 0.25 -o _f >> log;
     ../sbg/bin/ilmode 3irs.pdb --start 66 --end 76 --chain C -m 2 -i $i -a -1 -r 2  --drmsd 0.25 -o _b >> log;
     ../scripts/renum_tr.pl 3irs_f_traj.pdb  3irs_b_traj.pdb > mode_$i.pdb     
done
 </pre>
 
 <video  width="320px" height="175px"  src="https://user-images.githubusercontent.com/19269061/141282242-ac69849d-3ceb-4241-8f11-fcdb0ab5c0a4.mp4" autoplay="true" loop="true" controls="controls" >
  </video>


### morphing ###

<pre>
ilmode DHV17_3hszA1_ali.pdb --start 81 --end 93 --chain A -t DHV17_3ht0A2_ali.pdb -m 2 --skip_missingatoms -a 1 -C 1 --ns 5000 --flanks 1 --aliflank
s --drmsd 0.25 -o  _morph
vmd -m DHV17_3hszA1_ali.pdb DHV17_3ht0A2_ali.pdb  DHV17_3hszA1_ali_4m45A1_af1_traj.pdb 
</pre>
 
 
