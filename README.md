# ilmode
Normal Mode Analysis with constraints in internal coordinates


 
### Sampling #### 

Move a loop along a given mode direction till reach a give rmsd from the inital conformation: 
<pre>
<<<<<<< HEAD
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

=======
../sbg/bin/ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _mod1F -m 2 -s 0 -a  1  --rmsd 3.0 --drmsd 0.25 --verb 0
../sbg/bin/ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _mod1B -m 2 -s 0 -a -1  --rmsd 3.0 --drmsd 0.25 --verb 0    
../scripts/renum_tr.pl 3hsz_mod1F_traj.pdb 3hsz_mod1B_traj.pdb > 3hsz_mod1.pdb
>>>>>>> 317a3e873fa021eed106b2d60bb77434b8616fde
</pre>
The first command move forward (-a 1) the 81-93 loop until the rmsd from the inital position is > 3.0 Å (--rmsd 3.0). Every 0.25Å (--drmsd 0.25) away from the intial pose the moved loop coordinates are saved in the 3hsz_mod1F_traj.pdb trajectory file.  The second command do the same but backwards, and the final generate a forward-backward trajectory like this:  
 
<video  width="320px" height="175px"  src="https://user-images.githubusercontent.com/19269061/141765600-329048b6-d184-4ae3-ba90-5e459d7ad189.mp4" autoplay="true" loop="true" controls="controls" >
  </video>

<<<<<<< HEAD
../sbg/bin/ilmode 3hsz.pdb --start 81 --end 93 --chain A -m 2 -s 0 --skip_missingatoms -a 1 --ns 100 --rmsd 6 -o _sampling
=======
Note how the closure is fully mantained despite of the large motion. Please, activate the loop option (right click on play icon) for a better display. 
Here with a different loop.
>>>>>>> 317a3e873fa021eed106b2d60bb77434b8616fde

Here it is another example with a different loop in where we compute all the modes: 
<pre>
for ((i=1;i<=17;i++)); 
do
     echo "processing mode $i";
<<<<<<< HEAD
     ../sbg/bin/ilmode 3irs.pdb --start 66 --end 76 --chain C -m 2 -i $i -a  1 -r 2  --drmsd 0.25 -o _f >> log;
     ../sbg/bin/ilmode 3irs.pdb --start 66 --end 76 --chain C -m 2 -i $i -a -1 -r 2  --drmsd 0.25 -o _b >> log;
     ../scripts/renum_tr.pl 3irs_f_traj.pdb  3irs_b_traj.pdb > mode_$i.pdb     
=======
     ../sbg/bin/ilmode 3irs.pdb 66 76 --chain C -m 2 -i $i -a  1 -r 2 -C 1 -s 0 --drmsd 0.25 -o F >> log;
     ../sbg/bin/ilmode 3irs.pdb 66 76 --chain C -m 2 -i $i -a -1 -r 2 -C 1 -s 0 --drmsd 0.25 -o B >> log;
     renum_tr.pl  3irsF_traj.pdb  3irsB_traj.pdb > mode_$i.pdb     
>>>>>>> 317a3e873fa021eed106b2d60bb77434b8616fde
done
 </pre>
 And here we display all the results together:   
 <video  width="320px" height="175px"  src="https://user-images.githubusercontent.com/19269061/141282242-ac69849d-3ceb-4241-8f11-fcdb0ab5c0a4.mp4" autoplay="true" loop="true" controls="controls" >
  </video>

Alternatively to a single mode motion, you can move in the direcction defined by a random contribution of all the modes:   
<pre>
../sbg/bin/ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _allF -m 2 -s 1 -a  1  --rmsd 3.0 --drmsd 0.25 --verb 0
../sbg/bin/ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _allB -m 2 -s 1 -a -1  --rmsd 3.0 --drmsd 0.25 --verb 0 
../scripts/renum_tr.pl 3hsz_allF_traj.pdb 3hsz_allB_traj.pdb > 3hsz_all.pdb
</pre>

### Morphing ###

<pre>
ilmode DHV17_3hszA1_ali.pdb --start 81 --end 93 --chain A -t DHV17_3ht0A2_ali.pdb -m 2 --skip_missingatoms -a 1 -C 1 --ns 5000 --flanks 1 --aliflank
s --drmsd 0.25 -o  _morph
vmd -m DHV17_3hszA1_ali.pdb DHV17_3ht0A2_ali.pdb  DHV17_3hszA1_ali_4m45A1_af1_traj.pdb 
</pre>
 
 
