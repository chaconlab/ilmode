# ilmode
Normal Mode Analysis with constraints in internal coordinates


 
### sampling #### 

Move loop along a given mode till reach a give rmsd from the inital position: 
<pre>
ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _mod1F -m 2 -s 0 -a  1  --rmsd 3 --drmsd 0.25 --verb 0
ilmode 3hsz.pdb  81 93 --chain A -i 1 -o _mod1B -m 2 -s 0 -a -1  --rmsd 3 --drmsd 0.25 --verb 0    
../scripts/renum_tr.pl 3hsz_mod1F_traj.pdb 3hsz_mod1B_traj.pdb > 3hsz_mod1.pdb
</pre>

<video  width="320px" height="175px"  src="https://user-images.githubusercontent.com/19269061/141765600-329048b6-d184-4ae3-ba90-5e459d7ad189.mp4" autoplay="true" loop="true" controls="controls" >
  </video>
  
  

<pre>
for ((i=1;i<=17;i++)); 
do
     echo "processing mode $i";
     ../sbg/bin/ilmode DHV10_3irsC1_ali.pdb --start 66 --end 76 --chain C -m 2 -i $i -a  1 -r 2 -C 1 -s 0 --drmsd 0.25 -o _forward >> log;
     ../sbg/bin/ilmode DHV10_3irsC1_ali.pdb --start 66 --end 76 --chain C -m 2 -i $i -a -1 -r 2 -C 1 -s 0 --drmsd 0.25 -o _backward >> log;
     renum_tr.pl DHV10_3irsC1_ali_forward_traj.pdb   DHV10_3irsC1_ali_backward_traj.pdb > mode_$i.pdb     
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
 
 
