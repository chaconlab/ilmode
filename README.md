# ilmode
Normal Mode Analysis with constraints in internal coordinates


 
### sampling #### 
<pre>
ilmode_gcc DHV17_3hszA1_ali.pdb --start 81 --end 93  --chain A -m 2 -s 2 --skip_missingatoms -a 1 -C 1 --ns 100 --drmsd 0.25 --rmsd 12 -o _sampling
vmd -m DHV17_3hszA1_ali.pdb DHV17_3hszA1_ali_sampling_traj.pdb
</pre>


<pre>
for ((i=1;i<=17;i++)); 
do
     echo "processing mode $i";
     ../sbg/bin/ilmode DHV10_3irsC1_ali.pdb --start 66 --end 76 --chain C -m 2 -i $i -a  1 -r 2 -C 1 -s 0 --drmsd 0.25 -o _forward >> log;
     ../sbg/bin/ilmode DHV10_3irsC1_ali.pdb --start 66 --end 76 --chain C -m 2 -i $i -a -1 -r 2 -C 1 -s 0 --drmsd 0.25 -o _backward >> log;
     renum_tr.pl DHV10_3irsC1_ali_forward_traj.pdb   DHV10_3irsC1_ali_backward_traj.pdb > mode_$i.pdb     
done
 </pre>
 
 <video src="https://user-images.githubusercontent.com/19269061/141276068-1f77ef8e-a67f-464e-8a9a-526dca14c393.mp4" data-canonical-src="https://user-images.githubusercontent.com/19269061/141276068-1f77ef8e-a67f-464e-8a9a-526dca14c393.mp4" controls="controls" muted="muted" class="d-block rounded-bottom-2 width-fit" style="max-height:640px;">

  </video>
  


### morphing ###

<pre>
ilmode DHV17_3hszA1_ali.pdb --start 81 --end 93 --chain A -t DHV17_3ht0A2_ali.pdb -m 2 --skip_missingatoms -a 1 -C 1 --ns 5000 --flanks 1 --aliflank
s --drmsd 0.25 -o  _morph
vmd -m DHV17_3hszA1_ali.pdb DHV17_3ht0A2_ali.pdb  DHV17_3hszA1_ali_4m45A1_af1_traj.pdb 
</pre>
 
 
