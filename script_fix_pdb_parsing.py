# 2018, Clementi Group
# READ IN ANTON .DCD FILE AND ONE PDB AS TOPOLOGY.
# PDB and dcd will be saved in a non corrupted format
# Fixing the mdtraj bug manually and saving the fixed trajectory to file.


import numpy as np
import mdtraj as md
import time
import copy
import sys
import os
import glob


def fix_topology(topology):
  
   new_top = topology.copy()
   residues = {}
   
   for chain in new_top.chains:
      for residue in chain.residues:    # here extra residues are read in because of non contiguous atoms
         resname = str(residue)
         #print residue, resname
         if resname in list(residues.keys()):
            residues[resname].append(residue)
         else:
            residues[resname] = [residue]
         #print residues    # at this point we have a list like {'VAL9': [val9, val9]}
              # all residues are repeated the number of times they are listed in the pdb in non
              # continuous fashion
      
   for resname in list(residues.keys()):
      #print resname
      fragments = residues[resname]
      #print fragments    # this is the list entry like [val9, val9,...] each residue is listed multiple times (see above)
                    # the number of times each residue is listed is the number of fragments (see next line)
      #print 'fragments 0   ' + str(fragments[0])
      
      if len(fragments) > 1:        
         main_fragment = fragments[0]
         #print 'main_frament._atoms', main_fragment._atoms
         
         new_atom_list = []
         new_atom_list += main_fragment._atoms
         
         for i in range(1,len(fragments)):
            fragment = fragments[i]        # selecting atoms per fragment
            #print 'fragment   ' + str(i) + '   atoms=   ' + str(fragment._atoms)
            
            for atom in fragment.atoms:
               atom.residue = main_fragment
            
            new_atom_list += fragment._atoms   # appending the atom list to get the 'new_list' for that residue
            fragment._atoms = []               # setting no atom in this fragment
            fragment.chain = main_fragment.chain
         main_fragment._atoms = new_atom_list
         #print 'fine', main_fragment._atoms  
   
   # at this point, the order of atoms has been rearranged in such a way that all original fragments but 'main_fragment' 
   # are empty
   
   return new_top

def get_num(x):
    return ''.join(ele for ele in x if ele.isdigit())

def fix_traj(traj):
  time0 = time.time()
  new_traj = copy.deepcopy(traj)
  topology = new_traj.topology 
  
  new_top = fix_topology(topology)
  new_traj.topology = new_top 
  
  new_atom_sequence = [a for a in new_top.atoms]
  new_index_sequence = [a.index for a in new_top.atoms]
  
  for i in range(0, np.shape(traj.xyz)[0]):
    new_traj.xyz[i] = new_traj.xyz[i][new_index_sequence,:]
  
  for i in range(0, len(new_index_sequence)):
    new_atom_sequence[i].index = i
  
  new_top=new_traj.topology
  
  new_atom_names = [get_num(str(a.residue)).zfill(3)+'_'+str(a.name) for a in new_top.atoms]
  b=np.array(new_atom_names)
  new_index=np.argsort(b)
  #new_atom_names = [str(a.residue)+'-'+str(a.name) for a in topology.atoms]
  #print new_atom_names
  #new_index=sorted(range(len(new_atom_names)), key=new_atom_names.__getitem__)
  #print new_index
  #print new_atom_names
  #print b[new_index]
  
  #[new_atom_names[i] for i in new_index]
  #print [new_atom_sequence[i] for i in new_index]
  
  for i in range(0, np.shape(traj.xyz)[0]):
    new_traj.xyz[i] = new_traj.xyz[i][new_index,:]
  
  new_index_sequence = [a.index for a in new_top.atoms]
  #for i in range(0, len(new_index)):
  #  new_atom_sequence[i].index = i#new_index[i]
  
  tmp=copy.deepcopy([a for a in new_top.atoms])
  for i in range(0, len(new_index)):
    new_atom_sequence[i].element = tmp[new_index[i]].element
    new_atom_sequence[i].name = tmp[new_index[i]].name
  
  time1 = time.time()
  #print(time1 - time0)
  return new_traj


#----------------------------------------------------------------------------------------------
# read all trajectory data and fix them
print('starting fixing pdb structure')

#define all files input and output

file_input='system-protein.pdb'
file_output_fixed_pdb='system-protein-fixed.pdb'
file_output_fixed_gro='system-protein-fixed.gro'
file_output_fixed_netcdf='system-protein-fixed.netcdf'
file_output_fixed_mdcrd='system-protein-fixed.mdcrd'
#load
topol = md.load(file_input)
#convert into fixed strcuture
topol_new = fix_traj(topol)

#save as pdb
topol_new.save_pdb(file_output_fixed_pdb)

#reload pdb to correctly write gro file

topol2_new = md.load(file_output_fixed_pdb)

#save as GRO format
topol2_new.save_gro(file_output_fixed_gro)

# save as amber format
topol2_new.save_netcdf(file_output_fixed_netcdf)
# save as amber format
topol2_new.save_mdcrd(file_output_fixed_mdcrd)


#In case there are trajectories to fixed too
path_input = 'xxx'
#output format determined by filename extension
path_output_fixed = 'yyy'

traj = md.load(path_input, top = file_input)
# fixing the trajectory...
traj_new = fix_traj(traj)    
#saving
traj_new.save(path_output_fixed)

