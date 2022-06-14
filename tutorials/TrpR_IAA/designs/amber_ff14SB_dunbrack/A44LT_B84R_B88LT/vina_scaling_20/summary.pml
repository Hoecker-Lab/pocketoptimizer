load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/scaffold/amber_ff14SB/scaffold.pdb, WT_scaffold
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/ligand/amber_ff14SB/ligand.mol2, WT_ligand
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/0/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/0/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/1/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/1/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/2/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/2/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/3/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/3/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/4/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/4/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/5/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/5/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/6/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/6/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/7/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/7/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/8/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/8/receptor.pdb, designs
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/9/ligand.mol2, ligand_poses
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/9/receptor.pdb, designs
bg_color white
hide lines
cartoon loop
set stick_radius, 0.15
show stick, ligand_poses
select old_pocket, ((chain A and resi 44 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 84 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 88 and not (name C or name N or name O or name H or name HA))) and WT_scaffold
select new_pocket, ((chain A and resi 44 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 84 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 88 and not (name C or name N or name O or name H or name HA))) and designs
select mutable, ((chain A and resi 44) or (chain B and resi 88)) and designs
label mutable and name ca and design, "%s - %s" %(resn,resi)
show stick, old_pocket
show stick, new_pocket
show lines, WT_scaffold and not old_pocket
hide lines, WT_scaffold and name c+o+n
color gray90, designs and elem C
color gray70, WT_scaffold and elem C
color bluewhite, new_pocket and elem C
set_bond stick_radius, 0.25, mutable
set label_size, 20
select don, (elem N or elem O and (neighbor hydro))
select acc, (elem O or (elem N and not (neighbor hydro)))
color bluewhite, ligand_poses and elem C
color gray70, WT_ligand and elem C
dist HBA, (ligand_poses and acc), (designs and don), 3.2
dist HBD, (ligand_poses and don), (designs and acc), 3.2
delete don
delete acc
hide labels, HBA
hide labels, HBD
color lightorange, HBA
color lightorange, HBD
zoom new_pocket
hide (hydro)
deselect
