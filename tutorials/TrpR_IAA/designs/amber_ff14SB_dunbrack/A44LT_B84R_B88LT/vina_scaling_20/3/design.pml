bg_color white
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/3/receptor.pdb, design
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/designs/amber_ff14SB_dunbrack/A44LT_B84R_B88LT/vina_scaling_20/3/ligand.mol2, ligand
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/scaffold/amber_ff14SB/scaffold.pdb, WT_scaffold
hide lines
cartoon loop
color gray80, WT_scaffold and elem C
color gray90, design and elem C
set stick_radius, 0.15
show stick, ligand
hide cartoon ligand
select old_pocket, ((chain A and resi 44 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 84 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 88 and not (name C or name N or name O or name H or name HA))) and WT_scaffold
select new_pocket, ((chain A and resi 44 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 84 and not (name C or name N or name O or name H or name HA)) or (chain B and resi 88 and not (name C or name N or name O or name H or name HA))) and design
select mutable, ((chain A and resi 44) or (chain B and resi 88)) and design
label mutable and name ca and design, "%s - %s" %(resn,resi)
show stick, old_pocket
show lines, WT_scaffold and not old_pocket
hide lines, WT_scaffold and name c+o+n
show stick, new_pocket
set_bond stick_radius, 0.25, mutable
spectrum b, red_blue, new_pocket
color bluewhite, ligand and elem C
color splitpea, WT_scaffold and elem C
color gray80, WT_scaffold and not old_pocket and elem C
select don, (elem N or elem O and (neighbor hydro))
select acc, (elem O or (elem N and not (neighbor hydro)))
dist HBA, (ligand and acc), (design and don), 3.2
dist HBD, (ligand and don), (design and acc), 3.2
delete don
delete acc
hide labels, HBA
hide labels, HBD
color lightorange, HBA
color lightorange, HBD
set label_size, 20
load /agh/projects/jakob/PycharmProjects/PocketOptimizer2/docs/tutorials/TrpR_IAA/ligand/amber_ff14SB/ligand.mol2, WT_ligand
show stick, WT_ligand
color white, WT_ligand and elem C
zoom new_pocket
hide (hydro)
deselect
