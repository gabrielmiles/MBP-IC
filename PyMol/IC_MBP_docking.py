"""
Filename: OSM_MBP_docking.py
Date: 7/13/22
Author: Gabe Miles
Description: PyMol script for visualizing MBP/IC complexes.
Allows user to recreate positioning of complex in figures,
along with electrostatic and hydrophobicity maps. Includes
ray setting used for ray trace images.
"""

from pymol import cmd, stored

def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
	"""
	interfaceResidues -- finds 'interface' residues between two chains in a complex.
	
	PARAMS
		cmpx
			The complex containing cA and cB
		
		cA
			The first chain in which we search for residues at an interface
			with cB
		
		cB
			The second chain in which we search for residues at an interface
			with cA
		
		cutoff
			The difference in area OVER which residues are considered
			interface residues.  Residues whose dASA from the complex to
			a single chain is greater than this cutoff are kept.  Zero
			keeps all residues.
			
		selName
			The name of the selection to return.
			
	RETURNS
		* A selection of interface residues is created and named
			depending on what you passed into selName
		* An array of values is returned where each value is:
			( modelName, residueNumber, dASA )
			
	NOTES
		If you have two chains that are not from the same PDB that you want
		to complex together, use the create command like:
			create myComplex, pdb1WithChainA or pdb2withChainX
		then pass myComplex to this script like:
			interfaceResidues myComlpex, c. A, c. X
			
		This script calculates the area of the complex as a whole.  Then,
		it separates the two chains that you pass in through the arguments
		cA and cB, alone.  Once it has this, it calculates the difference
		and any residues ABOVE the cutoff are called interface residues.
			
	AUTHOR:
		Jason Vertrees, 2009.		
	"""
	# Save user's settings, before setting dot_solvent
	oldDS = cmd.get("dot_solvent")
	cmd.set("dot_solvent", 1)
	
	# set some string names for temporary objects/selections
	tempC, selName1 = "tempComplex", selName+"1"
	chA, chB = "chA", "chB"
	
	# operate on a new object & turn off the original
	cmd.create(tempC, cmpx)
	cmd.disable(cmpx)
	
	# remove cruft and inrrelevant chains
	cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
	
	# get the area of the complete complex
	cmd.get_area(tempC, load_b=1)
	# copy the areas from the loaded b to the q, field.
	cmd.alter(tempC, 'q=b')
	
	# extract the two chains and calc. the new area
	# note: the q fields are copied to the new objects
	# chA and chB
	cmd.extract(chA, tempC + " and (" + cA + ")")
	cmd.extract(chB, tempC + " and (" + cB + ")")
	cmd.get_area(chA, load_b=1)
	cmd.get_area(chB, load_b=1)
	
	# update the chain-only objects w/the difference
	cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
	
	# The calculations are done.  Now, all we need to
	# do is to determine which residues are over the cutoff
	# and save them.
	stored.r, rVal, seen = [], [], []
	cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

	cmd.enable(cmpx)
	cmd.select(selName1, 'none')
	for (model,resi,diff) in stored.r:
		key=resi+"-"+model
		if abs(diff)>=float(cutoff):
			if key in seen: continue
			else: seen.append(key)
			rVal.append( (model,resi,diff) )
			# expand the selection here; I chose to iterate over stored.r instead of
			# creating one large selection b/c if there are too many residues PyMOL
			# might crash on a very large selection.  This is pretty much guaranteed
			# not to kill PyMOL; but, it might take a little longer to run.
			cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

	# this is how you transfer a selection to another object.
	cmd.select(selName, cmpx + " in " + selName1)
	# clean up after ourselves
	cmd.delete(selName1)
	cmd.delete(chA)
	cmd.delete(chB)
	cmd.delete(tempC)
	# show the selection
	cmd.enable(selName)
	
	# reset users settings
	cmd.set("dot_solvent", oldDS)
	
	return rVal

cmd.extend("interfaceResidues", interfaceResidues)

def initial_orientation(path, config):
    """
    initial_orientation - positions complex in initial view

	PARAMS
		path - filename path for complex .pdb

		config - 1 or 2, based on docking configuration

	RETURNS
		Colored and orientated docking confirmation

	NOTES
		Use the docking files found in repository file 'original-models.'
		Be sure to use config option '1' for dock1, and '2' for dock2.
    """
    cmd.reinitialize()
    cmd.load(path)
    pathMtrx = path.split('/')
    filename = pathMtrx[len(pathMtrx) - 1]
    model = filename.replace(".pdb", "")
    cmd.set_name(model, "OM")
    cmd.show("surface")
    cmd.set("surface_color", "skyblue", "chain A")
    cmd.set("surface_color", "orange", "chain B")
    if config == '1':
        cmd.set_view((\
            0.547810793,    0.809414983,    0.211547419,\
            -0.178774685,    0.360282600,   -0.915552616,\
            -0.817279518,    0.463730723,    0.342069268,\
            0.000000000,    0.000000000, -257.350799561,\
            28.792110443,   50.535285950,   52.779388428,\
            202.897338867,  311.804260254,  -20.000000000 ))
    elif config == '2':
        cmd.set_view((\
            -0.739635348,   -0.367862731,   -0.563574135,\
            0.143549025,   -0.904357791,    0.401908547,\
            -0.657520115,    0.216365919,    0.721703112,\
            0.000000000,    0.000000000, -299.573425293,\
            67.182876587,   45.570850372,   51.338077545,\
            -32216.630859375, 32815.777343750,  -20.000000000 ))
cmd.extend('initial_orientation',initial_orientation)

def reset_view(config):
	"""
	reset_view - resets view to original loaded confirmation view

	PARAM
		config - 1 or 2, based on docking confirmation

	RETURNS
		original docking confirmation

	NOTES
		Ensure that correct confirmation is selected
	"""
	if config == '1':
		cmd.set_view((\
			0.547810793,    0.809414983,    0.211547419,\
			-0.178774685,    0.360282600,   -0.915552616,\
			-0.817279518,    0.463730723,    0.342069268,\
			0.000000000,    0.000000000, -257.350799561,\
			28.792110443,   50.535285950,   52.779388428,\
			202.897338867,  311.804260254,  -20.000000000 ))
	elif config == '2':
		cmd.set_view((\
			-0.739635348,   -0.367862731,   -0.563574135,\
			0.143549025,   -0.904357791,    0.401908547,\
			-0.657520115,    0.216365919,    0.721703112,\
			0.000000000,    0.000000000, -299.573425293,\
			67.182876587,   45.570850372,   51.338077545,\
			-32216.630859375, 32815.777343750,  -20.000000000 ))
cmd.extend('reset_view', reset_view)

def split_chains():
    """
    split_chains - from original orientation, splits MBP and IC, showing interacting surface

	PARAMS
		none
	
	RETURNS
		View of interacting surface
	
	NOTES
		If interacting residues are desired, run interfaceResidues before running this command.
		Be sure to save view before running APBS, plugin will shift view slightly after running.
    """
    cmd.select("MBP", "chain B")
    cmd.rotate("y", -90, "(MBP)")
    cmd.translate([-90, 0, 20], "(MBP)")

    cmd.select("OSM", "chain A")
    cmd.rotate("y", 90, "(OSM)")
    cmd.translate([0, 0, 20], "(OSM)")
    cmd.center()
    cmd.move("x", -10)
cmd.extend("split_chains", split_chains)

def ray_settings():
    """
    ray_settings - sets all necessary settings to recreate ray traces

	PARAMS
		none
	
	RETURNS
		PyMol set up for ray trace

    """
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.unset("depth_cue")
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("hash_max", 300)
    cmd.set("ambient", .4)
cmd.extend("ray_settings", ray_settings)

def hydrophobicity(selection='all'):
	"""
    PyMOL command to color protein molecules according to the Eisenberg hydrophobicity scale
	
	Source: http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html
	Amino acid scale: Normalized consensus hydrophobicity scale
	Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
	Reference: J. Mol. Biol. 179:125-142 (1984)

	Edited to show scale with purple by Gabe Miles
	
	Amino acid scale values:
	
	Ala:  0.620
	Arg: -2.530
	Asn: -0.780
	Asp: -0.900
	Cys:  0.290
	Gln: -0.850
	Glu: -0.740
	Gly:  0.480
	His: -0.400
	Ile:  1.380
	Leu:  1.060
	Lys: -1.500
	Met:  0.640
	Phe:  1.190
	Pro:  0.120
	Ser: -0.180
	Thr: -0.050
	Trp:  0.810
	Tyr:  0.260
	Val:  1.080
	
	Usage:
	hydrophobicity (selection)

	"""
	s = str(selection)
	print(s)
	cmd.set_color("color_ile3",[71, 67, 243])
	cmd.set_color("color_phe3",[88, 76, 244])
	cmd.set_color("color_val3",[103, 86, 245])
	cmd.set_color("color_leu3",[116, 95, 246])
	cmd.set_color("color_trp3",[128, 105, 247])
	cmd.set_color("color_met3",[139, 115, 248])
	cmd.set_color("color_ala3",[149, 124, 249])
	cmd.set_color("color_gly3",[159, 134, 250])
	cmd.set_color("color_cys3",[168, 144, 250])
	cmd.set_color("color_tyr3",[177, 153, 251])
	cmd.set_color("color_pro3",[186, 163, 252])
	cmd.set_color("color_thr3",[194, 173, 252])
	cmd.set_color("color_ser3",[202, 183, 253])
	cmd.set_color("color_his3",[210, 193, 253])
	cmd.set_color("color_glu3",[218, 203, 254])
	cmd.set_color("color_asn3",[226, 214, 254])
	cmd.set_color("color_gln3",[233, 224, 254])
	cmd.set_color("color_asp3",[241, 234, 255])
	cmd.set_color("color_lys3",[248, 245, 255])
	cmd.set_color("color_arg3",[255, 255, 255])
	cmd.color("color_ile3","("+s+" and resn ile)")
	cmd.color("color_phe3","("+s+" and resn phe)")
	cmd.color("color_val3","("+s+" and resn val)")
	cmd.color("color_leu3","("+s+" and resn leu)")
	cmd.color("color_trp3","("+s+" and resn trp)")
	cmd.color("color_met3","("+s+" and resn met)")
	cmd.color("color_ala3","("+s+" and resn ala)")
	cmd.color("color_gly3","("+s+" and resn gly)")
	cmd.color("color_cys3","("+s+" and resn cys)")
	cmd.color("color_tyr3","("+s+" and resn tyr)")
	cmd.color("color_pro3","("+s+" and resn pro)")
	cmd.color("color_thr3","("+s+" and resn thr)")
	cmd.color("color_ser3","("+s+" and resn ser)")
	cmd.color("color_his3","("+s+" and resn his)")
	cmd.color("color_glu3","("+s+" and resn glu)")
	cmd.color("color_asn3","("+s+" and resn asn)")
	cmd.color("color_gln3","("+s+" and resn gln)")
	cmd.color("color_asp3","("+s+" and resn asp)")
	cmd.color("color_lys3","("+s+" and resn lys)")
	cmd.color("color_arg3","("+s+" and resn arg)")
cmd.extend('hydrophobicity', hydrophobicity)


