# SALome UTilitairES.
#
# SALUTES is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SALUTES is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# I kindly ask you to keep this memoir of alterations tidy.
# 
# |---------------|-----------------------|------------------------------------|
# | DATE          | AUTHOR                | DESCRIPTION                        |
# |---------------|-----------------------|------------------------------------|
# | Oct, 14, 2023 | Helio C. Bortolon     | Original issue.                    |
# |               |                       |                                    |
# |---------------|-----------------------|------------------------------------|
#

# Setting up the environment:
from salome.geom import geomBuilder
geompy = geomBuilder.New()

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

import numpy as num

import salome_pluginsmanager

def MakeSubByProp(objname, sstype, ssrang, basename):
	
	# Selecting the geometrical object by its path in the study:
	objx = salome.myStudy.FindObjectByPath('/Geometry/' + objname).GetObject()
	
	# Selecting the sub-objects:
	lbase = geompy.ExtractShapes(objx, geompy.ShapeType[sstype], True)
	
	# Sweeping the sub-objects for the ones which fit the criteria:
	n = 1
	losub = []
	for xl in lbase:
		ssbp = geompy.BasicProperties(xl)
		if ((ssrang[0] < ssbp[0]) and (ssbp[0] < ssrang[1]) and \
			(ssrang[2] < ssbp[1]) and (ssbp[1] < ssrang[3]) and \
			(ssrang[4] < ssbp[2]) and (ssbp[2] < ssrang[5])):
			losub.append(xl)
			geompy.addToStudyInFather(objx, xl, basename + "{0:02d}".format(n))
			n = n + 1
	
	# Creating a group with the collected objects:
	gsobj = geompy.CreateGroup(objx, geompy.ShapeType[sstype])
	geompy.UnionList(gsobj, losub)
	
	# Inserting the group in the father object:
	geompy.addToStudyInFather(objx, gsobj, basename)
	return losub

salome_pluginsmanager.AddFunction('Salutes/MakeSubByProp',
                                    'Generates a sub-object by its properties.',
                                    MakeSubByProp)
                                    
# Separate specific subshapes based on an external criteria function.
def MakeSubByCrit(objname, sstype, criter, basename):
	
	# Selecting the geometrical target:
	objx = salome.myStudy.FindObjectByPath('/Geometry/' + objname).GetObject()
	
	# Selecting the sub-objects of the type of inerest:
	lbase = geompy.ExtractShapes(objx, geompy.ShapeType[sstype], True)
	
	# Sweeping the objects list:
	n = 1
	losub = []
	for xl in lbase:
		# Testing whether or not the sub-object obeys the external crieteria:
		ok = criter(xl)
		if ok == True:
			losub.append(xl)
			geompy.addToStudyInFather(objx, xl, basename + "{0:02d}".format(n))
			n = n + 1
	
	# Creating a group with the collected objects:
	gsobj = geompy.CreateGroup(objx, geompy.ShapeType[sstype])
	geompy.UnionList(gsobj, losub)
	
	# Inserting the group in the father object:
	geompy.addToStudyInFather(objx, gsobj, basename)
	return losub

salome_pluginsmanager.AddFunction('Salutes/MakeSubByCrit',
                                    'Generates a sub-object through an external criteria.',
                                    MakeSubByCrit)

# Make bolts from a list of surfaces representing nuts.
def MakeBoltsFromNuts(lonuts, length):
	'''
	MakeBoltsFromNuts:
	Input:
		lonuts: list of objects representing the nuts' footprints.
		length: how long will be the bolts.
	output:
		lobolts: list of objects representing the bolts.
	'''
	lobolts = []
	for xl in lonuts:
		nut = xl.GetName()
		# Reference ppsition and direction:
		cog = num.array(geompy.PointCoordinates(geompy.MakeCDG(xl)))
		nor = num.array(geompy.VectorCoordinates(geompy.GetNormal(xl)))
		nor = nor / num.linalg.norm(nor)
		# Vertices creation:
		top = geompy.MakeVertex(*cog)
		pta = cog - 0.45 * length * nor
		pta = geompy.MakeVertex(*pta)
		ptb = cog - 0.55 * length * nor
		ptb = geompy.MakeVertex(*ptb)
		bot = cog - length * nor
		bot = geompy.MakeVertex(*bot)
		# Edges creation:
		e1 = geompy.MakeEdge(top, pta)
		e2 = geompy.MakeEdge(pta, ptb)
		e3 = geompy.MakeEdge(ptb, bot)
		# Assembling the edges:
		bolt = geompy.MakeFuseList([e1, e2, e3])
		# Exploding in vertices and edges:
		lov = geompy.ExtractShapes(bolt, geompy.ShapeType['VERTEX'], True)
		loe = geompy.ExtractShapes(bolt, geompy.ShapeType['EDGE'], True)
		# Inserting in the study:
		geompy.addToStudy(bolt, 'b' + nut)
		geompy.addToStudyInFather(bolt, lov[0], 'va' + nut)
		geompy.addToStudyInFather(bolt, lov[1], 'pa' + nut)
		geompy.addToStudyInFather(bolt, lov[2], 'pb' + nut)
		geompy.addToStudyInFather(bolt, lov[3], 'vb' + nut)
		geompy.addToStudyInFather(bolt, loe[0], 'ea' + nut)
		geompy.addToStudyInFather(bolt, loe[1], 'ep' + nut)
		geompy.addToStudyInFather(bolt, loe[2], 'eb' + nut)
		lobolts.append(bolt)
	return lobolts

salome_pluginsmanager.AddFunction('Salutes/MakeBoltsFromNuts',
                                    'Generates bolts projecting from nuts imprints.',
                                    MakeBoltsFromNuts)

# Make bolts from two lists of surfaces by pairing them.
def MakeBoltsFromNutsPairs(lon1, lon2, basename1, basename2):
	'''
	MakeBoltsFromNutsPairs:
	Input:
		lon1, lon2: list of objects to search for the pairs.
		basename1, basname2: basic names to search for.
	output:
		lobolts: list of objects representing the bolts.
	'''
	
	# Sweeping the sub-objects for the targets to generate decision parameters:
	doa = {}
	dob = {}
	for ox in lon1:
		oxname = ox.GetName()
		cog = num.array(geompy.PointCoordinates(geompy.MakeCDG(ox)))
		nor = num.array(geompy.VectorCoordinates(geompy.GetNormal(ox)))
		nor = nor / num.linalg.norm(nor)
		doa[ox] = {'cog': cog, 'nor': nor, 'flag': True}
	
	for ox in lon2:
		oxname = ox.GetName()
		cog = num.array(geompy.PointCoordinates(geompy.MakeCDG(ox)))
		nor = num.array(geompy.VectorCoordinates(geompy.GetNormal(ox)))
		nor = nor / num.linalg.norm(nor)
		dob[ox] = {'cog': cog, 'nor': nor, 'flag': True}

	# Pairing the target surfaces and creating the relative bolts:
	lobolts = []
	
	for oa in doa.keys():
		for ob in dob.keys():
			if not (ob == oa):
				acog = doa[oa]['cog']
				anor = doa[oa]['nor']
				bcog = dob[ob]['cog']
				bnor = dob[ob]['nor']
				abve = (bcog - acog)
				abvu = abve / num.linalg.norm(abve)
				# If the faces are aligned and available:
				if (abs(abs(num.dot(anor, bnor)) - 1.) < 1e-8) and \
				   (abs(abs(num.dot(anor, abvu)) - 1.) < 1e-8) and \
				   (abs(abs(num.dot(bnor, abvu)) - 1.) < 1e-8) and \
					doa[oa]['flag'] and dob[ob]['flag']:
						# Flag them as treated and proceed the bolt construction:
						doa[oa]['flag'] = False
						dob[ob]['flag'] = False
						# The needed vertices:
						va = geompy.MakeVertex(*acog)
						vb = geompy.MakeVertex(*bcog)
						pa = acog + .45 * abve
						pb = acog + .55 * abve
						pa = geompy.MakeVertex(*pa)
						pb = geompy.MakeVertex(*pb)
						# The needed edges:
						ea = geompy.MakeEdge(va, pa)
						ep = geompy.MakeEdge(pa, pb)
						eb = geompy.MakeEdge(pb, vb)
						# Joining everything:
						bolt = geompy.MakeFuseList([ea, ep, eb])
						# Extracting sub-shapes:
						lov = geompy.ExtractShapes(bolt, geompy.ShapeType['VERTEX'], True)
						loe = geompy.ExtractShapes(bolt, geompy.ShapeType['EDGE'], True)
						# Inserting in the study:
						geompy.addToStudy(bolt, 'b' + oa.GetName() + ob.GetName())
						geompy.addToStudyInFather(bolt, lov[3], 'va' + oa.GetName())
						geompy.addToStudyInFather(bolt, lov[2], 'pa' + oa.GetName() + ob.GetName())
						geompy.addToStudyInFather(bolt, lov[1], 'pb' + oa.GetName() + ob.GetName())
						geompy.addToStudyInFather(bolt, lov[0], 'vb' + ob.GetName())
						geompy.addToStudyInFather(bolt, loe[2], 'ea' + oa.GetName() + ob.GetName())
						geompy.addToStudyInFather(bolt, loe[1], 'ep' + oa.GetName() + ob.GetName())
						geompy.addToStudyInFather(bolt, loe[0], 'eb' + oa.GetName() + ob.GetName())
						lobolts.append(bolt)
	return lobolts

salome_pluginsmanager.AddFunction('Salutes/MakeBoltsFromNutsPairs',
                                    'Generates bolts linking pair of nuts imprints.',
                                    MakeBoltsFromNutsPairs)

# Making circles concentric to existing circles.
def MakeConcCircle(objtarg, radius, basename):
	'''
	MakeConcCircle:
	Input:
		objtarg: base circle over which the new one will be made.
				Grab it with 'salome.myStudy.FindObjectByPath(<path>).GetObject().
		radius: the new circle radius.
		basename: basic name for the new circles.
	output:
		None. The created objects are inserted into the Geometry tree.
	'''
	
	# Choosing from an isolated part or a list of parts?
	kos = geompy.KindOfShape(objtarg)
	if kos[0].__str__() == 'CIRCLE':
		# Collectint information on the object:
		cog = geompy.MakeVertex(*kos[1:4])
		nor = geompy.MakeVectorDXDYDZ(*kos[4:7])
		rad = kos[7]
		# Making the concentric circle:
		circ = geompy.MakeCircle(cog, nor, radius)
		geompy.addToStudyInFather(objtarg, circ, basename)
	elif kos[0].__str__() == 'COMPOUND':
		# Open the list of objects:
		loso = geompy.ExtractShapes(objtarg, geompy.ShapeType['EDGE'], True)
		i = 1
		locirc = []
		# Sweeping he list:
		for ox in loso:
			# Collecting information on the object:
			kox = geompy.KindOfShape(ox)
			cog = geompy.MakeVertex(*kox[1:4])
			nor = geompy.MakeVectorDXDYDZ(*kox[4:7])
			rad = kox[7]
			# Making the concenctric circle:
			circ = geompy.MakeCircle(cog, nor, radius)
			locirc.append(circ)
			#geompy.addToStudy(circ, basename + str(i))
			i = i + 1
		# Creating a group to insert the set of circles:
		geompy.MakeFuseList(locirc, theName = basename)
		geompy.addToStudy(locirc, basename)

salome_pluginsmanager.AddFunction('Salutes/MakeConcCircle',
                                    'Generate a concentric circle.',
                                    MakeConcCircle)

# Meshing a group of bolts.
def MakeBoltsMeshes(lobolts):
	noseg = smesh.CreateHypothesis('NumberOfSegments')
	noseg.SetNumberOfSegments( 1 )
	reg1d = smesh.CreateHypothesis('Regular_1D')
	for bolt in lobolts:
		# Starting the mesh:
		boltmesh = smesh.Mesh(bolt)
		status = boltmesh.AddHypothesis(noseg)
		status = boltmesh.AddHypothesis(reg1d)
		isDone = boltmesh.Compute()
		# Getting the targets for the mesh groups from the geometry:
		bname = bolt.GetName()
		lovtx = geompy.ExtractShapes(bolt, geompy.ShapeType['VERTEX'], True)
		loedg = geompy.ExtractShapes(bolt, geompy.ShapeType['EDGE'], True)
		# Edge groups:
		boltall = boltmesh.GroupOnGeom(bolt, bname, SMESH.EDGE)
		bolt1 = boltmesh.GroupOnGeom(loedg[0], bname + '1', SMESH.EDGE)
		boltx = boltmesh.GroupOnGeom(loedg[1], bname + 'p', SMESH.EDGE)
		bolt2 = boltmesh.GroupOnGeom(loedg[2], bname + '2', SMESH.EDGE)
		# NOdes groups:
		boltall = boltmesh.GroupOnGeom(bolt, bname, SMESH.NODE)
		bv1 = boltmesh.GroupOnGeom(lovtx[0], bname + '_v1', SMESH.NODE)
		bv2 = boltmesh.GroupOnGeom(lovtx[3], bname + '_v2', SMESH.NODE)
		bp1 = boltmesh.GroupOnGeom(lovtx[1], bname + '_p1', SMESH.NODE)
		bp2 = boltmesh.GroupOnGeom(lovtx[2], bname + '_p2', SMESH.NODE)

salome_pluginsmanager.AddFunction('Salutes/MakeBoltsMeshes',
                                    'Make the meshes of bolts.',
                                    MakeBoltsMeshes)

def MakeJunctionK(center, length, angle, name):
	'''
	Input:
		center: list of coordinates.
		length: length of the arms.
		angle: angle of the K, in degrees.
		name: name of the junction.
	Output:
		None, the objects are inserted automatically in the tree.
	'''
	alpha = angle * math.pi / 180.
	
	p = geompy.MakeVertex(*center)
	pl1 = geompy.MakeVertexWithRef(p, length, 0, 0)
	pl2 = geompy.MakeVertexWithRef(p, length * math.cos(alpha), length * math.sin(alpha), 0)
	pl3 = geompy.MakeVertexWithRef(p, -length * math.cos(alpha), length * math.sin(alpha), 0)
	pl4 = geompy.MakeVertexWithRef(p, -length, 0, 0)

	l1 = geompy.MakeLineTwoPnt(p, pl1)
	l2 = geompy.MakeLineTwoPnt(p, pl2)
	l3 = geompy.MakeLineTwoPnt(p, pl3)
	l4 = geompy.MakeLineTwoPnt(p, pl4)
	
	kj = geompy.MakeCompound([l1, l2, l3, l4], theName = name)
	
	lol = geompy.ExtractShapes(kj, geompy.ShapeType['EDGE'], theName = 'L')

salome_pluginsmanager.AddFunction('Salutes/MakeJunctionK',
                                    'Make a floating K junction.',
                                    MakeJunctionK)
