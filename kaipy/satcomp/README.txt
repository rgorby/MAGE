Main functions and dictionary structures

scutils.py
	getCdasData:
		returns dictionary of satellite data, either:
			exactly as delivered by Cdaweb
			or, full days of all days needed, as delivered by Cdawev

scRCM.py
	getSpecieslambdata:
		returns:		
			result{
				'ilamc'  : lambda values (centers)
				'ilami'  : lambda interfacts (1 longer than ilamc)
				'lamscl' : scaling factor needed to calculate differential flux
			}

	getSCOmniDiffFlux:
		Get differential flux of given spacecraft dataset
		returns:
			ephdata : EMFISIS trajectory data
			dataset{
				'name' : spacecraft name
				'species'  : species name, 'ions' or 'electrons'
				'epoch' : list of epoch times
				'OmniDiffFlux' : omni-directional differential flux of species
				'energies' : energy channels. May be 1D (len(epoch)) or 2D (# channel bins, len(epoch))
			}
	getRCMTimes:
		Note: Needs both rcm.h5 and mhdrcm.h5 to get MJD times
		returns:
			rcmTimes{
				'Nt' : number of time steps
				'sIDs' : ordered step ids (int)
				'sIDstrs' : step ids to access rcm.h5 keys (str)
				'T' : time (from start of simulation)
				'MJD' : MJD's 
			}
	getRC_scTrack:
		Pull RCM data along a given spacecraft track
		sc = data from spacecraft track file (e.g. RBSPB.sc.h5)
		returns:
			result{
				'T' : sc['T']
				'MJD' : sc['MJDs']
				'MLAT' : sc['MLAT']
				'MLON' : sc['MLON']
				'vm' : rcm vm (= bVol**(-2/3))
				'xmin' : rcm xmin
				'ymin' : rcm ymin
				'zmin' : rcm zmin
				'eqmin' : sqrt(xmin**2 + ymin**2)
				'electrons' : speciesdata['electrons'] : {
								'energies' : energies in keV
								'eetas' : species etas
								'diffFlux' : omni-directional differential flux (bMirror adjusted)
				}
				'ions' : speciesdata['ions']
			}

	consolidateODFs:
		Consolidate RCM track and sc data to same energy grid
		returns:
			result{
				'energyGrid' : shared energy grid
				'sc' : {
					'name' : sc name
					'time' : sc epoch
					'diffFlux' : omni diff flux
				}
				'rcm' : {
					'time' : rcm MJDs
					'origEGrid' : original energy information
					'orig ODF' : original omni diff flux corresponding to 'origEGrid'
					'diffFlux' : new omni diff flux mapped to 'energyGrid'
				}
			}

	getIntentitiesVsL:
		Calculate rcm pressure on tkl (time, energy, L shell) grid
		returns:
			result{
				'T' : rcm times
				'MJD' : rcm MJDs
				'lambda' : lambda values
				'L_bins' : list of L shell bin values
				'energyGrid' : fixed energyGrid values
				'press_tkl' : rcm pressure on tkl grid
				'press_tl' : rcm pressure on tl grid (summed along energies)
			}

	getRCM_eqlatlon:
		Get (MHD)RCM variables in a nice package
		returns:
			result{
				'T' : rcm times
				'MJD' : rcm MJDs
				'MLAT' : rcm mlat values (1D)
				'MLON' : rcm mlon values (1D)
				'xmin' : rcm xmins
				'ymin' : rcm ymins
				'press' : rcm pressure
			}



