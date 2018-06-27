
ATSYMB        = atsymb.f90
CONSTANTS     = constants.f90

BBFS          = bbfs.f90
DIFF_AND_DISS = diff_and_diss.f90
DIFF_BIM      = diff_bim.f90
DISS          = diss.f90
INSTANT_TEMP  = instant_temp.f90
INTERN2       = intern2.f90
INTERN        = intern.f90
KMC_BIM	      = kmc_bim.f90
KMC           = kmc.f90
KMCEND        = kmcEND.f90
NM            = nm.f90
RESCALE       = rescale.f90
RRKM          = omegas.f90 constants_rrkm.f90 rrkm.f
ROTATE        = rotate.f90
SPRINT        = sprint.f
SPRINT2       = sprint2.f
SYMM0         = symm0.f
SYMM          = symm.f
TERMO         = termo.f90
TST_BIM       = tst_bim.f90
TST           = tst.f90


bin_PROGRAMS=bbfs.exe diff_and_diss.exe diff_bim.exe \
	     diss.exe instant_temp.exe intern2.exe intern.exe \
	     kmc_bim.exe kmc.exe kmcEND.exe nm.exe rescale.exe \
	     rrkm.exe rotate.exe sprint.exe sprint2.exe symm0.exe \
	     symm.exe termo.exe tst_bim.exe tst.exe     

bbfs_exe_SOURCES         = $(BBFS)
diff_and_diss_exe_SOURCES= $(CONSTANTS) $(DIFF_AND_DISS)
diff_bim_exe_SOURCES     = $(DIFF_BIM) $(CONSTANTS)
diss_exe_SOURCES         = $(ATSYMB) $(DISS) 
instant_temp_exe_SOURCES = $(INSTANT_TEMP) $(ATSYMB) $(CONSTANTS)
intern2_exe_SOURCES      = $(INTERN2) $(ATSYMB)
intern_exe_SOURCES       = $(INTERN) $(ATSYMB)
kmc_bim_exe_SOURCES      = $(KMC_BIM) 
kmc_exe_SOURCES          = $(KMC)
kmcEND_exe_SOURCES       = $(KMCEND) 
nm_exe_SOURCES           = $(NM) $(ATSYMB) $(CONSTANTS)
rescale_exe_SOURCES      = $(RESCALE) $(ATSYMB) $(CONSTANTS)
rrkm_exe_SOURCES         = $(RRKM) $(MOD3) $(MOD4)
rotate_exe_SOURCES       = $(ROTATE) $(ATSYMB) $(CONSTANTS)
sprint_exe_SOURCES       = $(SPRINT) 
sprint2_exe_SOURCES      = $(SPRINT2)
symm0_exe_SOURCES        = $(SYMM0)
symm_exe_SOURCES         = $(SYMM)
termo_exe_SOURCES        = $(TERMO) $(ATSYMB) $(CONSTANTS)
tst_bim_exe_SOURCES      = $(TST_BIM) $(CONSTANTS)
tst_exe_SOURCES          = $(TST) $(CONSTANTS)

