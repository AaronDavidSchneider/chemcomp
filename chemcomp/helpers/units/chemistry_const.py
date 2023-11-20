#   Adundances in terms of H (Asplund et al. 2009) + variation from Chiara plot
OHabu0 = 4.9e-4  # O/H abundance, solar value
CHabu0 = 2.69e-4  # C/H abundance, solar value
SiHabu0 = 3.24e-5  # Si/H adundance, solar value
FeHabu0 = 3.16e-5  # Fe/H adundance, solar value
SHabu0 = 1.32e-5  # S/H adundance, solar value
MgHabu0 = 3.98e-5  # Mg/H adundance, solar value
HeHabu0 = 0.085  # mass ratio of H in protoplanetary disk
AlHabu0 = 2.82e-6
TiHabu0 = 8.91e-8
KHabu0 = 1.07e-7
NaHabu0 = 1.74e-6
NHabu0 = 6.76e-5
VHabu0 = 8.59e-9

# for twopoppy comparison:
# OHabu0 = 8.91e-8  # O/H abundance, solar value
# CHabu0 = 0.0  # C/H abundance, solar value
# SiHabu0 = 0.0  # Si/H adundance, solar value
# FeHabu0 = 0.0  # Fe/H adundance, solar value
# SHabu0 = 0.0  # S/H adundance, solar value
# MgHabu0 = 0.0  # Mg/H adundance, solar value
# HeHabu0 = 0.00  # mass ratio of H in protoplanetary disk
# AlHabu0 = 0.0
# TiHabu0 = 8.91e-8
# KHabu0 = 0.0
# NaHabu0 = 0.0
# NHabu0 = 0.0
# VHabu0 = 0.0


#   Condensation temperatures from Öberg et al (2011) if not stated otherwise
#   code breaks if hirachy of temperatures is changedperature for CO
COtemp = 20  # condensation temperature for CO
N2temp = 20
CH4temp = 30  # condensation temperature for CH4
CO2temp = 70  # condensation temperature for CO2
H2Otemp = 150  # condensation temperature for H2O
Fe3O4temp = 371  # Lodders (2003)
FeStemp = 704  # Lodders (2003)
Mg2SiO4temp = 1354  # Lodders (2003)
Fe2O3temp = 1357  # Lodders (2003) for pure Fe
MgSiO3temp = 1500  # condensation temperature for MgSiO3
NH3temp = 90
H2Stemp = 150
Ctemp = 631
NaAlSi3O8temp = 958
KAlSi3O8temp = 1006
TiOtemp = 2000
Al2O3temp = 1653
VOtemp = 1423

MassO = 16.0
MassH = 1.0
MassFe = 56.0
MassMg = 24.3
MassSi = 28.0
MassS = 32.0
MassHe = 4.0
MassTi = 47.867
MassAl = 27
MassK = 39.0983
MassNa = 23
MassN = 14
MassV = 50.9415
MassC = 12.0  # C in terms of H

#   Masses in terms of H (atomic unit)
MassCO = MassC + MassO  #28.0  # CO mass in terms of H
MassCH4 = MassC + 4 * MassH  # CH4 in terms of H
MassCO2 = MassC + 2 * MassO   # CO2 in terms of H
MassH2O = MassO + 2 * MassH  # H20 in terms of H

MassFe2O3 = 2 * MassFe + 3 * MassO
MassFe3O4 = 3 * MassFe + 4 * MassO
MassFeS = MassFe + MassS
MassMgSiO3 = MassMg + MassSi + 3 * MassO  # MgSiO3 in terms of H
MassMg2SiO4 = 2 * MassMg + MassSi + 4 * MassO  # MgSiO3 in terms of H
MassNH3 = 3 * MassH + MassN
MassN2 = 2 * MassN
MassH2S = 2 * MassH  + MassS
MassNaAlSi3O8 = MassNa + MassAl + 3 * MassSi + 8 * MassO
MassKAlSi3O8 = MassK + MassAl + 3 * MassSi + 8 * MassO
MassTiO = MassTi + MassO
MassAl2O3 = 2 * MassAl + 3 * MassO
MassVO = MassV + MassO
