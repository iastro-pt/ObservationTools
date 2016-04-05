# Ficheiro de entrada deve ser um rdb (sem header) contendo alpha,
# delta, name, magnitude, e.g.:
# 00:10:02.176    -82:13:26.57    star1   6.0
#

import string

def OB(comments, name, alpha, delta,name2,texp,texp2,filter,mv):

    texp_red= '%.1f' %  (texp)
    texp_blue = '%.1f' % (texp2)

    OB=("""IMPEX.VERSION "2.0"
STTimeIntervals         ""
calibrationReq          ""
type                    "O"
InstrumentComments      "%s"
userComments            "mv=%s"
userPriority            "1"
LineNumber              "0"
name                    "%s"


objectClass         "  Star              "
comments                ""
ra                      "%s"
dec                     "%s"
epoch                   "2000.0"
equinox                 "2000"
propDec                 "0.000000"
propRA                  "0.000000"
diffRA                  "0.000000"
diffDec                 "0.000000"
LineNumber              "0"
TARGET.NAME             "%s"


air_mass                      "2.0"
fractional_lunar_illumination "1.0"
sky_transparency              "Variable, thin cirrus"
moon_angular_distance         "30"
seeing                        "2.0"
CONSTRAINT.SET.NAME           "No Name"


IPVersion                    "254.09"
instrument                   "UVES"
longDescription              ""
LineNumber                   "0"
OBSERVATION.DESCRIPTION.NAME "No Name"


ACQUISITION.TEMPLATE.NAME "UVES_red_acq_imsl"
TEL.TARG.OFFSETALPHA "0.0"
TEL.TARG.OFFSETDELTA "0.0"
TEL.AG.GUIDESTAR     "CATALOGUE"
TEL.GS1.ALPHA        "0"
TEL.GS1.DELTA        "0"
INS.SLIT1.NAME       "SLIC#3"
INS.FILT1.NAME       "%s"
INS.DPOL.MODE        "OFF"
INS.OPTI1.NAME       "OUT"
INS.OPTI2.NAME       "OVRSIZ"




TEMPLATE.NAME "UVES_red_obs_exp"
DET2.READ.SPEED  "225kHz,1x1,low"
DET2.WIN1.UIT1   "%s"
SEQ.NEXPO        "1"
SEQ.SOURCE       "POINT"
SEQ.NOFF         "1"
TEL.TARG.OFFSETX "0"
TEL.TARG.OFFSETY "0"
INS.REDEXP.MODE  "580"
INS.SLIT3.WID    "0.3"
DPR.TYPE         "OBJECT"


TEMPLATE.NAME "UVES_red_obs_exp"
DET2.READ.SPEED  "225kHz,1x1,low"
DET2.WIN1.UIT1   "%s"
SEQ.NEXPO        "1"
SEQ.SOURCE       "POINT"
SEQ.NOFF         "1"
TEL.TARG.OFFSETX "0"
TEL.TARG.OFFSETY "0"
INS.REDEXP.MODE  "860"
INS.SLIT3.WID    "0.3"
DPR.TYPE         "OBJECT"
""" % (comments, mv,name, alpha, delta, name2, filter,texp_blue, texp_red))
    return OB


file=open('list.dat')

nobj=0
time=0
for line in file.readlines():

    if line[0]=="#": continue

    fields=string.split(line,'\t')
   
    name=string.replace(fields[2],' ','')
    name2=string.replace(fields[2],' ','')
    alpha=fields[0]
    delta=fields[1]
    priority='1'
    mv=float(fields[3])
    comments="S/N=200-250 (at maximum)"
    
    if mv<7.0 and mv>=4.5:
      filter = "ND1"
    elif mv>=2.0 and mv<4.5:
      filter = "ND2"
    elif mv<2.0:
      filter = "ND3"
    elif mv>=7.0:
      filter = "FREE"
      
    fc=string.replace(string.strip(name)+'.jpg',' ','')
      
    texp0=60.0 #in sec for a mv=7.0 ETC SN~250 in V-band
    mv0=7.0
    texp=(10**(-(mv0-mv)/2.5)*texp0)

   # Nota: multipliquei por um pequeno factor para diminuir o
   # tempo de exposicao para bater certo com o tempo total que o p2pp dava

    if texp<2.0:
      texp2=texp
    elif texp>=2.0:
      texp2=texp-(texp-2.0)*0.3

    print name, ' -- ', mv, ' -- ', float(texp), ' -- ', float(texp2)
 
    nobj=nobj+1
    time=time+texp2
  
    ob=open(string.replace(string.strip(name)+'.obx',' ',''),'w')
    ob.write(OB(comments, name, alpha, delta,name2,texp2,texp2,filter,mv))
    ob.close()

print 'Total time:', (time+nobj*600)/3600.   


