import os
import shutil
#TODO this can be added to the generates of the OneELectron.dat file or they could be read from there. 
def nice_par(p):
    return "e" if p==1 else 'o' 

nc_orb = 28
no_orb = 2
ntot = nc_orb+no_orb
max_l = 6
max_lc = 2
max_le = 7
itake = 1
specnot = ["s","p","d","f","g","h","i"]
Ji, pi = [ int(pp) for pp in input("Initial angular momentum and parity: ").split(" ")]
Jf, pf = [ int(pp) for pp in input("Final angular momentum and parity: ").split(" ")]
#Ji = 1
#pi = -1
#Jf = 2
#pf = 1

boundflag = int(input("Lu Fano plot (0) or photoionization (1)? : "))

if(boundflag==0):
    dire = 'J%i_BS_%s'%(Jf,nice_par(pf))
else:
    name = input("Photionization identifier: ")
    dire = 'J%i_PI_%s_%s'%(Jf,nice_par(pf),name)
os.makedirs(dire, exist_ok=True)
#Initial state
init_ap = '''
!c -- Note that init should be 1 for the initial state ALWAYS. 
!c      The value of itake is the lone state whose photoabsorption will
!c      be calculated later in jjstream.  However Neig specifies 
!c      that in BoundToBound.f, the first Neig initial states of this
!c      symmetry will have their bound-to-bound oscillator strengths calculated
!c     Note that Neig, specified in the first line of input(unit 5) above, 
!c     gives the number of initial states for which oscillator strengths
!c     are desired.
'''
init = open("%s/RaInitialStateJ%i.dat"%(dire,Ji), mode='w')
init.write("%i %i %i %i %i   INPUT PARAMETERS: parity= +/- 1,  J, init, itake, Neig\n"%(pi, Ji, 1, itake, 2))
init.write(init_ap)
init.close()

#final state
#Determine the channels
js  = []
for l in range(max_l+1):
    for jj in set([abs(l+1/2),abs(l-1/2)]): js.append("%s %.1f"%(specnot[l],jj))

code = ""
epilogue = "\n\n\n\n"
nchan = 0
for lc in range(max_l+1):
    jc = set([abs(lc+1/2), abs(lc-1/2)])
    for jjc in jc:
        ind = js.index("%s %.1f"%(specnot[lc], jjc))
        epilogue += "%s %i - %i \n"%(js[ind], ind*(nc_orb+no_orb)+1,(ind+1)*(nc_orb+no_orb))

for lc in range(max_lc+1):
    for le in range(max_le+1): 
        jc = set([abs(lc+1/2), abs(lc-1/2)])
        je = set([abs(le+1/2), abs(le-1/2)])
        for jjc in jc:
            for jje in je:
                jadd = [jj for jj in range(int(abs(jjc-jje)),int(jjc+jje+1))]
                if((Jf in jadd) and ((-1)**((lc+le)%2)==pf)):
                    nchan += 1
                    ic = ntot*js.index('%s %.1f'%(specnot[lc],jjc))+1
                    tie = js.index('%s %.1f'%(specnot[le],jje))
                    ie = [tie*(nc_orb+no_orb)+nc_orb+1+i for i in range(no_orb)]
                    code+="%i %i %s %s\n"%(ic, ie[0], js[(ic-1)//ntot], js[tie])
                    for i in ie[1:]:
                        code+="%i %i\n"%(ic, i)

final = open("%s/RaFinalStateJ%i.dat"%(dire,Jf), mode='w')
final.write("%i %i %i **** parity= +/- 1,  J, initial(no != 1)\n"%(pf, Jf, 99999))
final.write("%i %i         **** number of open orbitals, number of channels\n"%(no_orb*nchan, nchan))
final.write(code)
final.write(epilogue)
final.close()

# After putting the state files in the folder we can then copy the excutables there and generate the
# energy file.

if(boundflag==0):
    elo, ehi, ne, nnew = input("Energy range as E low, E high, ne, channels to artificially open: ").split(" ")
    elo = float(elo)
    ehi = float(ehi)
    ne = int(ne)
    nnew = int(nnew)
else:
    elo, ehi, ne = str(input("Energy range as E low, E high, ne: ")).split(" ")
    elo = float(elo)
    ehi = float(ehi)
    ne = int(ne)

with open("%s/Einputs.dat"%(dire), "w") as file:
    if(boundflag==0):
        file.write("%.8e   %.8e   %i   %i\n"%(elo, ehi, ne, nnew))
    else:
        file.write("%.8e   %.8e   %i\n"%(elo, ehi, ne))

excutables = ["Makefile","jjwaveSAVE.f", 'jjr12.f', 'jjrmat.f', 'jjr12b.f', 'jjrmatb.f', 'jjstreamTheThr2.f', 'jjstream.f', 'Ra.gksub.f', 'gensub.f', 'matsub.f', 'OneElectron.dat']

for ex in excutables:
    shutil.copyfile(ex, dire+"/"+ex)

print("Run make in the folder to generate and use the excutables")