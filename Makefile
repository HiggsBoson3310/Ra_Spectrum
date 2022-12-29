
SUBS := Ra.gksub.f gensub.f matsub.f

OPT_FLAG:= -std=legacy -O3
FC:=gfortran

jjwaveSAVE.x: jjwaveSAVE.f jjr12.x jjr12b.x jjrmat.x jjrmatb.x jjstreamTheThr2.x jjstream.x
	$(FC) jjwaveSAVE.f $(SUBS) -o jjwaveSAVE.x $(OPT_FLAG)

jjr12.x:jjr12.f
	$(FC) jjr12.f $(SUBS) -o jjr12.x $(OPT_FLAG)

jjr12b.x:jjr12b.f
	$(FC) jjr12b.f $(SUBS) -o jjr12b.x $(OPT_FLAG)

jjrmat.x: jjrmat.f 
	$(FC) jjrmat.f $(SUBS) -o jjrmat.x $(OPT_FLAG)

jjrmatb.x:jjrmatb.f
	$(FC) jjrmatb.f $(SUBS) -o jjrmatb.x $(OPT_FLAG)

jjstreamTheThr2.x:jjstreamTheThr2.f
	$(FC) jjstreamTheThr2.f $(SUBS) -o jjstreamTheThr2.x $(OPT_FLAG)

jjstream.x:  jjstream.f
	$(FC) jjstream.f $(SUBS) -o jjstream.x $(OPT_FLAG)






	
