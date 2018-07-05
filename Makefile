CFLAGS=-fPIC

libtetradec.so: sub_dsp.o fmat_tet.o sdec_tet.o sub_cd.o tetra_op.o sub_sc_d.o fexp_tet.o cdec_tet.o fbas_tet.o libtetradec.o
	gcc -shared *.o -o libtetradec.so

clean:
	rm *.o
	rm *.so
