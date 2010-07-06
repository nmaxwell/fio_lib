include makefile.option

LIB_SRC	= 	vecmatop.cpp	bfio.cpp bfio_eval.cpp


LIB_OBJ	 = 	$(LIB_SRC:.cpp=.o)	

TST_SRC = 	tt.cpp

DEP     = 	$(LIB_SRC:.cpp=.d) $(TST_SRC:.cpp=.d)

libbfio.a:	${LIB_OBJ}
	$(AR) $(ARFLAGS) $@ $(LIB_OBJ)
	$(RANLIB) $@

tt:	libbfio.a tt.o
	${CXX} -o tt tt.o libbfio.a ${LDFLAGS}

-include $(DEP)

#------------------------------------------------------
tilde:
	rm -f *~

clean:
	rm -rf *~ *.d *.o *.a

tags:
	etags *.hpp *.cpp

