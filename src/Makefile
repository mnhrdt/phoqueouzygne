# configuration
CFLAGS  = -march=native -O3
LDLIBS  = -lm

# default target
s1aio   :

# bureaucracy
clean   : ; $(RM) *.o s1aio

# automatic dependences (necessary for directly-included files .c files)
.dep.mk : ; $(CC) -MM *.c > $@
-include .dep.mk
