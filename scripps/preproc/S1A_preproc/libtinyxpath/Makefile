# Makefile for lib_src 
LIB = libtinyxpath.a
INC = -I../include -I$(NETCDFHOME)/include
CXX = g++ 
CXXFLAGS = -Bstatic $(INC) $(M64) -w

SRCS =	action_store.cpp \
	lex_util.cpp \
	node_set.cpp \
	tinystr.cpp \
	tinyxml.cpp \
	tinyxmlerror.cpp \
	tinyxmlparser.cpp \
	tokenlist.cpp \
	xml_util.cpp \
	xpath_expression.cpp \
	xpath_processor.cpp \
	xpath_stack.cpp \
	xpath_static.cpp \
	xpath_stream.cpp \
	xpath_syntax.cpp \

OBJS= $(SRCS:.cpp=.o)

$(LIB) : $(OBJS)
	ar rv $(LIB) $?
	ranlib  $(LIB)

clean:
	\rm -f *.o 
