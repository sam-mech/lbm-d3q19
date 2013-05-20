CFLAGS=-g -Wall -pg -O2
LDLIBS=-lm
BIN=lbm

run: $(BIN)
	./$(BIN)


profile: clean $(BIN) run
	gprof $(BIN) > $(BIN)-profile.txt
	less $(BIN)-profile.txt

clean:
	$(RM) $(BIN)
	$(RM) *.o
