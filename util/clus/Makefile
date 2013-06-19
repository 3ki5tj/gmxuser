CC=gcc -O0 -g3 -Wall -Wextra -pedantic
deps=zcom2.h

bin=test1 test2 testload
all: $(bin)

test1: test1.c $(deps)
	$(CC) -o $@ $< -lm

test2: test2.c $(deps)
	$(CC) -o $@ $< -lm

testload: testload.c $(deps)
	$(CC) -o $@ $< -lm

clean:
	rm -f $(bin) *~ *.o x.txt clus.txt

.PHONY: all clean

