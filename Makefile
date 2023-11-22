
all: viterbi

clean:
	@rm -f viterbi

viterbi: viterbi.c
	${CC} -Wall -g -O0 -Wno-unused-function -DDEBUG=1 $^ -o $@

.PHONE=all clean
