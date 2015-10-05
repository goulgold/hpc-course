all:
	icpc trap.c -o trap -lpthread

mac:
	gcc trap.c -o trap -lpthread -std=c99
	gcc trap_busywait.c -o trap_busywait -lpthread -std=c99 -g
	gcc trap_mutex.c -o trap_mutex -lpthread -std=c99 -g

clean:
	rm -f trap
