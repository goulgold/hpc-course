all:
	icpc trap.c -o trap -lpthread
	icpc trap_busywait.c -o trap_busywait -lpthread
	icpc trap_mutex.c -o trap_mutex -lpthread

mac:
	gcc trap.c -o trap -lpthread -std=gnu99
	gcc trap_busywait.c -o trap_busywait -lpthread -std=gnu99 -g
	gcc trap_mutex.c -o trap_mutex -lpthread -std=gnu99 -g

clean:
	rm -f trap
	rm -f trap_mutex
	rm -f trap_busywait
