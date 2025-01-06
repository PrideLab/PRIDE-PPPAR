CC = gcc
CFLAGS = -std=c99 -O0
LIBS = -lm
OBJ_PATH = .obj

mhm: $(OBJ_PATH)/model.o $(OBJ_PATH)/rdresoh.o $(OBJ_PATH)/rdresoi.o
	$(CC) $(CFLAGS) -o mhm $(OBJ_PATH)/model.o $(OBJ_PATH)/rdresoh.o $(OBJ_PATH)/rdresoi.o $(LIBS)

$(OBJ_PATH):
	mkdir -p $(OBJ_PATH)

$(OBJ_PATH)/model.o: model.c | $(OBJ_PATH)
	$(CC) $(CFLAGS) -c model.c -o $(OBJ_PATH)/model.o

$(OBJ_PATH)/rdresoh.o: rdresoh.c | $(OBJ_PATH)
	$(CC) $(CFLAGS) -c rdresoh.c -o $(OBJ_PATH)/rdresoh.o

$(OBJ_PATH)/rdresoi.o: rdresoi.c | $(OBJ_PATH)
	$(CC) $(CFLAGS) -c rdresoi.c -o $(OBJ_PATH)/rdresoi.o

.PHONY: clean
clean:
	rm -f $(OBJ_PATH)/*.o mhm
