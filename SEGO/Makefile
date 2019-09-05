HEADER = __header__.h
TARGET = se

OBJS := $(patsubst %.cc,%.o,$(wildcard *.cc))
CC = g++
OPTION = -O3

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(OPTION) -lm 

%.o: %.cc $(HEADER)
	$(CC) $(CFLAGS) $(OPTION)  -c $<

clean:
	del -rf *.o
	del $(TARGET).exe



