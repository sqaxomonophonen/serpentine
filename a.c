#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "a.h"

void arghf(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	printf("\n");
	abort();
}

