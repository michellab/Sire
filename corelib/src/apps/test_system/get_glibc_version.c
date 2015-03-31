#include <stdio.h>
#include <gnu/libc-version.h>

int main (int argc, char **argv) 
{
    printf("%s\n", gnu_get_libc_version()); 
    return 0; 
}
