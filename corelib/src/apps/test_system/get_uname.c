
#include <sys/utsname.h>

#include <stdio.h>

int main(int argc, char **argv)
{
    struct utsname unameData;

    if (uname(&unameData) == 0)
    {
        printf("uname.sysname = %s\n", unameData.sysname);
        printf("uname.nodename = %s\n", unameData.nodename);
        printf("uname.release = %s\n", unameData.release);
        printf("uname.version = %s\n", unameData.version);
        printf("uname.machine = %s\n", unameData.machine);
    }
    else
    {
        printf("FAILED - uname not available\n");
        return -1;
    }

    return 0;
}
