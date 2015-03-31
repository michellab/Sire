
#include <mpi.h>

#include <iostream>

int main(int argc, char **argv)
{
    MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);

    if (not MPI::Is_initialized())
    {
        std::cout << "MPI is not initialized!\n";
        return -1;
    }

    if (MPI::Query_thread() != MPI_THREAD_MULTIPLE)
    {
        std::cout << "MPI multi-thread support is not available!\n";
        return -1;
    }

    if (MPI::Is_finalized())
    {
        std::cout << "Why has MPI finalized now???\n";
        return -1;
    }

    MPI::Finalize();

    if (not MPI::Is_finalized())
    {
        std::cout << "Why has MPI not finalized???\n";
        return -1;
    }

    std::cout << "This implementation of MPI can support Sire.\n";

    return 0;
}
