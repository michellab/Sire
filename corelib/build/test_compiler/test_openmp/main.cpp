

//// Test to see if I can run an OpenMP loop in a background thread...

//// (this crashes on OS X)

#include <QThread>

#include <cmath>
#include <iostream>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
  #define omp_get_num_threads() 1
#endif

class OpenMPThread : public QThread
{
public:
    OpenMPThread() : QThread()
    {}
    
    ~OpenMPThread()
    {}
    
protected:
    void run()
    {
        const int array_size = 10;
    
        double a[array_size];
        
        for (int i=0; i<array_size; ++i)
        {
            a[i] = i;
        }

        std::cout << omp_get_num_threads() << std::endl;
        
        //#pragma omp for
        for (int i=0; i<array_size; ++i)
        {
            a[i] = std::sqrt(a[i]);
        }
    }
};


int main(int argc, const char **argv)
{
    OpenMPThread thread;
    
    thread.start();
    
    thread.wait();

    static const int array_size = 10;

    double a[array_size];

    #pragma omp for
    for (int i=0; i<array_size; ++i)
    {
        a[i] = i;
        std::cout << "Hello world (" << a[i] << ")\n";
    }
    
    return 0;
}
