
#include "sireglobal.h"

#ifdef Q_OS_WIN
  #include "process_windows.cpp"
#else
  #include "process_unix.cpp"
#endif

