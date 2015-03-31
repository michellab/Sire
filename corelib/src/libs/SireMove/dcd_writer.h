#include "fastio.h"

#include <OpenMM.h>   // CONDITIONAL_INCLUDE

#include "SireMove/integratorworkspace.h"

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) fio_fwrite(((void *) buf), (size), 1, (fd))

using SireMove::AtomicVelocityWorkspace;

static int write_dcdheader(fio_fd fd, const char *remarks, int N, int ISTART, int NSAVC, double DELTA, int with_unitcell, int charmm);

static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, float *X, float *Y, float *Z, const double *unitcell, int charmm);

class DCD {

    public:
        DCD(const char * file_name, bool wrap_coord,int cutoff_type,std::vector<double> & center, int frequency, int atoms , double timestep,
        OpenMM::Context & context, SireMove::AtomicVelocityWorkspace & work_space):
        box_center(center), context_openmm(context), ws(work_space){
				filename=file_name;
				wrap=wrap_coord;
				flag_cutoff=cutoff_type;
				frequency_dcd=frequency;
				natoms = atoms;
				tstep = timestep;
				fio_open(file_name,FIO_WRITE, &fd);

        }

	~DCD(){};

	void write(int);
	void close(void);

	private:
		const char * filename;
		bool wrap;
		int flag_cutoff;
		OpenMM::Context & context_openmm;
		fio_fd fd;
		std::vector<double> box_center;
		AtomicVelocityWorkspace &ws;
	    int frequency_dcd;
	    int natoms;
	    double tstep;

};
