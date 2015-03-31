#include "dcd_writer.h"
#include "SireMol/moleculegroup.h"
#include <QDebug>

using namespace SireMol;

void DCD::close(void){
    fio_fclose(fd);
}

enum {
    NOCUTOFF = 0,
    CUTOFFNONPERIODIC = 1,
    CUTOFFPERIODIC = 2
};


void DCD::write(int current_frame){


    static double COG[3] = {0.0,0.0,0.0};


    static double COT[3] = {0.0,0.0,0.0};//center of translation

    const int nmols = ws.nMolecules();

    int infoMask = 0;

    infoMask = OpenMM::State::Positions;

    infoMask = infoMask + OpenMM::State::Velocities; 

    OpenMM::State state_openmm = context_openmm.getState(infoMask);

    //OpenMM vector coordinate
    std::vector<OpenMM::Vec3> positions_openmm(natoms);

    //OpenMM vector momenta
    std::vector<OpenMM::Vec3> velocities_openmm(natoms);

    positions_openmm = state_openmm.getPositions();
    velocities_openmm = state_openmm.getVelocities();
    
    double *box_dims;
    box_dims = new double[6];

    float *X;
    float *Y;
    float *Z;

    X = new float[natoms];
    Y = new float[natoms];
    Z = new float[natoms];

    OpenMM::Vec3 a;
    OpenMM::Vec3 b;
    OpenMM::Vec3 c;

    int box = 0;

    double delta = tstep/0.0488821;


    //qDebug() << "N atoms = " << natoms;
    //qDebug() << "time step = " << tstep << "ps";


    if(flag_cutoff == CUTOFFPERIODIC){
        box=1;
    }

    if(current_frame == 0)
        write_dcdheader(fd, "Created by OpenMM", natoms,0,frequency_dcd, delta, box,1);

    if(wrap == true && current_frame == 0){
        for(unsigned int i=0;i<positions_openmm.size();i++){

            COG[0] = COG[0] + positions_openmm[i][0]*(OpenMM::AngstromsPerNm); //X Cennter of Geometry
            COG[1] = COG[1] + positions_openmm[i][1]*(OpenMM::AngstromsPerNm); //Y Cennter of Geometry
            COG[2] = COG[2] + positions_openmm[i][2]*(OpenMM::AngstromsPerNm); //Z Cennter of Geometry
        }

        COG[0] = COG[0]/natoms;
        COG[1] = COG[1]/natoms;
        COG[2] = COG[2]/natoms;
        
        COT[0]=box_center[0]-COG[0];
        COT[1]=box_center[1]-COG[1];
        COT[2]=box_center[2]-COG[2];
        
        /*cout << "\n\nBox center[0] = " << box_center[0] << " Box center[1] = " << box_center[1] << " Box center[2] = " << box_center[2] << "\n\n";

        cout << "COG[0] = " << COG[0] << " COG[1] = " << COG[1] << " COG[2] = " << COG[2] << "\n\n";

        cout << "COT[0] = " << COT[0] << " COT[1] = " << COT[1] << " COT[2] = " << COT[2] << "\n\n";*/

    }



    if(flag_cutoff == CUTOFFPERIODIC){

        state_openmm.getPeriodicBoxVectors(a,b,c);

        box_dims[0]=a[0] * OpenMM::AngstromsPerNm;
        box_dims[1]=0.0;
        box_dims[2]=b[1] * OpenMM::AngstromsPerNm;
        box_dims[3]=0.0;
        box_dims[4]=0.0;
        box_dims[5]=c[2] * OpenMM::AngstromsPerNm;

    }

    int num_atoms_till_l=0;
    
    


    for (int l=0; l < nmols; ++l){

        const int natoms_mol = ws.nAtoms(l);

        double molecule_COG[3] = {0.0,0.0,0.0};
        double T[3] = {0.0,0.0,0.0};//Translation vector

        if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){

            for (int m=0; m < natoms_mol; m++){

                int openmmindex = num_atoms_till_l+m;


                molecule_COG[0] = molecule_COG[0] + positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;;
                molecule_COG[1] = molecule_COG[1] + positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;;
                molecule_COG[2] = molecule_COG[2] + positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;;

            }

            molecule_COG[0] = molecule_COG[0]/natoms_mol;
            molecule_COG[1] = molecule_COG[1]/natoms_mol;
            molecule_COG[2] = molecule_COG[2]/natoms_mol;

            molecule_COG[0] = molecule_COG[0] + COT[0];
            molecule_COG[1] = molecule_COG[1] + COT[1];
            molecule_COG[2] = molecule_COG[2] + COT[2];
            

            T[0] = molecule_COG[0] - box_dims[0]*round(molecule_COG[0]/box_dims[0]);
            T[1] = molecule_COG[1] - box_dims[2]*round(molecule_COG[1]/box_dims[2]);
            T[2] = molecule_COG[2] - box_dims[5]*round(molecule_COG[2]/box_dims[5]);

            T[0] = T[0] - molecule_COG[0];
            T[1] = T[1] - molecule_COG[1];
            T[2] = T[2] - molecule_COG[2];

        }//end if wrap


        for (int m=0; m < natoms_mol; m++){

            int openmmindex = num_atoms_till_l+m;

            X[openmmindex] = positions_openmm[openmmindex][0]*OpenMM::AngstromsPerNm;
            Y[openmmindex] = positions_openmm[openmmindex][1]*OpenMM::AngstromsPerNm;
            Z[openmmindex] = positions_openmm[openmmindex][2]*OpenMM::AngstromsPerNm;

            if((wrap == true) && (flag_cutoff == CUTOFFPERIODIC)){


                X[openmmindex] = X[openmmindex] + COT[0];
                Y[openmmindex] = Y[openmmindex] + COT[1];
                Z[openmmindex] = Z[openmmindex] + COT[2];


                X[openmmindex] = X[openmmindex] + T[0];
                Y[openmmindex] = Y[openmmindex] + T[1];
                Z[openmmindex] = Z[openmmindex] + T[2];

           }

        }

        num_atoms_till_l=num_atoms_till_l+natoms_mol;
        
    }//end for molecules

    if(flag_cutoff == CUTOFFPERIODIC)
        write_dcdstep(fd, current_frame +1, current_frame + frequency_dcd, natoms, X, Y, Z,box_dims, 1);
    else
        write_dcdstep(fd, current_frame +1, current_frame + frequency_dcd, natoms, X, Y, Z,NULL, 1);

}


/*
 * Write a header for a new dcd file
 * Input: fd - file struct opened for binary writing
 *        remarks - string to be put in the remarks section of the header.  
 *                  The string will be truncated to 70 characters.
 *        natoms, istart, nsavc, delta - see comments in read_dcdheader
 * Output: 0 on success, negative error code on failure.
 * Side effects: Header information is written to the dcd file.
 */
static int write_dcdheader(fio_fd fd, const char *remarks, int N, 
                    int ISTART, int NSAVC, double DELTA, int with_unitcell,
                    int charmm) {
  int out_integer;
  float out_float;
  char title_string[200];
  time_t cur_time;
  struct tm *tmbuf;
  char time_str[81];

  out_integer = 84;
  WRITE(fd, (char *) & out_integer, sizeof(int));
  strcpy(title_string, "CORD");
  WRITE(fd, title_string, 4);
  fio_write_int32(fd, 0);      /* Number of frames in file, none written yet   */
  fio_write_int32(fd, ISTART); /* Starting timestep                            */
  fio_write_int32(fd, NSAVC);  /* Timesteps between frames written to the file */
  fio_write_int32(fd, 0);      /* Number of timesteps in simulation            */
  fio_write_int32(fd, 0);      /* NAMD writes NSTEP or ISTART - NSAVC here?    */
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    out_float = DELTA;
    WRITE(fd, (char *) &out_float, sizeof(float));
    if (with_unitcell) {
      fio_write_int32(fd, 1);
    } else {
      fio_write_int32(fd, 0);
    }
  } else {
    WRITE(fd, (char *) &DELTA, sizeof(double));
  }
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  fio_write_int32(fd, 0);
  if (charmm) {
    fio_write_int32(fd, 24); /* Pretend to be Charmm version 24 */
  } else {
    fio_write_int32(fd, 0);
  }
  fio_write_int32(fd, 84);
  fio_write_int32(fd, 164);
  fio_write_int32(fd, 2);

  strncpy(title_string, remarks, 80);
  title_string[79] = '\0';
  WRITE(fd, title_string, 80);

  cur_time=time(NULL);
  tmbuf=localtime(&cur_time);
  strftime(time_str, 80, "REMARKS Created %d %B, %Y at %R", tmbuf);
  WRITE(fd, time_str, 80);

  fio_write_int32(fd, 164);
  fio_write_int32(fd, 4);
  fio_write_int32(fd, N);
  fio_write_int32(fd, 4);

  return DCD_SUCCESS;
}


/* 
 * Write a timestep to a dcd file
 * Input: fd - a file struct for which a dcd header has already been written
 *       curframe: Count of frames written to this file, starting with 1.
 *       curstep: Count of timesteps elapsed = istart + curframe * nsavc.
 *        natoms - number of elements in x, y, z arrays
 *        x, y, z: pointers to atom coordinates
 * Output: 0 on success, negative error code on failure.
 * Side effects: coordinates are written to the dcd file.
 */
static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, 
                  float *X, float *Y, float *Z, 
                  const double *unitcell, int charmm) {
  int out_integer;
  //int rc;

  if (charmm) {
    /* write out optional unit cell */
    if (unitcell != NULL) {
      out_integer = 48; /* 48 bytes (6 doubles) */
      fio_write_int32(fd, out_integer);
      WRITE(fd, unitcell, out_integer);
      fio_write_int32(fd, out_integer);
    }
  }

  /* write out coordinates */
  out_integer = N*4; /* N*4 bytes per X/Y/Z array (N floats per array) */
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) X, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Y, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);
  fio_write_int32(fd, out_integer);
  if (fio_fwrite((void *) Z, out_integer, 1, fd) != 1) return DCD_BADWRITE;
  fio_write_int32(fd, out_integer);

  /* update the DCD header information */
  fio_fseek(fd, NFILE_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curframe);
  fio_fseek(fd, NSTEP_POS, FIO_SEEK_SET);
  fio_write_int32(fd, curstep);
  fio_fseek(fd, 0, FIO_SEEK_END);

  return DCD_SUCCESS;
}

