#line 228 "/home/zhiyuan/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "mdmain.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "stdio.h"
#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "split.h"
#include "global.h"

/************************************************************************/

/* routine that implements the time-steps. Called by main routine and calls others */
double MDMAIN(long NSTEP, long NPRINT, long NSAVE, long NORD1, long ProcID)
{
    double XTT;
    long i;
    double POTA,POTR,POTRF;
    double XVIR,AVGT,TEN;
    double TTMV = 0.0, TKIN = 0.0, TVIR = 0.0;

    /*.......ESTIMATE ACCELERATION FROM F/M */
    INTRAF(&gl->VIR,ProcID);

    {
#line 43
	unsigned long	Error, Cycle;
#line 43
	long		Cancel, Temp;
#line 43

#line 43
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 43
	if (Error != 0) {
#line 43
		printf("Error while trying to get lock in barrier.\n");
#line 43
		exit(-1);
#line 43
	}
#line 43

#line 43
	Cycle = (gl->start).cycle;
#line 43
	if (++(gl->start).counter != (NumProcs)) {
#line 43
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 43
		while (Cycle == (gl->start).cycle) {
#line 43
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 43
			if (Error != 0) {
#line 43
				break;
#line 43
			}
#line 43
		}
#line 43
		pthread_setcancelstate(Cancel, &Temp);
#line 43
	} else {
#line 43
		(gl->start).cycle = !(gl->start).cycle;
#line 43
		(gl->start).counter = 0;
#line 43
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 43
	}
#line 43
	pthread_mutex_unlock(&(gl->start).mutex);
#line 43
};

    INTERF(ACC,&gl->VIR,ProcID);

    {
#line 47
	unsigned long	Error, Cycle;
#line 47
	long		Cancel, Temp;
#line 47

#line 47
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 47
	if (Error != 0) {
#line 47
		printf("Error while trying to get lock in barrier.\n");
#line 47
		exit(-1);
#line 47
	}
#line 47

#line 47
	Cycle = (gl->start).cycle;
#line 47
	if (++(gl->start).counter != (NumProcs)) {
#line 47
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 47
		while (Cycle == (gl->start).cycle) {
#line 47
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 47
			if (Error != 0) {
#line 47
				break;
#line 47
			}
#line 47
		}
#line 47
		pthread_setcancelstate(Cancel, &Temp);
#line 47
	} else {
#line 47
		(gl->start).cycle = !(gl->start).cycle;
#line 47
		(gl->start).counter = 0;
#line 47
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 47
	}
#line 47
	pthread_mutex_unlock(&(gl->start).mutex);
#line 47
};

    /* MOLECULAR DYNAMICS LOOP OVER ALL TIME-STEPS */

    for (i=1;i <= NSTEP; i++) {
        TTMV=TTMV+1.00;

        /* reset simulator stats at beginning of second time-step */

        /* POSSIBLE ENHANCEMENT:  Here's where one start measurements to avoid
           cold-start effects.  Recommended to do this at the beginning of the
           second timestep; i.e. if (i == 2).
           */

        /* initialize various shared sums */
        if (ProcID == 0) {
            long dir;
            if (i >= 2) {
                {
#line 65
	struct timeval	FullTime;
#line 65

#line 65
	gettimeofday(&FullTime, NULL);
#line 65
	(gl->trackstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 65
};
            }
            gl->VIR = 0.0;
            gl->POTA = 0.0;
            gl->POTR = 0.0;
            gl->POTRF = 0.0;
            for (dir = XDIR; dir <= ZDIR; dir++)
                gl->SUM[dir] = 0.0;
        }

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 76
	struct timeval	FullTime;
#line 76

#line 76
	gettimeofday(&FullTime, NULL);
#line 76
	(gl->intrastart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 76
};
        }

        {
#line 79
	unsigned long	Error, Cycle;
#line 79
	long		Cancel, Temp;
#line 79

#line 79
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 79
	if (Error != 0) {
#line 79
		printf("Error while trying to get lock in barrier.\n");
#line 79
		exit(-1);
#line 79
	}
#line 79

#line 79
	Cycle = (gl->start).cycle;
#line 79
	if (++(gl->start).counter != (NumProcs)) {
#line 79
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 79
		while (Cycle == (gl->start).cycle) {
#line 79
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 79
			if (Error != 0) {
#line 79
				break;
#line 79
			}
#line 79
		}
#line 79
		pthread_setcancelstate(Cancel, &Temp);
#line 79
	} else {
#line 79
		(gl->start).cycle = !(gl->start).cycle;
#line 79
		(gl->start).counter = 0;
#line 79
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 79
	}
#line 79
	pthread_mutex_unlock(&(gl->start).mutex);
#line 79
};
        PREDIC(TLC,NORD1,ProcID);
        INTRAF(&gl->VIR,ProcID);
        {
#line 82
	unsigned long	Error, Cycle;
#line 82
	long		Cancel, Temp;
#line 82

#line 82
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 82
	if (Error != 0) {
#line 82
		printf("Error while trying to get lock in barrier.\n");
#line 82
		exit(-1);
#line 82
	}
#line 82

#line 82
	Cycle = (gl->start).cycle;
#line 82
	if (++(gl->start).counter != (NumProcs)) {
#line 82
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 82
		while (Cycle == (gl->start).cycle) {
#line 82
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 82
			if (Error != 0) {
#line 82
				break;
#line 82
			}
#line 82
		}
#line 82
		pthread_setcancelstate(Cancel, &Temp);
#line 82
	} else {
#line 82
		(gl->start).cycle = !(gl->start).cycle;
#line 82
		(gl->start).counter = 0;
#line 82
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 82
	}
#line 82
	pthread_mutex_unlock(&(gl->start).mutex);
#line 82
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 85
	struct timeval	FullTime;
#line 85

#line 85
	gettimeofday(&FullTime, NULL);
#line 85
	(gl->intraend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 85
};
            gl->intratime += gl->intraend - gl->intrastart;
        }


        if ((ProcID == 0) && (i >= 2)) {
            {
#line 91
	struct timeval	FullTime;
#line 91

#line 91
	gettimeofday(&FullTime, NULL);
#line 91
	(gl->interstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 91
};
        }

        INTERF(FORCES,&gl->VIR,ProcID);

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 97
	struct timeval	FullTime;
#line 97

#line 97
	gettimeofday(&FullTime, NULL);
#line 97
	(gl->interend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 97
};
            gl->intertime += gl->interend - gl->interstart;
        }

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 102
	struct timeval	FullTime;
#line 102

#line 102
	gettimeofday(&FullTime, NULL);
#line 102
	(gl->intrastart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 102
};
        }

        CORREC(PCC,NORD1,ProcID);

        BNDRY(ProcID);

        KINETI(gl->SUM,HMAS,OMAS,ProcID);

        {
#line 111
	unsigned long	Error, Cycle;
#line 111
	long		Cancel, Temp;
#line 111

#line 111
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 111
	if (Error != 0) {
#line 111
		printf("Error while trying to get lock in barrier.\n");
#line 111
		exit(-1);
#line 111
	}
#line 111

#line 111
	Cycle = (gl->start).cycle;
#line 111
	if (++(gl->start).counter != (NumProcs)) {
#line 111
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 111
		while (Cycle == (gl->start).cycle) {
#line 111
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 111
			if (Error != 0) {
#line 111
				break;
#line 111
			}
#line 111
		}
#line 111
		pthread_setcancelstate(Cancel, &Temp);
#line 111
	} else {
#line 111
		(gl->start).cycle = !(gl->start).cycle;
#line 111
		(gl->start).counter = 0;
#line 111
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 111
	}
#line 111
	pthread_mutex_unlock(&(gl->start).mutex);
#line 111
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 114
	struct timeval	FullTime;
#line 114

#line 114
	gettimeofday(&FullTime, NULL);
#line 114
	(gl->intraend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 114
};
            gl->intratime += gl->intraend - gl->intrastart;
        }

        TKIN=TKIN+gl->SUM[0]+gl->SUM[1]+gl->SUM[2];
        TVIR=TVIR-gl->VIR;

        /*  check if potential energy is to be computed, and if
            printing and/or saving is to be done, this time step.
            Note that potential energy is computed once every NPRINT
            time-steps */

        if (((i % NPRINT) == 0) || ( (NSAVE > 0) && ((i % NSAVE) == 0))){

            if ((ProcID == 0) && (i >= 2)) {
                {
#line 129
	struct timeval	FullTime;
#line 129

#line 129
	gettimeofday(&FullTime, NULL);
#line 129
	(gl->interstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 129
};
            }

            /*  call potential energy computing routine */
            POTENG(&gl->POTA,&gl->POTR,&gl->POTRF,ProcID);
            {
#line 134
	unsigned long	Error, Cycle;
#line 134
	long		Cancel, Temp;
#line 134

#line 134
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 134
	if (Error != 0) {
#line 134
		printf("Error while trying to get lock in barrier.\n");
#line 134
		exit(-1);
#line 134
	}
#line 134

#line 134
	Cycle = (gl->start).cycle;
#line 134
	if (++(gl->start).counter != (NumProcs)) {
#line 134
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 134
		while (Cycle == (gl->start).cycle) {
#line 134
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 134
			if (Error != 0) {
#line 134
				break;
#line 134
			}
#line 134
		}
#line 134
		pthread_setcancelstate(Cancel, &Temp);
#line 134
	} else {
#line 134
		(gl->start).cycle = !(gl->start).cycle;
#line 134
		(gl->start).counter = 0;
#line 134
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 134
	}
#line 134
	pthread_mutex_unlock(&(gl->start).mutex);
#line 134
};

            if ((ProcID == 0) && (i >= 2)) {
                {
#line 137
	struct timeval	FullTime;
#line 137

#line 137
	gettimeofday(&FullTime, NULL);
#line 137
	(gl->interend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 137
};
                gl->intertime += gl->interend - gl->interstart;
            }

            POTA=gl->POTA*FPOT;
            POTR=gl->POTR*FPOT;
            POTRF=gl->POTRF*FPOT;

            /* compute some values to print */
            XVIR=TVIR*FPOT*0.50/TTMV;
            AVGT=TKIN*FKIN*TEMP*2.00/(3.00*TTMV);
            TEN=(gl->SUM[0]+gl->SUM[1]+gl->SUM[2])*FKIN;
            XTT=POTA+POTR+POTRF+TEN;

            if ((i % NPRINT) == 0 && ProcID == 0) {
                fprintf(six,"     %5ld %14.5lf %12.5lf %12.5lf  \
                %12.5lf\n %16.3lf %16.5lf %16.5lf\n",
                        i,TEN,POTA,POTR,POTRF,XTT,AVGT,XVIR);
            }
        }

        /* wait for everyone to finish time-step */
        {
#line 159
	unsigned long	Error, Cycle;
#line 159
	long		Cancel, Temp;
#line 159

#line 159
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 159
	if (Error != 0) {
#line 159
		printf("Error while trying to get lock in barrier.\n");
#line 159
		exit(-1);
#line 159
	}
#line 159

#line 159
	Cycle = (gl->start).cycle;
#line 159
	if (++(gl->start).counter != (NumProcs)) {
#line 159
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 159
		while (Cycle == (gl->start).cycle) {
#line 159
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 159
			if (Error != 0) {
#line 159
				break;
#line 159
			}
#line 159
		}
#line 159
		pthread_setcancelstate(Cancel, &Temp);
#line 159
	} else {
#line 159
		(gl->start).cycle = !(gl->start).cycle;
#line 159
		(gl->start).counter = 0;
#line 159
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 159
	}
#line 159
	pthread_mutex_unlock(&(gl->start).mutex);
#line 159
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 162
	struct timeval	FullTime;
#line 162

#line 162
	gettimeofday(&FullTime, NULL);
#line 162
	(gl->trackend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 162
};
            gl->tracktime += gl->trackend - gl->trackstart;
        }
    } /* for i */

    return(XTT);

} /* end of subroutine MDMAIN */
