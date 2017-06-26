/* PDN solver 
 * based on superLU 
 * zhiyuan yang 
 */
#include "slu_mt_ddefs.h"
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include "PDNConfig.h"



double **calculate_PDN_Midx(int **C4map, int **TSVmap, int Wnum, int Lnum, int layer_num, int *MidxSize)
{
	int var_num = Wnum * Lnum * layer_num + 2; 
	int i, j, l, k;

	/* now I calculate the number of non-zero entries and build Mid
 	*/ 
  	int count = 0; // count the number of non-zeros 
  	for (l = 0; l < layer_num; l ++){
        for (i = 0; i < Wnum; i ++){
            for (j = 0; j < Lnum; j ++){
                if (j+1 < Lnum){
                    count ++; 
                }
                if (j > 0){
                    count ++; 
                }
                if (i+1 < Wnum){
                    count ++; 
                }
                if (i > 0){
                    count ++;
                }

                if (l > 0){
                    if (TSVmap[i][j] == 1){
                        count ++; 
                    }
                    if (l+1 < layer_num && TSVmap[i][j] == 1){
                        count ++;
                    }
                    count ++;
                }
                else{
                    if (TSVmap[i][j] == 1){
                        count ++;
                    }
                    if (C4map[i][j] == 1){
                        count ++;

                        count ++; 
                    }
                    else{
                        count ++;
                    }
                }
            }
        }
    }
  	for (i = 0; i < Wnum; i++)
  		for (j = 0; j < Lnum; j++)
  			if (C4map[i][j] == 1)
  				count ++;
  	count = count + 4; 

  	double **Midx; 
  	// allocate space for Midx 
  	if ( !(Midx = (double **)malloc((count) * sizeof(double *))) ) printf("Malloc fails for Midx[].\n");
  	for (i = 0; i < count; i++)
    	if ( !(Midx[i] = (double *)malloc((3) * sizeof(double))) ) printf("Malloc fails for Midx[%d][].\n", i);

    // initialize the value
    for (i = 0; i < count; i ++)
    	for (j = 0; j < 3; j ++)
    		Midx[i][j] = 0;



    /// debug print the TSV map 
    /*printf("------- TSV map -------\n");
    for (i = 0; i < Wnum; i ++){
        for (j = 0; j < Lnum; j ++){
            printf("%d\t", TSVmap[i][j]);
        }
        printf("\n");
    }
    printf("-----------------------\n");
*/
    // filling the Midx array 
    int lat_con, vert_con, c4count; 
    int idx = 0; 
    for (l = 0; l < layer_num; l ++){
    	for (i = 0; i < Wnum; i ++){
    		for (j = 0; j < Lnum; j ++){
    			lat_con = 0; vert_con = 0; 
    			if (j+1 < Lnum){
    				Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j; 
    				Midx[idx][1] = l * Wnum*Lnum + i * Lnum + j + 1; 
    				Midx[idx][2] = -1 / R_GRID; 
    				lat_con ++;
    				idx ++; 
    			}
    			if (j > 0){
    				Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    				Midx[idx][1] = l * Wnum*Lnum + i * Lnum + j - 1; 
    				Midx[idx][2] = -1 / R_GRID; 
    				lat_con ++; 
    				idx ++; 
    			}
    			if (i+1 < Wnum){
    				Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    				Midx[idx][1] = l * Wnum*Lnum + (i+1) * Lnum + j;
    				Midx[idx][2] = -1 / R_GRID; 
    				lat_con ++;
    				idx ++; 
    			}
    			if (i > 0){
    				Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    				Midx[idx][1] = l * Wnum*Lnum + (i-1) * Lnum + j; 
    				Midx[idx][2] = -1 / R_GRID; 
    				lat_con ++; 
    				idx ++;
    			}

    			if (l > 0){
    				if (TSVmap[i][j] == 1){
    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][1] = (l-1) * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][2] = -1 / R_PDNTSV; 
    					vert_con ++;
    					idx ++; 
    				}
    				if (l+1 < layer_num && TSVmap[i][j] == 1){
    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j; 
    					Midx[idx][1] = (l+1) * Wnum*Lnum + i * Lnum + j; 
    					Midx[idx][2] = -1 / R_PDNTSV; 
    					vert_con ++; 
    					idx ++;
    				}
    				Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j; 
    				Midx[idx][1] = l * Wnum*Lnum + i * Lnum + j; 
    				Midx[idx][2] = 1 / R_GRID * lat_con + 1 / R_PDNTSV * vert_con;  
    				idx ++;
    			}
    			else{
    				if (TSVmap[i][j] == 1){
    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][1] = (l+1) * Wnum*Lnum + i * Lnum + j; 
    					Midx[idx][2] = -1 / R_PDNTSV;
    					vert_con ++;
    					idx ++;
    				}
    				if (C4map[i][j] == 1){
    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][1] = l * Wnum*Lnum + i * Lnum + j; 
    					Midx[idx][2] = 1 / R_GRID * lat_con + 1/ R_PDNTSV * vert_con + 1 / R_BUMP; 
    					idx ++;

    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][1] = layer_num * Wnum*Lnum; 
    					Midx[idx][2] = -1 / R_BUMP;
    					idx ++; 
    				}
    				else{
    					Midx[idx][0] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][1] = l * Wnum*Lnum + i * Lnum + j;
    					Midx[idx][2] = 1 / R_GRID * lat_con + 1 / R_PDNTSV * vert_con; 
    					idx ++;
    				}
    			}
    		}
    	}
    }
    /// for layer = 0 ///
    c4count = 0;
    for (i = 0; i < Wnum; i ++){
    	for (j = 0; j < Lnum; j ++){
    		if (C4map[i][j] == 1){
    			c4count ++; 
    			Midx[idx][0] = layer_num * Wnum*Lnum; 
    			Midx[idx][1] = i * Lnum + j;
    			Midx[idx][2] = -1 / R_BUMP;
    			idx ++;
    		}
    	}
    }
    /// for the off-chip resistance ///
    Midx[idx][0] = layer_num * Wnum*Lnum; 
    Midx[idx][1] = layer_num * Wnum*Lnum + 1;
    Midx[idx][2] = -1 / R_PKG;
    idx ++;

    Midx[idx][0] = layer_num * Wnum*Lnum; 
    Midx[idx][1] = layer_num * Wnum*Lnum;
    Midx[idx][2] = 1 / R_PKG + 1 / R_BUMP * c4count; 
    idx ++;

    Midx[idx][0] = layer_num * Wnum*Lnum + 1; 
    Midx[idx][1] = layer_num * Wnum*Lnum; 
    Midx[idx][2] = -1 / R_PKG;
    idx ++;

    Midx[idx][0] = layer_num * Wnum*Lnum + 1; 
    Midx[idx][1] = layer_num * Wnum*Lnum + 1; 
    Midx[idx][2] = 1 / R_PKG + 1 / R_PCB; 

    //////////////// for debug ////////////////
    if (idx != count-1)
    	printf("ERROR: count = %d;\t idx = %d\n", count, idx);
    else
        printf("Dimension Matches: count = %d\n", count);

    /////// sort the Midx ////////
    int idx_e, idx_s = 0; 
    int ind = 0;
    double row_t, col_t, val_t;
    for (i = 0; i < count; i ++){
    	if (Midx[i][0] > ind){
    		idx_e = i-1; 
    		ind = Midx[i][0]; 

    		for (l = idx_s; l < idx_e; l ++){
    			for (k = idx_s; k < idx_e; k ++){
    				if (Midx[k][1] > Midx[k+1][1]){
    					row_t = Midx[k][0]; col_t = Midx[k][1]; val_t = Midx[k][2];
                    	Midx[k][0] = Midx[k+1][0]; Midx[k][1] = Midx[k+1][1]; Midx[k][2] = Midx[k+1][2];
                    	Midx[k+1][0] = row_t; Midx[k+1][1] = col_t; Midx[k+1][2] = val_t; 
    				}
    			}
    		}
            idx_s = i;
    	}

    }

    //// print out the results 
    /*for (i = 0; i < count; i ++)
        printf("%.1f\t%.1f\t%.6f\n", Midx[i][0], Midx[i][1], Midx[i][2]);
    */

    *MidxSize = count; 
    return Midx;

}




double ***calculate_steady_PDN(double ***powerM, double **Midx, int count, int Wnum, int Lnum, int layer_num, double svdd)
{
    /* convert the values to the SuperMatrix format 
 */ 

	int i, j, l;

    SuperMatrix   A, L, U, B;
    double   *a;
    int_t      *asub, *xa;
    int_t      *perm_r; /* row permutations from partial pivoting */
    int_t      *perm_c; /* column permutation vector */
    SCPformat *Lstore;
    NCPformat *Ustore;
    int_t      nrhs, ldx, info, m, n, nnz, b;
    int_t      nprocs; /* maximum number of processors to use. */
    int_t      panel_size, relax, maxsup;
    int_t      permc_spec;
    trans_t  trans;
    double   *rhs;
    superlu_memusage_t   superlu_memusage;

    nrhs              = 1;
    trans             = NOTRANS;
    nprocs             = 6;
    b                 = 1;
    panel_size        = sp_ienv(1);
    relax             = sp_ienv(2);
    maxsup            = sp_ienv(3);

    /* Initialize matrix A. */
    m = n = Wnum*Lnum*(layer_num) + 2;
    nnz = count;
    if ( !(a = doubleMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for a[]."); // I cannot free the space
    if ( !(asub = intMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for asub[]."); // I cannot free the space
    if ( !(xa = intMalloc(n+1)) ) SUPERLU_ABORT("Malloc fails for xa[]."); // I cannot free the space

    /* assign values to the arrays: a, asub and xa */ 
    int row = -1;
    for (i = 0; i < count; i ++)
    {
        if (Midx[i][0] > row)
        {
            row = Midx[i][0]; // enter a new column 
            xa[row] = i; // index of the first item of each row
        }
        a[i] = Midx[i][2]; 
        asub[i] = (int) Midx[i][1]; // column index of each item 
    }
    xa[row+1] = count; 


    printf("Building the sparse matrix ...\n");
    printf("Dimension of the G matrix is %d x %d\n", m, n);
    printf("Number of non-zero entries is %d\n", nnz); 
    

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    //printf("m = %d\n", m);

   // dPrint_CompCol_Matrix("A", &A);
    /* Create right-hand side matrix B. */
    if ( !(rhs = doubleMalloc(m)) ) SUPERLU_ABORT("Malloc fails for rhs[].");
    // assign values to B 
    for (i = 0; i < m; i++) // initialize rhs to 0
        rhs[i] = 0; 

    for (l = 0; l < layer_num; l ++){
    	for (i = 0; i < Wnum; i ++){
    		for (j = 0; j < Lnum; j ++){
    			rhs[l * Wnum*Lnum + i * Lnum + j] = -1 * powerM[i][j][l] / svdd; 
                //printf("rhs[%d] = -1 * powerM[%d][%d][%d] (%.8f) / %.2f = %.8f\n", l*Wnum*Lnum+i*Lnum+j, i, j, l, powerM[i][j][l], svdd, rhs[l * Wnum*Lnum + i * Lnum + j]);
    		}
    	}
    }
    rhs[layer_num * Wnum*Lnum + 1] = svdd / R_PCB; 

    /*printf("rhs:\n");
    for (i = 0; i < m; i ++)
        printf("%.8f\n", rhs[i]);
    */

    // free the space 
    for (i = 0; i < Wnum; i++)
    {
        for (j = 0; j < Lnum; j++)
        {
            free(powerM[i][j]);
        }
        free(powerM[i]); 
    }
    free(powerM);

    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);


    //dPrint_Dense_Matrix("B", &B);

    if ( !(perm_r = intMalloc(m)) ) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) SUPERLU_ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */     
    permc_spec = 0;
    get_perm_c(permc_spec, &A, perm_c);


    printf("Finish building the sparse matrix\n");
    printf("------------------------------------------------------------\n\n");


    /* Solve the linear system. */
    pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);
    

    printf("Finish solving the linear equation\n");

    //dPrint_Dense_Matrix("B", &B);

    /* extract the Temperature from B */ 
    DNformat *Astore = (DNformat *) B.Store; 
    double *Vt; // vector stores the temperature for all grids 
    double ***V; // matrix stores the temperature for active layers 
    
    Vt = (double *) Astore->nzval; 
    

    // allocate space for T 
    if ( !(V = (double ***)malloc(Wnum * sizeof(double **))) ) printf("Malloc fails for V[].\n");
    for (i = 0; i < Wnum; i++)
    {
        if ( !(V[i] = (double **)malloc(Lnum * sizeof(double *))) ) printf("Malloc fails for V[%d][].\n", i);
        for (j = 0; j < Lnum; j++)
        {
            if ( !(V[i][j] = (double *)malloc(layer_num * sizeof(double))) ) printf("Malloc fails for V[%d][%d][].\n", i,j);
        }
    }
    for (l = 0; l < layer_num; l ++)
        for (i = 0; i < Wnum; i ++)
            for (j = 0; j < Lnum; j++)
                V[i][j][l] = Vt[l * Wnum*Lnum + i * Lnum + j];

    printf("Finish converting the voltage matrix\n");
    printf("Free the space...\n");
    
    if ( info == 0 ) {
    //dinf_norm_error(nrhs, &B, xact); /* Inf. norm of the error */

        Lstore = (SCPformat *) L.Store;
        Ustore = (NCPformat *) U.Store;
        printf("#NZ in factor L = " IFMT "\n", Lstore->nnz);
        printf("#NZ in factor U = " IFMT "\n", Ustore->nnz);
        printf("#NZ in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - L.ncol);
    
        superlu_dQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);
        printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
        superlu_memusage.for_lu/1024/1024, 
        superlu_memusage.total_needed/1024/1024,
        superlu_memusage.expansions);

    }



    /* De-allocate storage */
    // free the arrays defined by myself


    SUPERLU_FREE (rhs);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    printf("finish SUPERLU_FREE\n");
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
    /* De-allocate other storage */
    //free(K); free(H); free(layerP); free(Tt);

    printf("================= FINISH STEADY PDN SOLVER ===============\n\n");

    return V;
}