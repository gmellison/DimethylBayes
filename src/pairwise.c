/*
 * =====================================================================================
 *
 *       Filename:  pairwise.c
 *
 *    Description: G 
 *
 *        Version:  1.0
 *        Created:  11/16/2023 12:59:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
/**
 * @author      : greg (greg@$HOSTNAME)
 * @file        : pairwise
 * @created     : Thursday Nov 16, 2023 12:59:45 EST
 */

#include "bayes.h"
#include "pairwise.h"
#include "utils.h"
#include "command.h"
#include "mcmc.h"
#include "model.h"
#include "likelihood.h"
#include "model.h"


#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#define LIKEPW_EPSILON        1.0e-300
#define LIKEPW_EPSILON        1.0e-300


// global variables for pairwise likelihoods:
// global variables for pw likelihoods:
MrBFlt nucFreqs[4];
MrBFlt doubletFreqs[4][4];

MrBFlt relRates[6];
int    tempNumStates;
int    numTriples;
int I=4;
int J=4;

int defFreqs = NO;
int defDoublets = NO;

MrBFlt alpha=0.1;
MrBFlt dist=0.1;

int *pairwiseCounts;
int **tripleCounts;
int numTrips;
int numPairs;

/* globals declared here */
int defPairwise=NO;
int defTriples=NO;
int allocPairwise=NO;
int allocTriples=NO;

/*  
int usePairwise=NO;
int useTriples=NO;
int useFullForAlpha=NO;
 */


// local prototypes:
int              toIdx(int x);
void             SetupQMat(MrBFlt **Q);
MrBFlt*          AllocateDists(int nPairs);
int***           AllocateDoubletCounts(int nPairs);
void             FreeDoubletCounts(int ***doubletCounts, int np);
int              CalcPairwiseDistsPolyTree(PolyTree *t, PairwiseDists *pd);
int              pairIdx(int i, int j, int n);
int              tripletIdx(int i, int j, int k, int n);
MrBFlt           CalcTripletLogLike_JC(PairwiseDists *pd);
int              FreeDists(MrBFlt* dists);
int              FreePairwiseDists(PairwiseDists* pd);
double           CalcPairwiseLogLike_Full(PairwiseDists *pd);
int              FreeTriples(void);
int              FreePairwise(void);


#define  dIdx( i,j, dim_j ) i*dim_j + j
#define  tIdx( k, i, j, dim_i, dim_j ) k*dim_i*dim_j + i*dim_j + j

int pairIdx(int i, int j, int n) {
    int k;

    if (i == j || i < 0 || j < 0 || i >= n || j >= n) {
        MrBayesPrint("Error in pair indexing; i=%d, j=%d,n=%d \n", i,j,n);
        return(-1);
    }

    if (i > j) {
        k = i;
        i = j;
        j = k;
    } 

    return(n * (n-1)/2 - (n-i)*(n-i-1)/2 + (j-i-1));
}

int tripletIdx(int i, int j, int k, int n) {

    int l;

    if (i == j || j == k || i == k || i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n || n < 3) {
        MrBayesPrint("Error in triplet indexing; i=%d, j=%d,k=%d,n=%d \n", i,j,k,n);
        return(-1);
    }

    /*  if i > j, swap i,j */
    if (i > j) {
        l = i;
        i = j;
        j = l;
    } 

    /*  if i > k, swap i,k so i is for sure the minimum */
    if (i > k) {
        l = i;
        i = k;
        k = l;
    } 

    /*  now if j > k, swap j,k */
    if (j > k) {
        l = j;
        j = k;
        k = l;
    } 

    return(n*(n-1)*(n-2)/6 - (n-i-1)*(n-i-2)*(n-i-3)/6 - (n-j-1)*(n-j-2)/2 - (n-k-1) - 1) ;
}


/*  
 *
 *
 *
 */

int DoEstQPairwise(void) {

    MrBayesPrint("%s Running 'Estqpairwise' command \n", spacer);

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }
           
    if (defFreqs == NO) 
    {
        MrBayesPrint("%s Getting Nucleotide Frequences \n", spacer);
        CountFreqs(); // sets nucleotide frequencies in global object 'nucFreqs'
    }

    if (defDoublets == NO) 
    {
        MrBayesPrint("%s Getting Doublet Counts \n", spacer);
        CountDoublets(numTaxa);
    }

    MrBayesPrint("%s Estimating Q Pairwise  \n", spacer);

    int i,j,k;
    int t1, t2;
    int r1, r2, ridx;
    MrBFlt tau;
    MrBFlt *ratesEst;

    MrBFlt **V;
    MrBFlt **Vinv;      
    MrBFlt **LaLog;     
                       
    MrBComplex **Vc;    
    MrBComplex **Vcinv; 
    MrBFlt la[4]; 
    MrBFlt laC[4];
                       
    MrBFlt **F;            
    MrBFlt **Q;            
    MrBFlt **Temp;         

    int numPairs = numTaxa*(numTaxa-1)/2;
    ratesEst = malloc(6*numPairs*sizeof(MrBFlt));

    for (t1=0; t1<(numTaxa-1); t1++) {
        for (t2=t1+1; t2<numTaxa; t2++) {
 
            // use mb machinery to compute eigens
            // set up matrices 
            V     = AllocateSquareDoubleMatrix(4);
            Vinv  = AllocateSquareDoubleMatrix(4);
            LaLog = AllocateSquareDoubleMatrix(4);

            Vc    = AllocateSquareComplexMatrix(4); 
            Vcinv = AllocateSquareComplexMatrix(4);

            F = AllocateSquareDoubleMatrix(4);
            Q = AllocateSquareDoubleMatrix(4);
            Temp = AllocateSquareDoubleMatrix(4);
       
            tau=0.0;

            k = pairIdx(t1,t2,numTaxa);

            // set up matrix of empirical transitions 
            //   probabilities, 'Q' (will be modified in place)  
            for (i=0;i<4;i++) {
                for (j=0;j<4;j++){
                    F[i][j] = 1.0 * pairwiseCounts[tIdx(k,i,j,I,J)] / (2*numChar*nucFreqs[j]);
                }
            }

            // Symmetrize F:
            for (i=0;i<4;i++) {
                for (j=i;j<4;j++){
                    F[i][j] = (F[i][j]+F[j][i])/2.0;
                    if (i != j)
                            F[j][i] = F[i][j];
                }
            }


            int isComplex=GetEigens(4,F,la,laC,V,Vinv,Vc,Vcinv);
            (void)isComplex;
            
            // diagonal matrix of log^{lambda_i}, lambdas are eigvals of Q
            for (i=0;i<4;i++) {
                LaLog[i][i] = log(la[i]);
                for (j=0;j<i;j++) {
                    LaLog[i][j]=0.0;
                    LaLog[j][i]=0.0;
                }
            }
        
            // calculate the matrix log: log(e^Qt) = V log[La] V^-1) 
            MultiplyMatrices(4,V,LaLog,Temp); 
            MultiplyMatrices(4,Temp,Vinv,Q); // Q is ptr to resulting matrix
       
            // set diagonal so rowsums are 0, mult by inverse Dpi mat:
            double rowsum;
            for (i=0;i<4;i++) {
                rowsum=0.0;
                for (j=0;j<4;j++) {
                    if (j!=i) rowsum+=Q[i][j];
                }
                Q[i][i]=-1*rowsum;
            }
        
            for (i=0;i<4;i++)
                tau += -1.0 * Q[i][i] * nucFreqs[i];
        
            MultiplyMatrixByScalar(4, Q, 1.0/tau, Q);
        
            MrBayesPrint("%s Tau Hat (%d,%d): %f \n", spacer, t1,t2, tau);
        
            // get the individual rates: 
            for (r1=1; r1<4; r1++) {
                for (r2=0; r2<r1; r2++) {
                    ridx = pairIdx(r1,r2,4);
                    ratesEst[dIdx(k,ridx,6)] = Q[r1][r2];
                }
            }

            MrBayesPrint("%s Estimated Relative Rates (pair %d): ", spacer, k);
            for (i=0;i<6;i++)
                MrBayesPrint(" %f", ratesEst[dIdx(k,i,6)]);
            MrBayesPrint("\n");

            // free all matrices 
            FreeSquareDoubleMatrix(V);
            FreeSquareDoubleMatrix(Vinv);
            FreeSquareDoubleMatrix(Q);
            FreeSquareComplexMatrix(Vc);
            FreeSquareComplexMatrix(Vcinv);
            FreeSquareDoubleMatrix(Temp);
            FreeSquareDoubleMatrix(F);
        }
    }


    /*  average the rates across pairs */
    MrBFlt r6;
    MrBFlt ratesOut[6] = {0.0}; 

    
    /*  Normalize the rates:  */
    for (k=0; k<numPairs; k++) {

        r6=ratesEst[dIdx(k,5,6)];
        for (i=0;i<6;i++) {
            ratesEst[dIdx(k,i,6)]/=r6;
            ratesOut[i]+=ratesEst[dIdx(k,i,6)]/(1.0*numPairs);
        }
    }

    MrBayesPrint("%s Estimated Relative Rates (Averaged): ", spacer);
    for (i=0;i<6;i++)
        MrBayesPrint(" %f", ratesOut[i]);
    MrBayesPrint("\n");

    free(ratesEst);

    return(NO_ERROR);
}



double CalcPairwiseLogLike_Reduced(PairwiseDists *pd) {

    int i,j,k;
    int ri;
    MrBFlt tp;
    MrBFlt ll=0.0;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    /* Compute average pairwise dist:  */
    MrBFlt dist = 0.0;
    for (k=0;k<pd->nPairs;k++) 
        dist += (pd->dists)[k];

    dist = dist/(1.0*pd->nPairs);


    /*  pool doublet counts:  */
    int **doubletCountsPooled;
    doubletCountsPooled = AllocateSquareIntegerMatrix(4);

    for (k=0; k<pd->nPairs; k++) {
        for (i=0; i<4; i++) 
            for (j=0; j<4; j++) 
                doubletCountsPooled[i][j] += pairwiseCounts[tIdx(k,i,j,I,J)];
    }
  
    // compute gtr transition probabilities
    // set up matrices for taking the matrix exponent
    MrBFlt **V         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Vinv      = AllocateSquareDoubleMatrix(4);
    MrBFlt **LaExp     = AllocateSquareDoubleMatrix(4);
    MrBFlt **TempMat   = AllocateSquareDoubleMatrix(4);

    MrBComplex **Vc    = AllocateSquareComplexMatrix(4);    // eigenvectors
    MrBComplex **Vcinv = AllocateSquareComplexMatrix(4);    // inverse eigvect matrix
    MrBFlt *la         = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // eigenvalues
    MrBFlt *laC        = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // 

    MrBFlt **Q         = AllocateSquareDoubleMatrix(4);
    MrBFlt **Qtausr    = AllocateSquareDoubleMatrix(4);

    // allocate array of 4x4 matrices -- one for each rate category
    MrBFlt **pTrans[100];
    for (int ri=0; ri<4; ri++)
            pTrans[ri]=AllocateSquareDoubleMatrix(4);

    SetupQMat(Q);

    // calculate the site pattern probability for each 
    //   dg4 rate category
    for (ri=0; ri<4; ri++) {
    
        // probability transition matrix for site rate i:
        MultiplyMatrixByScalar(4, Q, dist*dg4[ri], Qtausr);  

        int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
        if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!! \n");
           
        // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
        for (i=0;i<4;i++) {
            LaExp[i][i] = exp(la[i]);
            for (j=0;j<i;j++) {
                LaExp[i][j]=0.0;
                LaExp[j][i]=0.0;
            }
        }
           
        // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
        MultiplyMatrices(4,V,LaExp,TempMat); 
        MultiplyMatrices(4,TempMat,Vinv,pTrans[ri]);
    }

 
    for (i=0; i<4; i++) 
        for (j=0; j<4; j++) 
            MrBayesPrint("%s Q(%d,%d) %f \n", spacer, i,j, pTrans[1][i][j]);
    
    // now sum over all the doublet patterns
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
         
            // reset tp, and sum up the conditional 
            //   trans probs over discrete rate categories
            tp=0.0; 
            for (int ri=0; ri<4; ri++) {
                tp += 0.25 * pTrans[ri][i][j];
            }

            // calculate triple probability by summing over 
            //   conditional probs given rate categories 
            ll += doubletCountsPooled[i][j] * log(tp);
            
        }
    }

    FreeSquareDoubleMatrix(TempMat);
    for (int ri=0;ri<4;ri++)
        FreeSquareDoubleMatrix(pTrans[ri]);


    return ll;

}

double CalcPairwiseLogLike_Full(PairwiseDists *pd) {

    int i,j;
    int ri;
    int t1, t2, k;
    MrBFlt tp;
    MrBFlt ll=0.0;
    MrBFlt dist;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    MrBFlt **Q = AllocateSquareDoubleMatrix(4);
    SetupQMat(Q);

    for (t1=1; t1<pd->nTaxa; t1++)
        {
        for (t2=0; t2<t1; t2++)
            {

            k = pairIdx(t1,t2,pd->nTaxa);    
            dist = pd->dists[k];

            // compute gtr transition probabilities
            // set up matrices for taking the matrix exponent
            MrBFlt **V         = AllocateSquareDoubleMatrix(4);
            MrBFlt **Vinv      = AllocateSquareDoubleMatrix(4);
            MrBFlt **LaExp     = AllocateSquareDoubleMatrix(4);
            MrBFlt **TempMat   = AllocateSquareDoubleMatrix(4);

            MrBComplex **Vc    = AllocateSquareComplexMatrix(4);    // eigenvectors
            MrBComplex **Vcinv = AllocateSquareComplexMatrix(4);    // inverse eigvect matrix
            MrBFlt *la         = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // eigenvalues
            MrBFlt *laC        = (MrBFlt*)malloc(sizeof(MrBFlt)*4); // 

            MrBFlt **Qtausr    = AllocateSquareDoubleMatrix(4);

            // allocate array of 4x4 matrices -- one for each rate category
            MrBFlt **pTrans[100];
            for (int ri=0; ri<4; ri++)
                    pTrans[ri]=AllocateSquareDoubleMatrix(4);

            // calculate the site pattern probability for each 
            //   dg4 rate category
            for (ri=0; ri<4; ri++) {
            
                // probability transition matrix for site rate i:
                MultiplyMatrixByScalar(4, Q, dist*dg4[ri], Qtausr);  

                int isComplex=GetEigens(4,Qtausr,la,laC,V,Vinv,Vc,Vcinv);
                if (isComplex) MrBayesPrint("Complex Eigens in Eidendecomp!! \n");
                   
                // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
                for (i=0;i<4;i++) {
                    LaExp[i][i] = exp(la[i]);
                    for (j=0;j<i;j++) {
                        LaExp[i][j]=0.0;
                        LaExp[j][i]=0.0;
                    }
                }
                   
                // diagonal matrix of e^{lambda_i}, lambdas are eigvals of Q
                MultiplyMatrices(4,V,LaExp,TempMat); 
                MultiplyMatrices(4,TempMat,Vinv,pTrans[ri]);
            }
       
            // now sum over all the doublet patterns
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                 
                    // reset tp, and sum up the conditional 
                    //   trans probs over discrete rate categories
                    tp=0.0; 
                    for (int ri=0; ri<4; ri++) {
                        tp += 0.25 * pTrans[ri][i][j];
                    }

                    // calculate triple probability by summing over 
                    //   conditional probs given rate categories 
                    MrBayesPrint("%s k=%d, i=%d, j=%d, tIdx: %d, Count: %d, log(p): %f \n", 
                            spacer,
                            k,i,j,
                            tIdx(k,i,j,I,J), 
                            pairwiseCounts[tIdx(k,i,j,I,J)],  
                            log(tp));
                    ll += pairwiseCounts[tIdx(k,i,j,I,J)] * log(tp);
                    
                }

            }

            FreeSquareDoubleMatrix(TempMat);
            for (int ri=0;ri<4;ri++)
                FreeSquareDoubleMatrix(pTrans[ri]);
            }
        }


    return ll;

}


int DoPairwiseLogLike(void) {

    MrBayesPrint("Running 'Pairwiseloglike' command.\n");

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }

    if (defFreqs==NO) 
    {
        CountFreqs();
    }

    if (defDoublets==NO) 
    {
        CountDoublets(numTaxa);
        MrBayesPrint("Done Counting Doublets \n");
    }

    if (numUserTrees==0) {
        MrBayesPrint ("%s   No user trees defined! \n", spacer);
        return (ERROR);
    }

    int l;
    MrBFlt ll;
    PairwiseDists *pd; 

    MrBayesPrint ("%s   Calculating Pw likelihood for %d defined user trees. \n", spacer, numUserTrees);

    for (l=0; l<numUserTrees; l++) {
        pd = AllocatePairwiseDists();
        InitPairwiseDistsPolyTree(userTree[l],pd);

        ll=CalcPairwiseLogLike_Full(pd);
        MrBayesPrint("Pairwise Log-Likelihood: %f \n", ll);
        FreePairwiseDists(pd);
    }


    return(NO_ERROR);

}

/* Inits a PairwiseDists struct from a polytree  
 *
 *  
 */

PairwiseDists* AllocatePairwiseDists(void) {
    PairwiseDists *pd;
    pd=(PairwiseDists*)SafeCalloc(1,sizeof(PairwiseDists));
    return(pd);
}

/*
void InitPairwiseDists(Tree *t, PairwiseDists *pd) 
{
    pd->nTaxa=(t->nNodes-t->nIntNodes);
    pd->nPairs=(pd->nTaxa)*(pd->nTaxa-1)/2;
    pd->dists=AllocateDists(pd->nTaxa);
    CalcPairwiseDists_ReverseDownpass(t, pd);
}
*/

void InitPairwiseDistsPolyTree(PolyTree *pt, PairwiseDists *pd) 
{
    pd->nTaxa=(pt->nNodes-pt->nIntNodes);
    pd->nPairs=(pd->nTaxa)*(pd->nTaxa-1)/2;
    pd->dists=AllocateDists(pd->nTaxa);
    CalcPairwiseDistsPolyTree(pt, pd);
}

int FreePairwiseDists(PairwiseDists* pd) 
{
    FreeDists(pd->dists);
    free(pd);
    return(NO_ERROR);
}

int CalcPairwiseDists_ReverseDownpass(Tree *t, int division, int chain)
{
    int         i,j,a,d,k;
    TreeNode    *p;
    double      x;
    MrBFlt      *dists, **distsTemp;
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    dists = m->pwDists[chain];
    int numExtNodes = t->nNodes - t->nIntNodes;

    /*  We'll calculate (in distsTemp) all node dists including internal nodes  */
    distsTemp=(MrBFlt**)malloc(t->nNodes*sizeof(MrBFlt*)); 
    for (k=0;k<t->nNodes;k++)
            distsTemp[k]=(MrBFlt*)malloc(t->nNodes*sizeof(MrBFlt));

    /*  make sure dists are init to 0  */
    for (i=0; i<t->nNodes; i++) 
        for (j=0; j<t->nNodes; j++)
                distsTemp[i][j]=0.0;

    /* loop over all nodes  */
    for (i=(t->nNodes)-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        a = p->anc->index;   
        d = p->index;
        x = p->length;

        /*  start with distance from node to ancestor  */
        distsTemp[d][a] = x;
        distsTemp[a][d] = x;

        /* now revisit previously visited nodes, updating distances  */ 
        for (j=i+1; j<t->nNodes; j++)
            {
                k=t->allDownPass[j]->index;
                if (k==a) 
                    continue;

                /*  dist from current node to prior node =  
                 *      dist from anc to prior node + dist from anc to current node  */
                distsTemp[d][k] = distsTemp[a][k] + x;
                distsTemp[k][d] = distsTemp[k][a] + x;

            }
        }

    /*  set the distances in the output array */
    for (i=0; i<(numExtNodes-1); i++) {
        for (j=i+1; j<numExtNodes; j++) {
            dists[pairIdx(i,j,numExtNodes)]=distsTemp[i][j];
        }
    }

    /*  free temp array and return pointer to taxa pairwise distances */
    for (k=0;k<t->nNodes;k++)
        free(distsTemp[k]);
    free(distsTemp);

    return(NO_ERROR);
}



int CalcPairwiseDistsPolyTree(PolyTree *t, PairwiseDists *pd)
{
    int         i,j,a,d,k;
    PolyNode    *p;
    double      x;

    int numExtNodes = t->nNodes - t->nIntNodes;

    /*  We'll calculate (in distsTemp) all node dists including internal nodes  */
    MrBFlt **distsTemp;
    distsTemp=(MrBFlt**)malloc(t->nNodes*sizeof(MrBFlt*)); 
    for (k=0;k<t->nNodes;k++)
            distsTemp[k]=(MrBFlt*)malloc(t->nNodes*sizeof(MrBFlt));

    /*  make sure dists are init to 0  */
    for (i=0; i<t->nNodes; i++) 
        for (j=0; j<t->nNodes; j++)
                distsTemp[i][j]=0.0;

    /* loop over all nodes  */
    for (i=1; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        a = p->anc->index;   
        d = p->index;
        x = p->length;

        /*  start with distance from node to ancestor  */
        distsTemp[d][a] = x;
        distsTemp[a][d] = x;

        /* now revisit previously visited nodes, updating distances  */ 
        for (j=i-1; j>=0; j--)
            {
                k=(&t->nodes[j])->index;
                if (k==a) 
                { 
                    continue;
                }

                /*  dist from current node to prior node =  
                 *      dist from anc to prior node + dist from anc to current node  */
                distsTemp[d][k] = distsTemp[a][k] + x;
                distsTemp[k][d] = distsTemp[k][a] + x;

            }
        }

    /*  set the distances in the output array */
    for (i=0; i<(numExtNodes-1); i++) {
        for (j=i+1; j<numExtNodes; j++) {
            pd->dists[pairIdx(i,j,pd->nTaxa)]=distsTemp[i][j];
        }
    }

    /*  free temp array and return pointer to taxa pairwise distances */
    for (k=0;k<t->nNodes;k++)
        free(distsTemp[k]);
    free(distsTemp);

    return(NO_ERROR);
}


void CountFreqs(void) {

    int s,i,j,nucIdx;
    int nucCounts[4] = {0};

    for (s=0;s<numChar;s++) {
        for (i=0;i<numTaxa;i++) {
            if (matrix[pos(i,s,numChar)]==GAP)
                continue;

            nucIdx=toIdx(matrix[pos(i,s,numChar)]);
            nucCounts[nucIdx]++;
        }
    }

    for (j=0;j<4;j++) {
        nucFreqs[j]=((MrBFlt)nucCounts[j])/(numChar*numTaxa);
    }

    defFreqs=YES;
}


int CountDoublets(int nTaxa) {

    int s,i,j,k;
    int id1, id2;
 
    int I=4, J=4; 

    int nPairs = nTaxa*(nTaxa-1)/2;

    if (defMatrix == NO) 
    {
        MrBayesPrint("%s Matrix not Defined! \n", spacer);
    }

    if (defDoublets == YES) 
    {
        MrBayesPrint("%s Recomputing Doublets \n", spacer);
        free(pairwiseCounts);
    }

    pairwiseCounts = malloc(4*4*nPairs*sizeof(int));
    if (pairwiseCounts == NULL)
        {
        MrBayesPrint("%s Error Allocating Doublet Counts \n",spacer);
        return(ERROR);
        }

    // initialize doublet counts to 0
    for (k=0; k<nPairs; k++) 
        for (i=0; i<4; i++) 
             for (j=0; j<4; j++) {
                pairwiseCounts[tIdx(k,i,j,I,J)]=0;
             }

    for (i=1;i<nTaxa;i++) 
        {
        for (j=0;j<i;j++) 
            {
            k=pairIdx(i,j,nTaxa);

            for (s=0;s<numChar;s++) 
                {
                if (matrix[pos(i,s,numChar)]==GAP | matrix[pos(j,s,numChar)]==GAP)
                    continue;

                /* nucleotides at position x of sequences i & j   */        
                id2=toIdx(matrix[pos(j,s,numChar)]);
                id1=toIdx(matrix[pos(i,s,numChar)]);

                /* increment the count of that nucleotide pair, at pair k  */
                pairwiseCounts[tIdx(k,id1,id2,I,J)]++;
                }
            }
        }

    defDoublets = YES;

    return(NO_ERROR);
}

int toIdx(int x){
    if (x == 1) return 0;
    else if (x == 2) return 1;
    else if (x == 4) return 2;
    else if (x == 8) return 3;
    else if (x == GAP) return 5;
    else {
        MrBayesPrint("Problem in to Idx (maybe bad input: x=%d?)\n", x);
        return ERROR;
    }
}

int triplePos(int i, int j, int k) {
    return i*16 +j*4 + k;
}

void SetupQMat(MrBFlt **Q) {

    int i, j;

    if (defFreqs == NO) 
        {
        CountFreqs();
        }

    Q[0][1] = Q[1][0] = relRates[0];
    Q[0][2] = Q[2][0] = relRates[1];
    Q[0][3] = Q[3][0] = relRates[2];
    Q[1][2] = Q[2][1] = relRates[3];
    Q[1][3] = Q[3][1] = relRates[4];
    Q[2][3] = Q[3][2] = relRates[5];


    for (i=0;i<4;i++) {
        for (j=0;j<4;j++) {
            Q[i][j] = Q[i][j] * nucFreqs[j];
        }
    }

    // set diagonal elements: 
    double rowsum;
    for (i=0;i<4;i++) {
        rowsum=0.0;
        for (j=0;j<4;j++) {
            if (j!=i) rowsum+=Q[i][j];
        }
        Q[i][i]=-1*rowsum;
    }
} 

int DoPwSetParm (char *parmName, char *tkn)
{
    /*  char        *tempStr;  */
    MrBFlt      tempD;
    int j;

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Autoclose (autoClose) **********************************************************/
        if (!strcmp(parmName, "Dist"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                dist = tempD;
                MrBayesPrint ("%s   Setting distance for pw likelihood calcs to %lf\n", spacer, dist);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }        /* set Nowarnings (noWarn) **********************************************************/
        else if (!strcmp(parmName, "Alpha"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                alpha = tempD;
                MrBayesPrint ("%s   Setting alpha for pw likelihood calcs to %lf\n", spacer, alpha);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                return (ERROR);
                }
            }        /* set Nowarnings (noWarn) **********************************************************/


        /* set Revmatpr (revMatPr) *********************************************************/
        else if (!strcmp(parmName, "Relrates"))
            {

            if (expecting == Expecting(EQUALSIGN)) 
                {
                expecting = Expecting(LEFTPAR);
                }

            else if (expecting == Expecting(LEFTPAR))
                {
                if (1) // TODO: add relrates to isargvalid (IsArgValid(tkn, tempStr) == NO_ERROR)

                    {
                    relRates[0] = relRates[1] = 1.0;
                    relRates[2] = relRates[3] = 1.0;
                    relRates[4] = relRates[5] = 1.0;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid Relrates argument\n", spacer);
                    return (ERROR);
                    }

                expecting  = Expecting(NUMBER);
                tempNumStates = 0;
                }

            else if (expecting == Expecting(NUMBER))
                {
                /* find out what type of prior is being set */
                /* find and store the number */
                sscanf (tkn, "%lf", &tempD);
                tempNum[tempNumStates++] = tempD;

                if (tempNumStates == 1)
                    expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                else if (tempNumStates < 6)
                    expecting  = Expecting(COMMA);
                else
                    expecting = Expecting(RIGHTPAR);
                }

            else if (expecting == Expecting(COMMA))
                {
                expecting  = Expecting(NUMBER);
                }

            else if (expecting == Expecting(RIGHTPAR))
                {
                for (j=0; j<6; j++)
                    {
                    if (tempNumStates == 1)
                        relRates[j] = tempNum[0] / (MrBFlt) 6.0;
                    else
                        relRates[j] = tempNum[j];
                    }

                MrBayesPrint ("%s   Setting Relrates to (%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf,%1.2lf)\n", spacer, 
                    relRates[0], relRates[1], relRates[2],
                    relRates[3], relRates[4], relRates[5]);

                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);

                }
            else
                return (ERROR);

            } // closes parsing relrates param

        /* set Quitonerror (quitOnError) **************************************************/
        else
            {
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


MrBFlt* AllocateDists(int nDists) {

    MrBFlt* dists;
    dists=(MrBFlt*)malloc(nDists*sizeof(MrBFlt));
    return (dists);
}

int FreeDists(MrBFlt* dists) {

    free(dists);
    return (NO_ERROR);
}

int*** AllocateDoubletCounts(int nPairs) {

    int i;
    int ***doubletCounts=(int***)malloc(nPairs*sizeof(int**));
    for (i=0; i<nPairs; i++) {
        doubletCounts[i]=AllocateSquareIntegerMatrix(4);
    }
    return (doubletCounts);
}


void FreeDoubletCounts(int ***doubletCounts, int np) {
    int i;
    for (i=0; i<np; i++) {
        FreeSquareIntegerMatrix(doubletCounts[i]);
    }
    free(doubletCounts);
}

int DoTripletLogLike(void) {

    if (defMatrix == NO)
    {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
    }

    if (numUserTrees == 0)
    {
        MrBayesPrint ("%s   No user trees defined. \n", spacer);
        return (ERROR);
    }

    if (defFreqs==NO) 
    {
        CountFreqs();
    }

    if (defTriples==NO)
    {
        CountTriplets();
    }

    PairwiseDists* pd;
    int ti;
    MrBFlt ll=0.0;

    for (ti=0; ti<numUserTrees; ti++) {
        pd = AllocatePairwiseDists();
        InitPairwiseDistsPolyTree(userTree[ti],pd);

        ll = CalcTripletLogLike_JC(pd);
        MrBayesPrint("Triplet Quasi Log-Likelihood for Tree %d: %f", ti, ll);
        FreePairwiseDists(pd);
    }

    return(NO_ERROR);

}



/* PrintNodes: Print a list of tree nodes, pointers and length */
void PrintPairwiseDists (PairwiseDists *pd)
{
    int i,j;
    int n;

    n = pd->nTaxa;
    /* printf ("tip1\ttip2\tdist\n"); */
    for (i=1; i<n; i++)
        {
        for (j=0; j<i; j++)
            {
            printf ("%d\t%d\t%f\n",
              i,j,
              pd->dists[pairIdx(i,j,n)]);
            }
        }

    printf ("\n");
}

MrBFlt CalcTripletLogLike_JC(PairwiseDists *pd) {

    double al=0.5;
    int i,j,k,ri;
    int seq1, seq2, seq3;
    int pairidx, tripidx, nucidx;
     

    MrBFlt dg4[4];
    DiscreteGamma(dg4,al,al,4,1);

    MrBFlt pii[3][4];
    MrBFlt pij[3][4];

    MrBFlt tripleProbs[64];
    MrBFlt transitionProbs[3][64];

    MrBFlt pwdist[3];
    MrBFlt cndist[3];
 
    MrBFlt probByRateCat[4];
    MrBFlt ll=0.0;
 
    /* 
     * loop over alignment sequence triples 
     */ 
    for (seq1=0; seq1<(pd->nTaxa-2); seq1++) {
        for (seq2=seq1+1; seq2<(pd->nTaxa-1); seq2++) {
            for (seq3=seq2+1; seq3<pd->nTaxa; seq3++) {

                tripidx = tripletIdx(seq1,seq2,seq3,pd->nTaxa);

                /*
                 * get the distances between triple taxa (pw and cn)
                 */
                pwdist[0]=pd->dists[pairIdx(seq1,seq2,pd->nTaxa)];
                pwdist[1]=pd->dists[pairIdx(seq1,seq3,pd->nTaxa)];
                pwdist[2]=pd->dists[pairIdx(seq2,seq3,pd->nTaxa)];
                MrBayesPrint("Triple Dists: %f, %f, %f \n", pwdist[0], pwdist[1], pwdist[2]);

                cndist[0] = (pwdist[0] + pwdist[1] - pwdist[2])/2.0;
                cndist[1] = (pwdist[0] + pwdist[2] - pwdist[1])/2.0;
                cndist[2] = (pwdist[1] + pwdist[2] - pwdist[0])/2.0;
                MrBayesPrint("Centernode Dists: %f, %f, %f \n", cndist[0], cndist[1], cndist[2]);

                /*  
                 *  calc JC transition probs for each centernode dist
                 */
                for (ri=0; ri<4; ri++) {
                    for (pairidx=0; pairidx<3; pairidx++) {
                        MrBFlt etr = exp(-(4.0/3) * cndist[pairidx] * dg4[ri]);
                        pii[pairidx][ri]=(1.0/4) + (3.0/4) * etr;
                        pij[pairidx][ri]=(1.0/4) - (1.0/4) * etr;

                        /*
                        MrBayesPrint("  Pij[%d][%d]: %f", pairidx, ri, pij[pairidx][ri]);
                        MrBayesPrint("  Pii[%d][%d]: %f", pairidx, ri, pii[pairidx][ri]);
                         */

                        // set up array of nucleotide transition probabilities                                
                        for (i=0;i<4;i++) {
                            for (j=0;j<4;j++) {
                                if (i==j) {
                                    transitionProbs[pairidx][tIdx(i,j,ri,4,4)]=pii[pairidx][ri];
                                } else {
                                    transitionProbs[pairidx][tIdx(i,j,ri,4,4)]=pij[pairidx][ri];
                                }
                            }
                        }
                    }
                }

                /*  
                *  calc JC transition probs for each centernode dist
                */
                for (i=0; i<4; i++) {
                    for (j=0; j<4; j++) {
                        for (k=0; k<4; k++) {

                            nucidx = tIdx(i,j,k,4,4);

                            // no need to calculate probs for not present in data
                            if (tripleCounts[tripidx][nucidx] == 0) continue;  

                            // ensure triple site pattern prob is init to 0:
                            tripleProbs[nucidx] = 0.0;

                            // calculate the site pattern probability for each 
                            //   dg4 rate category
                            for (ri=0; ri<4; ri++) {

                                // reset helper rate cat array at current index
                                probByRateCat[ri]=0.0;

                                // sum over the 4 possible internal nucleotides, 
                                // given rate cat ri:   
                                for (int nucx=0; nucx<4; nucx++) {
                                        probByRateCat[ri]+=
                                                    nucFreqs[nucx]
                                                    *transitionProbs[0][tIdx(nucx,i,ri,4,4)]
                                                    *transitionProbs[1][tIdx(nucx,j,ri,4,4)]
                                                    *transitionProbs[2][tIdx(nucx,k,ri,4,4)];
                                }

                                // calculate triple probability by summing over 
                                //   conditional probs given rate categories 
                                tripleProbs[nucidx]+=0.25*probByRateCat[ri];
                            }

                            ll += tripleCounts[tripidx][nucidx] * log(tripleProbs[nucidx]);
                        }
                    }
                }
            }
        }
    }
    return(ll);
}

  
MrBFlt CalcTripletLogLike_JC_Reduced(PairwiseDists *pd, MrBFlt alpha) {


    MrBFlt tripleProbs[64], transitionProbs[64];

    int i,j,k, ri;
    int ti;
    int spIdx;

    MrBFlt dist;

    MrBFlt dg4[4];
    DiscreteGamma(dg4,alpha,alpha,4,1);

    MrBFlt pii[4];
    MrBFlt pij[4];

    /*  Calculate avg dist: 
     */
    dist = 0.0;
    for (k=0;k<pd->nPairs;k++) 
        dist += (pd->dists)[k];

    dist = dist/(2.0*pd->nPairs);


    /*  pool triplet counts:  */
    int tripleCountsPooled[64];

    for (ti=0; ti<pd->nPairs; ti++)
        for (i=0; i<4; i++) 
            for (j=0; j<4; j++) 
                for (k=0;k<4;k++) 
                    tripleCountsPooled[tIdx(i,j,k,4,4)] += tripleCounts[ti][tIdx(i,j,k,4,4)];
    
 
    for (int i=0; i<4; i++) {
        MrBFlt etr = exp(-(4.0/3) * dist* dg4[i]);
        pii[i]=(1.0/4) + (3.0/4) * etr;
        pij[i]=(1.0/4) - (1.0/4) * etr;
    }

    // set up array of transition probabilities
    for (i=0;i<4;i++) {
            for (j=0;j<4;j++) {
                    for (int ri=0;ri<4;ri++) {
                            if (i==j) {
                                transitionProbs[triplePos(i,j,ri)]=pii[ri];
                            } else {
                                transitionProbs[triplePos(i,j,ri)]=pij[ri];
                            }
                    }
            }
    }

    MrBFlt probsByRateCat[4];
    MrBFlt ll=0.0;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {

                spIdx=tIdx(i,j,k,4,4);

                // no need to calculate probs for not present in data
                if (tripleCounts[spIdx] == 0) continue;  

                // ensure triple site pattern prob is init to 0:
                tripleProbs[spIdx] = 0.0;

                // calculate the site pattern probability for each 
                //   dg4 rate category
                for (ri=0; ri<4; ri++) {

                    // reset helper rate cat array at current index
                    probsByRateCat[ri]=0.0;

                    // sum over the 4 possible internal nucleotides, 
                    // given rate cat ri:   
                    for (int nucx=0; nucx<4; nucx++) {
                            probsByRateCat[ri]+=
                                        nucFreqs[nucx]
                                        *transitionProbs[triplePos(nucx,i,ri)]
                                        *transitionProbs[triplePos(nucx,j,ri)]
                                        *transitionProbs[triplePos(nucx,k,ri)];
                    }

                    // calculate triple probability by summing over 
                    //   conditional probs given rate categories 
                    tripleProbs[spIdx]+=0.25*probsByRateCat[ri];
                }

                ll += tripleCountsPooled[spIdx] * log(tripleProbs[spIdx]);
            }
        }
    }
    return (ll);
}



/*-----------------------------------------------------------------
|
|   TiProbs_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsPairwise_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    
    int         i, j, k, p, index;
    MrBFlt      *dists; 
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange;
    CLFlt       *tiP;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    /* find pw dists and transition probabilities */
    dists = m->pwDists[chain];

    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
   
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    for (p=0; p<numPairs; p++)
        {

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];

        length = dists[p];

        /* fill in values */
        index=0;
        for (k=0; k<m->numRateCats; k++)
            {

            t = length * baseRate * catRate[k];

            /* numerical errors will ensue if we allow very large or very small branch lengths,
            which might occur in relaxed clock models */
            if (t < TIME_MIN)
                {
                /* Fill in identity matrix */
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = 1.0;
                        else
                            tiP[index++] = 0.0;
                        }
                    }
                }
            else if (t > TIME_MAX)
                {
                /* Fill in stationary matrix */
                for (i=0; i<4; i++)
                    for (j=0; j<4; j++)
                        tiP[index++] = 0.25;
                }
            else
                {
                /* calculate probabilities */
                pChange   = (CLFlt) (0.25 - 0.25 * exp(-(4.0/3.0)*t));
                pNoChange = (CLFlt) (0.25 + 0.75 * exp(-(4.0/3.0)*t));

                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = pNoChange;
                        else
                            tiP[index++] = pChange;

                        }
                    }
                }
            }
        } /*  end loop over pairs */

        /*  Pw transition probabilities are done, now let's fill in the doublet probs: */
    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   TiProbsPairwise_Gen: update transition probabilities for 4by4
|       nucleotide model with nst == 6 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsPairwise_Gen (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    register int    i, j, k, n, p, s, index;
    MrBFlt          t, *catRate, baseRate, *eigenValues, *cijk, *bs,
                    EigValexp[64], sum, *ptr, theRate,
                    length;
    CLFlt           *tiP;
    ModelInfo       *m;
    MrBFlt          *dists; 

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    n = m->numModelStates; 

    /* find pw dists and transition probabilities */
    dists = m->pwDists[chain];

    /* get base rate */
    baseRate = GetRate (division, chain);
   
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
   
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    /* get eigenvalues and cijk pointers */
    eigenValues = m->cijks[m->cijkIndex[chain]];
    cijk        = eigenValues + (2 * n);

    for (p=0; p<numPairs; p++)
        {

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        length = dists[p];

        /* fill in values */
        index=0;
        for (k=0; k<m->numRateCats; k++)
            {

            t = length * baseRate * catRate[k];

            /* numerical errors will ensue if we allow very large or very small branch lengths,
            which might occur in relaxed clock models */
            if (t < TIME_MIN)
                {
                /* Fill in identity matrix */
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        if (i == j)
                            tiP[index++] = 1.0;
                        else
                            tiP[index++] = 0.0;
                        }
                    }
                }
            else if (t > TIME_MAX)
                {
                /* Get base freq */
                bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

                /* Fill in stationary matrix */
                for (i=0; i<4; i++)
                    for (j=0; j<4; j++)
                        tiP[index++] = (CLFlt) bs[j];
                }
            else
                {
                /* We actually need to do some work... */
                for (s=0; s<4; s++)
                    EigValexp[s] =  exp(eigenValues[s] * t);

                ptr = cijk;
                for (i=0; i<4; i++)
                    {
                    for (j=0; j<4; j++)
                        {
                        sum = 0.0;
                        for (s=0; s<4; s++)
                            sum += (*ptr++) * EigValexp[s];
                        tiP[index++] = (CLFlt) ((sum < 0.0) ? 0.0 : sum);
                        }
                    }
                }
            }
        } /*  end loop over pairs */
        /*  Pw transition probabilities are done, now let's fill in the doublet probs: */
    return (NO_ERROR);
}

/*-----------------------------------------------------------------
|
|   DoubletProbabilities_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/

int DoubletProbs_JukesCantor(int division, int chain)
{

    
    int         i, j, k, p, index, dpIdx;
    CLFlt       *tiP, *doubP;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    for (p=0; p<numPairs; p++)
        {

        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        /*  reset doublet probs to 0: */
        dpIdx = 0;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                doubP[dpIdx++] = 0.0;
        
        index=0; 
        for (k=0; k<m->numRateCats; k++) 
            {
            dpIdx=0;
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    doubP[dpIdx++] += 0.25 * tiP[index++]/((MrBFlt)m->numRateCats);
            }
        }
    return(NO_ERROR);
}


int DoubletProbs_Gen(int division, int chain)
{

    int         i, j, k, p, index, dpIdx;
    CLFlt       *tiP, *doubP;   
    ModelInfo   *m;
    MrBFlt      *bs;

    m = &modelSettings[division];
    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    for (p=0; p<numPairs; p++)
        {
        tiP = m->tiProbsPw[m->pwIndex[chain][p]];
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        /*  reset doublet probs to 0: */
        dpIdx = 0;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                doubP[dpIdx++] = 0.0;
        
        index=0; 
        for (k=0; k<m->numRateCats; k++) 
            {
            dpIdx=0;
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    doubP[dpIdx++] += bs[i] * tiP[index++]/((MrBFlt)m->numRateCats);
                    }
                }
            }
        }

    return(NO_ERROR);
}




int CountPairwise(void) {

    int             i,j,k,c,id1,id2;
    
    /*  For now, only inplemented for a single partition 
     *  of DNA type data */
    /*  
    MrBayesPrint("Counting Pairwise\n"); 
    if (mp->dataType!=DNA && mp->dataType != RNA)
        { 
        MrBayesPrint("%s Can't count pairwise site-patterns for non DNA data  \n", spacer);
        return(ERROR);
        }
    */

    MrBayesPrint("Counting Pairwise\n");
    if (defMatrix == NO) 
        {
        MrBayesPrint("%s Matrix needs to be defined before counting doublets  \n", spacer);
        return(ERROR);
        }  
    
    if (numDefinedPartitions > 1)    
        { 
        MrBayesPrint("%s Pairwise count likelihood only implemented for a single partition  \n", spacer);
        return(ERROR);
        }

    if (memAllocs[ALLOC_PAIRWISE] == YES) 
        {
        free(pairwiseCounts);
        pairwiseCounts=NULL;
        memAllocs[ALLOC_PAIRWISE] = NO;
        }

    numPairs=(int)numTaxa*(numTaxa-1)/2;
    MrBayesPrint("%s Using pairwise likelihood with %d pairs. \n", spacer, numPairs);

    /* first allocate pairwise doublet counts:  */
    pairwiseCounts=(int*)SafeMalloc(numPairs * 16 * sizeof(int));
    if (!pairwiseCounts)
        {
        MrBayesPrint("%s Problem allocating pairwise counts! \n", spacer);
        free(pairwiseCounts);
        return(ERROR);
        }
    memAllocs[ALLOC_PAIRWISE]=YES;

    /* now count the doublets across taxa pairs */
    for (i=0; i<(numTaxa-1); i++)
        {
        for (j=i+1; j<numTaxa; j++)
            {
            k=pairIdx(i,j,numTaxa);
            for (c=0;c<numChar;c++)
                {
                if (matrix[pos(i,c,numChar)]==GAP | matrix[pos(j,c,numChar)]==GAP)
                    continue;

                 /* nucleotides at position x of sequences i & j   */        
                id2=toIdx(matrix[pos(j,c,numChar)]);
                id1=toIdx(matrix[pos(i,c,numChar)]);

                /* increment the count of that nucleotide pair, at pair k  */
                pairwiseCounts[tIdx(k,id1,id2,4,4)]++;

               }
            }
        }

    defPairwise=YES;
    return (NO_ERROR);
}

int PrepareHybridStep(int chain)
{
    ModelInfo   *m;
    Tree        *tree;
    TreeNode    *p;
    int         i,d;

    d = 0;
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);

    m->upDateCijk=YES;


    for (i=0; i<tree->nIntNodes; i++) 
        {
        p = tree->intDownPass[i];
        p->left->upDateTi=YES;
        p->right->upDateTi=YES;
        p->upDateCl=YES; 
        }

    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   Likelihood_Pairwise: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int Likelihood_Pairwise (int division, int chain, MrBFlt *lnL)
{
    /* calculate Jukes Cantor likelihood, using pw counts and tps. */
    
    int         i, j, idx, p, nijk;
    MrBFlt      like;
    CLFlt       *doubP, pijk;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    (*lnL)=0.0;

    for (p=0; p<numPairs; p++)
        {

        /* find transition probabilities */
        doubP = m->doubletProbs[m->pwIndex[chain][p]];

        idx=0;
        for (i=0; i<4; i++) 
            {
            for (j=0; j<4; j++) 
                {
                like = 0.0;
                nijk=pairwiseCounts[tIdx(p,i,j,4,4)];
                pijk=doubP[idx++];

                if (0) {
                MrBayesPrint("doubl prob: %f \n", pijk);
                MrBayesPrint("count : %d \n", nijk);
                }

                if (nijk == 0)
                    like+=0;
                else 
                    like+=nijk*log(pijk);         
 
            /* check against LIKE_EPSILON (values close to zero are problematic) */
                if (pijk < LIKEPW_EPSILON)
                    {
                    (*lnL) = MRBFLT_NEG_MAX;
                    abortMove = YES;
                    return ERROR;
                    }

                (*lnL)+=like;
                }
            }
        }

    return (NO_ERROR);
}


MrBFlt LogLikePairwise(int chain) 
{
    ModelInfo  *m;
    Tree       *tree;
    MrBFlt     chainLnLike;

    int d=0; /*  for now, only implemented for a single DNA partition */
    
    chainLnLike = 0.0;
    
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);

    if (m->upDateCijk == YES)
        {
        if (UpDateCijk(d, chain)== ERROR)
            {
            (m->lnLike[2*chain+state[chain]]) = MRBFLT_NEG_MAX; /* effectively abort the move */
            return (MRBFLT_NEG_MAX);
            }
        m->upDateAll = YES;
        }


    CalcPairwiseDists_ReverseDownpass(tree,d,chain);
    TiProbsPairwise_Gen(d,chain);
    DoubletProbs_Gen(d,chain);

    Likelihood_Pairwise(d,chain,&(m->lnLike[2*chain+state[chain]]));
    chainLnLike += m->lnLike[2*chain+state[chain]];
    return(chainLnLike);
}


int FreePairwise(void) 
{
    free(pairwiseCounts);
    return (NO_ERROR);
}









 
/*******************************************
 *                                         *
 *                                         *
 * Triplet based quasi-likelihood stuff    *
 *                                         *
 *                                         *
 *******************************************/
                                           
int FreeTriples(void) 
{
    int i;

    for (i=0; i<numTrips; i++)
        if (tripleCounts[i] != NULL)
            free(tripleCounts[i]);
        else 
            return (ERROR);

    if (tripleCounts != NULL)
        free(tripleCounts);
    else 
        return (ERROR);
    
    return (NO_ERROR);
}




int CountTriplets(void) {

    int             i,j,k,c,id1,id2,id3,tripIdx;
    int             *counts;
    
    /*  For now, only inplemented for a single partition 
     *  of DNA type data */
    /*  
    MrBayesPrint("Counting Pairwise\n"); 
    if (mp->dataType!=DNA && mp->dataType != RNA)
        { 
        MrBayesPrint("%s Can't count pairwise site-patterns for non DNA data  \n", spacer);
        return(ERROR);
        }
    */

    if (defMatrix == NO) 
        {
        MrBayesPrint("%s Matrix needs to be defined before counting triplets  \n", spacer);
        return(ERROR);
        }  
    
    if (numDefinedPartitions > 1)    
        { 
        MrBayesPrint("%s Pairwise count likelihood only implemented for a single partition  \n", spacer);
        return(ERROR);
        }

    if (memAllocs[ALLOC_TRIPLES] == YES) 
        {
        FreeTriples();
        tripleCounts=NULL;
        memAllocs[ALLOC_TRIPLES] = NO;
        }


    numTrips=(int)numTaxa*(numTaxa-1)*(numTaxa-2)/6;
    MrBayesPrint("NumTriplets: %d \n", numTrips);

    /* first allocate pairwise doublet counts:  */
    tripleCounts=(int**)SafeMalloc(numTrips * sizeof(int*));
    for (i=0; i<numTrips; i++)
            tripleCounts[i]=(int*)SafeMalloc(64*sizeof(int));

    if (!pairwiseCounts)
        {
        MrBayesPrint("%s Problem allocating pairwise counts! \n", spacer);
        free(pairwiseCounts);
        return(ERROR);
        }

    /* now count the doublets across taxa pairs */
    tripIdx=0;
    for (i=0; i<(numTaxa-2); i++)
        {
        for (j=i+1; j<(numTaxa-1); j++)
            {
            for (k=j+1; k<numTaxa; k++)
                {
                counts=tripleCounts[tripIdx++];
                for (c=0;c<numChar;c++)
                    {
                    if (matrix[pos(i,c,numChar)]==GAP | matrix[pos(j,c,numChar)]==GAP | matrix[pos(k,c,numChar)]==GAP)
                        continue;

                    /* nucleotides at position x of sequences i & j   */        
                    id1=toIdx(matrix[pos(i,c,numChar)]);
                    id2=toIdx(matrix[pos(j,c,numChar)]);
                    id3=toIdx(matrix[pos(k,c,numChar)]);

                    /* increment the count of that nucleotide pair, at pair k  */
                    counts[tIdx(id1,id2,id3,4,4)]++;
                    }
                }
            }
        }

    return (NO_ERROR);
}


int CalcTripletCnDists(int division, int chain)
{
    int         i,j,k,tripIdx;
    MrBFlt      *pwDists, *cnDists, pwd1, pwd2, pwd3; 
    ModelInfo   *m;  

    m=&modelSettings[division];
    pwDists=m->pwDists[chain];
    cnDists=m->tripleCnDists[chain];

    tripIdx=0;
    for (i=0; i<(numTaxa-2); i++)
        {
        for (j=i+1; j<(numTaxa-1); j++)
            {
            pwd1=pwDists[pairIdx(i,j,numTaxa)];
            for (k=j+1; k<numTaxa; k++)
                {
                pwd2=pwDists[pairIdx(i,k,numTaxa)];
                pwd3=pwDists[pairIdx(j,k,numTaxa)];

                cnDists[tripIdx++]=(pwd1+pwd2-pwd3)/(2.0);
                cnDists[tripIdx++]=(pwd1+pwd3-pwd2)/(2.0);
                cnDists[tripIdx++]=(pwd2+pwd3-pwd1)/(2.0);
                }
            }
        }

    
    return(NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   TiProbsTriplet_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbsTriplet_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    int         i, j, k, l, p, index, tripIdx;
    MrBFlt      *tripleDists; 
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange;
    CLFlt       *tiP;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    /* find pw dists and transition probabilities */
    tripleDists = m->tripleCnDists[chain];

    /* get base rate */
    baseRate = GetRate (division, chain);

    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
   
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    tripIdx=0;
    for (p=0; p<numTrips; p++)
        {

        for (l=0;l<3;l++)
            {
            
            length = tripleDists[tripIdx];
            tiP = m->tripleTiProbs[m->tripDistIndex[chain][tripIdx]];
            tripIdx++;

            index=0;
            for (k=0; k<m->numRateCats; k++)
                {

                t = length * baseRate * catRate[k];

                /* numerical errors will ensue if we allow very large or very small branch lengths,
                which might occur in relaxed clock models */
                if (t < TIME_MIN)
                    {
                    /* Fill in identity matrix */
                    for (i=0; i<4; i++)
                        {
                        for (j=0; j<4; j++)
                            {
                            if (i == j)
                                tiP[index++] = 1.0;
                            else
                                tiP[index++] = 0.0;
                            }
                        }
                    }
                else if (t > TIME_MAX)
                    {
                    /* Fill in stationary matrix */
                    for (i=0; i<4; i++)
                        for (j=0; j<4; j++)
                            tiP[index++] = 0.25;
                    }
                else
                    {
                    /* calculate probabilities */
                    pChange   = (CLFlt) (0.25 - 0.25 * exp(-(4.0/3.0)*t));
                    pNoChange = (CLFlt) (0.25 + 0.75 * exp(-(4.0/3.0)*t));

                    for (i=0; i<4; i++)
                        {
                        for (j=0; j<4; j++)
                            {
                            if (i == j)
                                tiP[index++] = pNoChange;
                            else
                                tiP[index++] = pChange;

                            }
                        }
                    }
                }
            } /*  end loop over pairs */
        }

    return (NO_ERROR);
}

/*-----------------------------------------------------------------
|
|   TiProbsTriplet_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TripletProbs_JukesCantor (int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    int         i, j, k, l, p, w, spIdx, tripLengthIdx, x,y,z, idx1, idx2, idx3;
    MrBFlt      *tripleDists; 
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange, numCats;
    CLFlt       *tiP1, *tiP2, *tiP3, *tripProb, *tripProbTemp;   
    ModelInfo   *m;
    int         indexStep, indexStart;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];
    numCats=((CLFlt)m->numRateCats);

    /*  tripProbTemp=(CLFlt*)SafeMalloc(4 * 4 * 4 * sizeof(CLFlt)); */

    indexStep=4*4;

    CLFlt dummy = 0.0; 
    tripLengthIdx=0;
    for (p=0; p<numTrips; p++)
        {
        indexStart=0;
        tripProb=m->tripleProbs[m->tripIndex[chain][p]];

        /*  initialize triplet probabilities to 0: */
        for (k=0; k<m->tripleProbsLength; k++)
                tripProb[k]=0.0;

        tiP1 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
        tiP2 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
        tiP3 = m->tripleTiProbs[m->tripDistIndex[chain][tripLengthIdx++]];
 
        for (k=0; k<m->numRateCats; k++)
            {

            spIdx=0;
            idx1=indexStart;   
            for (x=0; x<4; x++)
                {
                idx2=indexStart; 
                for (y=0; y<4; y++)
                    {
                    idx3=indexStart;
                    for (z=0; z<4; z++)
                        {
                        for (w=0; w<4; w++)
                            {
                                /*  internal triplet nucleotide is w */
                                /*  need transition probabilities w->x,w->y,w->z  */
                                tripProb[spIdx] += 0.25 * tiP1[idx1+w] * tiP2[idx2+w] * tiP3[idx3+w] / numCats;
                                /*
                                if (chain ==0 && x == 0 && y == 0 && z == 2 && p == 0) {

                                        MrBayesPrint("Site Pattern AAC|w=%d, triplet 1: ---- \n ", w);
                                        MrBayesPrint("  Triple Probs | rate cat = %d:\n ", k);
                                        MrBayesPrint("  tiP1[idx1+w]: %f,  tiP2[idx2+w]: %f tiP3[idr3+w]: %f  \n ", tiP1[idx1+w], tiP2[idx2+w], tiP3[idx3+w] );

                                        dummy += 0.25 * tiP1[idx1+w] * tiP2[idx2+w] * tiP3[idx3+w] / numCats ;
                                        if (w == 3) MrBayesPrint("  triplet prob for rate cat %d: %f \n\n", k, dummy);
                                }  */
                            }
                        spIdx++;
                        idx3+=4;
                        } /*  end of z loop */
                    idx2+=4;
                    } /*  end of y loop  */
                idx1+=4;
                } /*  end of x loop */

            indexStart+=indexStep;
            } /* end rate category loop  */

        /*   tripProbTemp[tripIdx]=;   */
        } /*  end loop over triplets */

    /*  free(tripProbTemp);  */
    return (NO_ERROR);
}

int Likelihood_Triples (int division, int chain, MrBFlt *lnL)
{
    /* calculate Jukes Cantor likelihood, using pw counts and tps. */
    
    int         i, j, k, idx, p, nijk;
    MrBFlt      like;
    CLFlt       *tripP, pijk;   
    ModelInfo   *m;

    /* MrBFlt  *bs;  don't need base freqs since this is JC submodel...*/
    m = &modelSettings[division];

    (*lnL) = 0.0;

    for (p=0; p<numTrips; p++)
        {

        /* find transition probabilities */
        tripP = m->tripleProbs[m->tripIndex[chain][p]];

        idx=0;
        for (i=0; i<4; i++) 
            {
            for (j=0; j<4; j++) 
                {
                for (k=0; k<4; k++)
                    {
                    like = 0.0;
                    nijk=tripleCounts[p][idx];
                    pijk=tripP[idx++];

                    if (nijk == 0)
                        like+=0; 
                    else 
                        like+=nijk*log(pijk);         
 
            /* check against LIKE_EPSILON (values close to zero are problematic) */
                    if (pijk < LIKEPW_EPSILON)
                        {
                        (*lnL) = MRBFLT_NEG_MAX;
                        abortMove = YES;
                        return ERROR;
                        }

                    *lnL+=like;
                    }
                }
            }
        }

    return (NO_ERROR);
}


MrBFlt LogLikeTriplet(int chain)
{
    ModelInfo  *m;
    Tree       *tree;
    MrBFlt     *lnL;
    int d=0; /*  for now, only implemented for a single DNA partition */

    lnL =  &(m->lnLike[2 * chain + state[chain]]);
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);

    CalcTripletCnDists(d, chain); 
    TiProbsTriplet_JukesCantor(d, chain); 
    TripletProbs_JukesCantor(d,chain);
    TIME(Likelihood_Triples(d,chain,lnL),CPULilklihood); 
}


MrBFlt LogLikeTriplet_Alpha(int chain)
{
    ModelInfo  *m;
    Tree       *tree;
    MrBFlt     lnL;
    int d=0; /*  for now, only implemented for a single DNA partition */

    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);

    CalcTripletCnDists(d, chain); 
    TiProbsTriplet_JukesCantor(d, chain); 
    TripletProbs_JukesCantor(d,chain);
    TIME(Likelihood_Triples(d,chain,&lnL),CPULilklihood); 
    return(lnL);
}







