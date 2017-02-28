#ifndef KNOTANALYSIS_H_INCLUDED
#define KNOTANALYSIS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
	double * chain;
	int size;
	int state;
	int state1;
	int state2;
	double start;
	double start1;
	double start2;
	double end;
	double end1;
	double end2;
	double length;
	double length1;
	double length2;
	double position;
}jKN;

jKN* jKN_alloc( double * chain,const int size)
{
	jKN* knot=(jKN*)malloc(sizeof(jKN));
	knot->chain=chain;
	knot->size=size;
	knot->state=-1;
	knot->state1=-1;
	knot->state2=-1;
	knot->start=-1.0;
	knot->start1=-1.0;
	knot->start2=-1.0;
	knot->end=-1.0;
	knot->end1=-1.0;
	knot->end2=-1.0;
	knot->length=-1.0;
	knot->length1=-1.0;
	knot->length2=-1.0;
	knot->position=-1.0;
	return knot;
}

#ifndef J_VECTORMATRIX_INCLUDED
#define J_VECTORMATRIX_INCLUDED

int** matrix_int_calloc(const int dimension1,const int dimension2)
{
	int i,j;
	int **Matrix;
	Matrix=(int**)malloc(dimension1*sizeof(int*));
	for(i=0; i<dimension1; i++) {
		Matrix[i]=(int*)malloc(dimension2*sizeof(int));
		for(j=0; j<dimension2; j++) {
			Matrix[i][j]=0;
		}
	}
	return Matrix;
}

void matrix_int_free(int **Matrix,const int dimension1)
{
	int i;
	for(i=0; i<dimension1; i++) {
		free(Matrix[i]);
	}
	free(Matrix);
}

double** matrix_double_calloc(const int dimension1,const int dimension2)
{
	int i,j;
	double **Matrix;
	Matrix=(double**)malloc(dimension1*sizeof(double*));
	for(i=0; i<dimension1; i++) {
		Matrix[i]=(double*)malloc(dimension2*sizeof(double));
		for(j=0; j<dimension2; j++) {
			Matrix[i][j]=0.0;
		}
	}
	return Matrix;
}

void matrix_double_free(double **Matrix,const int dimension1)
{
	int i;
	for(i=0; i<dimension1; i++) {
		free(Matrix[i]);
	}
	free(Matrix);
}

double* vector_double_alloc(const int dimension)
{
	double * Vector;
	Vector=(double*)malloc(dimension*sizeof(double));
	return Vector;
}

void vector_double_free(double*Vector)
{
	free(Vector);
}

int * vector_int_alloc(const int dimension)
{
	int * Vector;
	Vector=(int*)malloc(dimension*sizeof(int));
	return Vector;
}

void vector_int_free(int * Vector)
{
	free(Vector);
}

#endif

int KnotAnalyse(const double * chain,int size,int projection,double *tvalue,double *tstart, double *tend)
{
	int ni,nj,nk,Nbnd,dim1,dim2,dim3,crosscnt=0,si0j0_j1,si0j1_j0,sj0j1_i0,sj0j1_i1,si0j0_i1,si0j1_i1,maxcrss,R3cnt;
	int **crssngs,**crssndx,R1,R2,R3,*p1,*p2,*p3,*s1,*s2,*s3,swap,p1o=0,p2o=0,p3o=0,s1o=0,s2o=0,s3o=0,stuckcnt=0;
	double xi1,xi0,yi1,yi0,xj1,xj0,yj1,yj0,mi0j0,mi0j1,mj0j1,mi0i1,*tvals,xcr,ycr,icrl,jcrl,zi0,zi1,zj0,zj1,zcri,zcrj;
	double cx,cy,cz,Rxz,Ryz,cos,sin,*chainR;
	Nbnd=size-1; // number of bonds, = number of points -1
	maxcrss=(int)(((Nbnd/2)-1)*(Nbnd-1)); // maximum possible crossings
	crssngs=matrix_int_calloc(maxcrss,3); // storage for crossings as they are found
	p1=vector_int_alloc(Nbnd);
	p2=vector_int_alloc(Nbnd);
	p3=vector_int_alloc(Nbnd);
	s1=vector_int_alloc(Nbnd);
	s2=vector_int_alloc(Nbnd);
	s3=vector_int_alloc(Nbnd);
	tvals  =vector_double_alloc(maxcrss); // stores 'length along chain' values
	chainR =vector_double_alloc(3*size);
	/// RE-ORIENTATE CHAIN ALONG END TO END AXIS ///
	cx=chain[3*size-3]-chain[0]; // difference in x-coord
	cz=chain[3*size-1]-chain[2]; // difference in z-coord
	Rxz=sqrt(cx*cx+cz*cz);
	cos=cz/Rxz;
	sin=cx/Rxz;
	for(ni=0;ni<size;ni++) { // rotate original chain around y-axis
		chainR[3*ni+0]=cos*chain[3*ni+0]-sin*chain[3*ni+2];
		chainR[3*ni+1]=chain[3*ni+1]; // leave y-coord
		chainR[3*ni+2]=sin*chain[3*ni+0]+cos*chain[3*ni+2];
	}
	cy=chainR[3*size-2]-chainR[1]; // new difference in y-coord
	cz=chainR[3*size-1]-chainR[2]; // new difference in z-coord
	Ryz=sqrt(cy*cy+cz*cz);
	cos=cz/Ryz;
	sin=cy/Ryz;
	for(ni=0;ni<size;ni++) { // further rotation around x-axis
		chainR[3*ni+1]=cos*chainR[3*ni+1]-sin*chainR[3*ni+2];
		chainR[3*ni+2]=sin*chainR[3*ni+1]+cos*chainR[3*ni+2];
	}
	dim3=projection; // should definitely not be z-axis (=2)
	dim1=(dim3+1)%3;
	dim2=(dim3+2)%3;
	/// FIND ALL CROSSINGS ///
	for(ni=0;ni<Nbnd-2;ni++) {
		xi0=chainR[3*ni+dim1];
		yi0=chainR[3*ni+dim2];
		zi0=chainR[3*ni+dim3];
		xi1=chainR[3*ni+3+dim1];
		yi1=chainR[3*ni+3+dim2];
		zi1=chainR[3*ni+3+dim3];
		for(nj=ni+2;nj<Nbnd;nj++) {
			xj0=chainR[3*nj+dim1];
			yj0=chainR[3*nj+dim2];
			zj0=chainR[3*nj+dim3];
			xj1=chainR[3*nj+3+dim1];
			yj1=chainR[3*nj+3+dim2];
			zj1=chainR[3*nj+3+dim3];
			mj0j1=(yj0-yj1)/(xj0-xj1); // gradient of j0->j1 line
			sj0j1_i0=(((mj0j1)*(xi0-xj0)-(yi0-yj0))>0); // which side is i0 compared to j0->j1 line ?
			sj0j1_i1=(((mj0j1)*(xi1-xj0)-(yi1-yj0))>0); // which side is i1 compared to j0->j1 line ?
			if(sj0j1_i0==sj0j1_i1) continue; // they need to be opposite sides
			mi0j0=(yi0-yj0)/(xi0-xj0);
			si0j0_j1=(((mi0j0)*(xj1-xi0)-(yj1-yi0))>0); // which side is j1 compared to i0->j0 line ?
			si0j0_i1=(((mi0j0)*(xi1-xi0)-(yi1-yi0))>0); // which side is i1 compared to i0->j0 line ?
			if(si0j0_i1!=si0j0_j1) continue; // they need to be on the same side
			mi0j1=(yi0-yj1)/(xi0-xj1);
			si0j1_j0=(((mi0j1)*(xj0-xi0)-(yj0-yi0))>0);	// which side is j0 compared to i0->j1 line ?
			si0j1_i1=(((mi0j1)*(xi1-xi0)-(yi1-yi0))>0);	// which side is i1 compared to i0->j1 line ?
			if(si0j1_i1!=si0j1_j0) continue; // they need to be on the same side
			mi0i1=(yi0-yi1)/(xi0-xi1);
			xcr=((mi0i1*xi0-yi0)-(mj0j1*xj0-yj0))/(mi0i1-mj0j1); // x-coord of bond cross
			ycr=mi0i1*(xcr-xi0)+yi0; // y coord of bond cross
			icrl=sqrt(((xcr-xi0)*(xcr-xi0)+(ycr-yi0)*(ycr-yi0))/((xi1-xi0)*(xi1-xi0)+(yi1-yi0)*(yi1-yi0))); // fractional distance to crossing points from ith bead
			jcrl=sqrt(((xcr-xj0)*(xcr-xj0)+(ycr-yj0)*(ycr-yj0))/((xj1-xj0)*(xj1-xj0)+(yj1-yj0)*(yj1-yj0))); // fractional distance to crossing points from jth bead
			zcri=zi0+(zi1-zi0)*icrl; // zi value at the crossing point
			zcrj=zj0+(zj1-zj0)*jcrl; // zj value at the crossing point
			crssngs[crosscnt][2]=(zcri>zcrj?1:-1); // 1 if first bond above second, -1 if first bond below
			tvals[2*crosscnt  ]=ni+icrl; // record position along chain of primary bond cross
			tvals[2*crosscnt+1]=nj+jcrl; // record position along chain of secondary bond cross
			crosscnt++;
		}
	}
	crssndx=matrix_int_calloc(2*(dim1=crosscnt)+2,3); // cross-reference, ordered storage for crossings after sorting, use re-use dim1 as original crosscnt to free crssndx later
//	printf("crosscnt: %d\n",crosscnt);
	if(crosscnt==0) goto knotclearup;
	if(crosscnt>maxcrss) {
		printf("Warning: Too many crossings, increase buffer size.");
		crosscnt=-2;
		goto knotclearup;
	}
	/// RANK BOND CROSSINGS ///
//	printf("Checkpoint: A\n");

	for(ni=0;ni<crosscnt;ni++) {
		crssngs[ni][0]+=1; // index starts at 1
		crssngs[ni][1]+=2; // by def, secondary bond is further along than primary, no need to compare
		for(nj=ni+1;nj<crosscnt;nj++) { // step through all other t-values and increment ranks
			if(tvals[2*ni]  >tvals[2*nj])   crssngs[ni][0]++; // compare position of primary bond to all others
			else crssngs[nj][0]++;
			if(tvals[2*ni]  >tvals[2*nj+1]) crssngs[ni][0]++;
			else crssngs[nj][1]++;
			if(tvals[2*ni+1]>tvals[2*nj])   crssngs[ni][1]++; // compare position of secondary bond to all others
			else crssngs[nj][0]++;
			if(tvals[2*ni+1]>tvals[2*nj+1]) crssngs[ni][1]++;
			else crssngs[nj][1]++;
		}
	}
	/// MAKE KNOT MATRIX AND INDEX CROSSINGS///
//	printf("Checkpoint: B\n");
	for(ni=0;ni<crosscnt;ni++) {
		crssndx[crssngs[ni][0]][0]= crssngs[ni][1]; // which bond does primary bond cross?
		crssndx[crssngs[ni][0]][1]= crssngs[ni][2]; // does it go under or over?
		crssndx[crssngs[ni][0]][2]= crssngs[ni][0]; // record original position
		crssndx[crssngs[ni][1]][0]= crssngs[ni][0]; // same for secondary bond
		crssndx[crssngs[ni][1]][1]=-crssngs[ni][2]; //
		crssndx[crssngs[ni][1]][2]= crssngs[ni][1]; //
	}
	maxcrss=crosscnt; // reassign maxcrss as total number of intial crossings
	while(crosscnt!=0){
		R1=R2=R3=0;
//		printf("Checkpoint: C\n");
		/// REIDEMEISTER 1 ///
		for(ni=1;ni<2*crosscnt+1;ni++) { // step through all crossings
			if(crssndx[ni][0]==ni+1) { // identify twists
				for(nj=1;nj<=2*crosscnt-2;nj++) { // remove that crossing
					if(nj>ni-1) { // slide primary crossings 'up' if beyond removed setion
						crssndx[nj][0]=crssndx[nj+2][0]; // copy secondary rank down
						crssndx[nj][1]=crssndx[nj+2][1]; // copy under/over marker down
						crssndx[nj][2]=crssndx[nj+2][2]; // copy the original cross rank down
					}
					if(crssndx[nj][0]>ni+1) crssndx[nj][0]-=2; // reduce rank of secondary crossings after removed section by two
				}
				crosscnt--; // decrement number of crossings in current 'diagram'
				stuckcnt=0; // reset the stuck counter if the knot matrix changes size
				R1=1; // flag move performed
				break;
			}
		}
		if(R1==1) continue;
//		printf("Checkpoint: D\n");
		/// REIDEMEISTER 2 ///
		for(ni=1;ni<2*crosscnt;ni++) { // step through all crossings
			if(abs(crssndx[ni+1][0]-crssndx[ni][0])==1 && crssndx[ni][1]==crssndx[ni+1][1]) { // identify adjacent crossings both under or both over
				nk=(crssndx[ni+1][0]<crssndx[ni][0]?crssndx[ni+1][0]:crssndx[ni][0]); // lowest secondary rank (always>ni)
				for(nj=1;nj<=2*crosscnt-4;nj++) { // step through all crossings again
					if     (nj>nk-3) { // if rank higher than both pairs removed
						crssndx[nj][0]=crssndx[nj+4][0]; // slide primary crossings four 'up'
						crssndx[nj][1]=crssndx[nj+4][1];
						crssndx[nj][2]=crssndx[nj+4][2];
					}
					else if(nj>ni-1) { // rank only higher than one pair
						crssndx[nj][0]=crssndx[nj+2][0]; // slide primary crossings two 'up'
						crssndx[nj][1]=crssndx[nj+2][1];
						crssndx[nj][2]=crssndx[nj+2][2];
					}
					if     (crssndx[nj][0]>nk+1) crssndx[nj][0]-=4; // lower secondary crossings ranks
					else if(crssndx[nj][0]>ni+1) crssndx[nj][0]-=2;
				}
				crosscnt-=2; // double decrement number of diagram crossings
				stuckcnt=0;
				R2=1; // flag move performed
				break;
			}
		}
		if(R2==1) continue;
		R3cnt=0;
//		printf("Checkpoint: E\n");
		/// REIDEMEISTER 3 ///
		for(ni=1;ni<2*crosscnt;ni++) { // step through all but last crossings and find available R3 moves
//			printf("Checkpoint: Ea\n");
			if(crssndx[ni][1]==crssndx[ni+1][1]) { // look for arcs with two adjacent crossings on same side
//				printf("Checkpoint: Eb\n");
				for(nj=-1;nj<2;nj+=2) {
//					printf("Checkpoint: Ec\n");
					for(nk=-1;nk<2;nk+=2) {
//						printf("Checkpoint: Ed\n");
						if(crssndx[ni][0]==2*crosscnt&&nj==1) continue; // dont look outside of matrix!
//						printf("Checkpoint: Ed1\n");
//						printf("ni nj nk: %d %d %d\n",ni,nk,nk);
//						printf("crssndx: %d %d\n",crssndx[ni][0],crssndx[ni+1][0]);
						if(crssndx[crssndx[ni][0]+nj][0]==crssndx[ni+1][0]+nk) { // do they share a neightbour?
//							printf("Checkpoint: Ee\n");
							p1[R3cnt]=ni;
							s1[R3cnt]=crssndx[ni][0];
							p2[R3cnt]=ni+1;
							s2[R3cnt]=crssndx[ni+1][0];
							p3[R3cnt]=s1[R3cnt]+nj;
							s3[R3cnt]=s2[R3cnt]+nk;
							if(p1[R3cnt]==p1o&&s1[R3cnt]==s1o&&p2[R3cnt]==p2o&&s2[R3cnt]==s2o&&p3[R3cnt]==p3o&&s3[R3cnt]==s3o) continue; // dont 'undo' last R3 move
							R3cnt++;
						}
					}
				}
			}
		}
//		printf("Checkpoint: E1\n");
		if(R3cnt>0) {
			for(ni=0;ni<R3cnt;ni++) { // look for good moves
				if(abs(s1[ni]-s2[ni])==1||abs(s3[ni]-p1[ni])==1||abs(p3[ni]-p2[ni])==1) { // will an R1 move become available?
					R3cnt=ni; // if so perform corresponding R3 move
					R3=1; // flag R3 move found
					break;
				}
				if(crssndx[2*s3[ni]-s2[ni]][0]==2*p1[ni]-p2[ni]&&crssndx[2*s3[ni]-s2[ni]][1]==crssndx[s2[ni]][1]) { // will an R2 move become available
					R3cnt=ni;
					R3=1;
					break;
				}
				if(crssndx[2*p3[ni]-s1[ni]][0]==2*p2[ni]-p1[ni]&&crssndx[2*p3[ni]-s1[ni]][1]==crssndx[s1[ni]][1]) {
					R3cnt=ni;
					R3=1;
					break;
				}
				if(crssndx[2*s2[ni]-s3[ni]][0]==2*s1[ni]-p3[ni]&&crssndx[2*s2[ni]-s3[ni]][1]==crssndx[s3[ni]][1]) {
					R3cnt=ni;
					R3=1;
					break;
				}
			}
//			printf("Checkpoint: E2\n");
			if(R3==0) {
				R3cnt=rand()%R3cnt; // choose random move if no 'good' ones
			}
			p1o=crssndx[p1[R3cnt]][0]=s3[R3cnt]; // perform chosen R3 move
			p2o=crssndx[p2[R3cnt]][0]=p3[R3cnt];
			p3o=crssndx[p3[R3cnt]][0]=p2[R3cnt];
			s1o=crssndx[s1[R3cnt]][0]=s2[R3cnt];
			s2o=crssndx[s2[R3cnt]][0]=s1[R3cnt];
			s3o=crssndx[s3[R3cnt]][0]=p1[R3cnt];
			swap=crssndx[p2[R3cnt]][1]; // map correct under/over values
			crssndx[p2[R3cnt]][1]=crssndx[p1[R3cnt]][1];
			crssndx[p1[R3cnt]][1]=swap;
			swap=crssndx[s1[R3cnt]][1];
			crssndx[s1[R3cnt]][1]=crssndx[p3[R3cnt]][1];
			crssndx[p3[R3cnt]][1]=swap;
			swap=crssndx[s3[R3cnt]][1];
			crssndx[s3[R3cnt]][1]=crssndx[s2[R3cnt]][1];
			crssndx[s2[R3cnt]][1]=swap;
			stuckcnt++; // count how many consecutive R3 moves being performed
		}
		if(R1==0&&R2==0&&R3==0) break; // no moves can be identified
		if(stuckcnt>50) { // if only performing R3 moves, (presumed cyclically), break
			printf("\nKnot Warning: Stuck performing R3 moves\n");
			break;
		}
	}
//	printf("Checkpoint: F\n");
	if(crosscnt>0) { // check orignal position of first bond in reduced scheme
		for(ni=0;ni<maxcrss;ni++) { // search for entry in list of crossings ordered via first encountered, ie compare ranks
			if(crssndx[1][2]         ==crssngs[ni][0]) *tstart=tvals[2*ni];
			if(crssndx[1][2]         ==crssngs[ni][1]) *tstart=tvals[2*ni+1];
			if(crssndx[2*crosscnt][2]==crssngs[ni][0]) *tend  =tvals[2*ni];
			if(crssndx[2*crosscnt][2]==crssngs[ni][1]) *tend  =tvals[2*ni+1];
		}
		*tvalue=(*tend-*tstart)/Nbnd; // relative length of chain involved in actual knot
	}
	knotclearup:
//	printf("Checkpoint: G\n");
	matrix_int_free(crssngs,maxcrss);
	matrix_int_free(crssndx,2*dim1+1);
	vector_double_free(tvals);
	vector_double_free(chainR);
	return crosscnt;
}

void KnotScan(jKN *knot)
{
//	printf("Pt:1\n");
	knot->state1=KnotAnalyse(knot->chain,knot->size,1,&(knot->length1),&(knot->start1),&(knot->end1)); // analyse knot from two projections
//	printf("Pt:2\n");
	knot->state2=KnotAnalyse(knot->chain,knot->size,0,&(knot->length2),&(knot->start2),&(knot->end2));
//	printf("Pt:3\n");
	knot->state=(knot->state1<=knot->state2?knot->state1:knot->state2); // take lowest knot state as algorithm usually overestimates
//	printf("Pt:4\n");
	if(knot->state>0) {
		knot->length=0.5*(knot->length1+knot->length2); // take average length
		knot->start=0.5*(knot->start1+knot->start2); // take average start
		knot->end=0.5*(knot->end1+knot->end2); // take average end
		knot->position=0.5*(knot->start+knot->end); // find average position
	} else {
		knot->start=-1.0;
		knot->start1=-1.0;
		knot->start2=-1.0;
		knot->end=-1.0;
		knot->end1=-1.0;
		knot->end2=-1.0;
		knot->length=-1.0;
		knot->length1=-1.0;
		knot->length2=-1.0;
		knot->position=-1.0;
	}
}

#endif
