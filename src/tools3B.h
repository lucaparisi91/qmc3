#ifndef TOOLS3B_H
#define TOOLS3B_H

#define LOOP3B( NPARTICLES,  BODY )	\
{ \
  int kj=0;  \
  int ki2=0; \
  int ki = 0; \
  for(int k=0;k<NPARTICLES;k++) \
    { \
      int ji=0;\
      for (int j=0;j<k;j++ & kj++) \
	{ \
	 ki=ki2; \
	  for(int i=0;i<j;i++ & ji++ & ki++) \
	    { \
              BODY \
	    }  \
	}\
      ki2+=k; \
      }	      \
} \

#define LOOP3B_DIS( N1, N2 , N3,  BODY )		\
{ \
  int kj=0;  \
  int ki2=0; \
  int ki = 0; \
  for(int k=0;k<N3;k++) \
    { \
      int ji=0;\
      for (int j=0;j<N2;j++ & kj++) \
	{ \
	 ki=ki2; \
	  for(int i=0;i<N1;i++ & ji++ & ki++) \
	    { \
              BODY \
	    }  \
	}\
      ki2+=N1; \
      }	      \
} \


#endif
