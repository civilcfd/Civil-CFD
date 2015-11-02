/* vector_macros.h */

#ifndef _VECTOR_MACROS_H
#define _VECTOR_MACROS_H

#define vector_add(a,b,c) \
  (a)[0] = (b)[0] + (c)[0]; \
  (a)[1] = (b)[1] + (c)[1]; \
  (a)[2] = (b)[2] + (c)[2];

#define vector_subtract(a,b,c) \
  (a)[0] = (b)[0] - (c)[0]; \
  (a)[1] = (b)[1] - (c)[1]; \
  (a)[2] = (b)[2] - (c)[2];

#define vector_multiply(a,b,c) \
  (a)[0] = (b)[0] * c; \
  (a)[1] = (b)[1] * c; \
  (a)[2] = (b)[2] * c;

#define cross_product(a,b,c) \
  (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
  (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
  (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define inner_product(v,q) \
    ((v)[0] * (q)[0] + \
    (v)[1] * (q)[1] + \
    (v)[2] * (q)[2])
    
#define vector_magnitude(v) \
    (sqrt(pow((v)[0],2) + pow((v)[1],2) + pow((v)[2],2)))

#endif

