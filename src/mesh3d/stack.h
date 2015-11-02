/* stack.h
 *
 * header to stack.c  
 * implements the point struct, used for the stack */

#ifndef STACK_H
#define STACK_H

struct point_data {
  long int i, j, k;
  struct point_data *next;
  struct point_data *prev;
  struct point_data *first;
};

struct point_data *stack_push(struct point_data *stack, long int i, long int j, long int k);

struct point_data *stack_pop(struct point_data *stack);

#endif
