/* stack.c
 *
 * implements a very simple stack to flood fill the mesh */
#include <stdio.h>
#include <stdlib.h>

#include "stack.h"

struct point_data *last=NULL;

struct point_data *stack_push(struct point_data *stack, long int i, long int j, long int k) { 
  
  struct point_data *p;

  p = malloc(sizeof(struct point_data));
  p->i = i;
  p->j = j;
  p->k = k;

  p->first = stack;
  p->next  = NULL;

  if(last == NULL) {
    stack = p;
    stack->first = stack;
    last  = stack;
    p->prev = NULL;
  }
  else {
    p->prev = last;
    last->next = p;
    last = p;
  }
  
  return stack;
}

struct point_data *stack_pop(struct point_data *stack)
{
  struct point_data *p;

  if(last == NULL) return NULL;

  p = last;

  last = last->prev;
  
  if(last != NULL)
    last->next = NULL;

  return p;
}
