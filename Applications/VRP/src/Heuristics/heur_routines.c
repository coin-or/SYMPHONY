#include <malloc.h>
#include <pvm3.h>

#include "vrp_const.h"
#include "pvm_error.h"
#include "heur_routines.h"

/*===========================================================================*/

/*------------------------------------------------------------------*\
| This function receives the data that the function broadcast sends  |
| and returns the tid of the sending process                         |
\*------------------------------------------------------------------*/

int receive(heur_prob *p)
{  
   int r_bufid, info, bytes, msgtag, parent;
   
   PVM_FUNC(r_bufid, pvm_recv(-1, VRP_DATA));
   PVM_FUNC(info, pvm_bufinfo(r_bufid, &bytes, &msgtag, &parent));
   PVM_FUNC(info, pvm_upkint(&(p->dist.wtype), 1, 1));
   PVM_FUNC(info, pvm_upkint(&(p->vertnum), 1, 1));
   PVM_FUNC(info, pvm_upkint(&(p->depot), 1, 1));
   PVM_FUNC(info, pvm_upkint(&p->capacity, 1, 1));
   p->demand = (int *) calloc (p->vertnum, sizeof(int));
   PVM_FUNC(info, pvm_upkint(p->demand, p->vertnum, 1));
   p->edgenum = p->vertnum*(p->vertnum-1)/2;
   if (p->dist.wtype){ /* not EXPLICIT */
      p->dist.coordx = (float *) calloc(p->vertnum, sizeof(float));
      p->dist.coordy = (float *) calloc(p->vertnum, sizeof(float));
      PVM_FUNC(info, pvm_upkfloat(p->dist.coordx, (int)p->vertnum, 1));
      PVM_FUNC(info, pvm_upkfloat(p->dist.coordy, (int)p->vertnum, 1));
      if ((p->dist.wtype == _EUC_3D) || (p->dist.wtype == _MAX_3D) || 
	  (p->dist.wtype == _MAN_3D)){
	 p->dist.coordz = (float *) calloc(p->vertnum, sizeof(float));
	 PVM_FUNC(info, pvm_upkfloat(p->dist.coordz, (int)p->vertnum, 1));
      }
   }
   else{ /* EXPLICIT */
      p->dist.cost = (int *) malloc (p->edgenum*sizeof(int));
      PVM_FUNC(info, pvm_upkint(p->dist.cost, (int)p->edgenum, 1));
   }
   PVM_FUNC(info, pvm_freebuf(r_bufid));
   return(parent);
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*\
| This function sends the final tour back to the parent                       |
\*---------------------------------------------------------------------------*/

void send_tour(_node *tour, int cost, int numroutes, int algorithm,
	       double cpu_time, int parent, int vertnum, int routes,
	       route_data *route_info)
{ 
   int s_bufid, info;
   
   PVM_FUNC(s_bufid, pvm_initsend(PvmDataRaw));
   PVM_FUNC(info, pvm_pkbyte((char *)tour, (vertnum)*sizeof(_node), 1));
   PVM_FUNC(info, pvm_pkint(&cost, 1, 1));
   PVM_FUNC(info, pvm_pkint(&numroutes, 1, 1));
   PVM_FUNC(info, pvm_pkint(&algorithm, 1, 1));
   PVM_FUNC(info, pvm_pkdouble(&cpu_time, 1, 1));
   if (routes)
      PVM_FUNC(info, pvm_pkbyte((char *)route_info, 
				(numroutes+1)*sizeof(route_data), 1));
   PVM_FUNC(info, pvm_send(parent, HEUR_TOUR));
   PVM_FUNC(info, pvm_freebuf(s_bufid));
   
   return;
}

/*===========================================================================*/

void free_heur_prob(heur_prob *p)
{
   if (p){
      if (p->cur_tour){
	 FREE(p->cur_tour->tour);
	 FREE(p->cur_tour->route_info);
	 FREE(p->cur_tour);
      }
      FREE(p->demand);
      FREE(p->dist.coordx);
      FREE(p->dist.coordy);
      FREE(p->dist.coordz);
      FREE(p->dist.cost);
      FREE(p);
   }
}

/*===========================================================================*/

void free_lb_prob(lb_prob *p)
{
   if (p){
      FREE(p->demand);
      FREE(p->dist.coordx);
      FREE(p->dist.coordy);
      FREE(p->dist.coordz);
      FREE(p->dist.cost);
      FREE(p);
   }
}

