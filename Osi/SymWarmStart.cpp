
#include "SymWarmStart.hpp"

//#############################################################################

SymWarmStart::SymWarmStart(warm_start_desc *& ws, 
			   bool dominateWarmStart)
{
   int i, num;
   cut_data *temp;

   if(!ws){
      cout<<"SymWarmStart(): Given warmStart is empty!"<<endl;
   }
   else{
      warmStart_ = new warm_start_desc;     

      if(!dominateWarmStart){
	 memcpy(warmStart_, ws, sizeof(warm_start_desc));
	 num = warmStart_->cut_num;
	 warmStart_->cuts = new cut_data* [num];
	 for(i = 0; i<num; i++){
	    temp = new cut_data;
	    memcpy(temp, ws->cuts[i], sizeof(cut_data));
	    temp->coef = new char[temp->size];
	    memcpy(temp->coef, ws->cuts[i]->coef, sizeof(char)*temp->size);
	    warmStart_->cuts[i] = temp;
	 }
	 warmStart_->rootnode = new bc_node;	 
	 copyTree_(warmStart_->rootnode, ws->rootnode);
      }
      else{
	 warmStart_ = ws;
	 ws = 0;
      }  
   }
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(char * f)
{

   ifstream is(f);

   if(!is){
      cout<<"SymWamrStart(): Can not open the given file!"<<endl;
   }
   else{
      warmStart_ = new warm_start_desc;            
      readWarmStartFromFile_(is);
   }

   is.close();
}

/*===========================================================================*/
/*===========================================================================*/

SymWarmStart::SymWarmStart(const SymWarmStart & symWS)
{
   int i, num;
   warm_start_desc * wsCopy;
   cut_data * temp;
   SymWarmStart * sWS = const_cast<SymWarmStart *>(&symWS);

   if(&symWS){
      wsCopy = const_cast<warm_start_desc *>(sWS->getWarmStartDesc());
      warmStart_ = new warm_start_desc;      
      memcpy(warmStart_, wsCopy, sizeof(warm_start_desc));
      num = warmStart_->cut_num;
      warmStart_->cuts = new cut_data* [num];
	 for(i = 0; i<num; i++){
	    temp = new cut_data;
	    memcpy(temp, wsCopy->cuts[i], sizeof(cut_data));
	    temp->coef = new char[temp->size];
	    memcpy(temp->coef, wsCopy->cuts[i]->coef, 
		   sizeof(int)*temp->size);	   
	    warmStart_->cuts[i] = temp;
	 }

	 warmStart_->rootnode = new bc_node;	 
	 copyTree_(warmStart_->rootnode, wsCopy->rootnode);
   }
}

/*===========================================================================*/
/*===========================================================================*/
SymWarmStart::~SymWarmStart()
{

   int i, temp;
   warm_start_desc * ws = warmStart_;
   if(ws) {
      if(ws->rootnode) {
	 trimTree(ws->rootnode);
      }
      if(ws->cuts){
	 temp = ws->cut_num;
	 for(i = 0; i < temp; i++){
	    if(ws->cuts[i]){
	       if(ws->cuts[i]->coef){
		  delete [] ws->cuts[i]->coef;
	       }
	       delete ws->cuts[i];
	    }
	 }
	 delete [] ws->cuts;
      }
      delete ws;
   }
   warmStart_ = 0;
}

/*===========================================================================*/
/*===========================================================================*/

CoinWarmStart * SymWarmStart::clone () const
{
   return new SymWarmStart(*this);
}

/*===========================================================================*/
/*===========================================================================*/

const warm_start_desc * SymWarmStart::getWarmStartDesc() const
{
   if(warmStart_){
      return warmStart_;
   }
   else{
      cout<<"getWarmStart(): No loaded warm start desc. to return!"<<endl;
      return 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::trimTree(bc_node *& node)
{
   int i;
   
   if(node){
      for(i = 0; i<node->bobj.child_num; i++){
	 trimTree(node->children[i]);
      }
      deleteNode_(node);        
      node = 0;
   }
}

/*===========================================================================*/
/*===========================================================================*/

bool SymWarmStart::writeToFile(char * f)
{
   
   int i, j, temp;
   if(!warmStart_){
      cout<<"There is no loaded warmStart to write!"<<endl;
      return false;
   }
   else{
      ofstream os(f);
      os.precision(6);
      os.fill('#');
      os.width(75);
      os<<""<<endl<<" BOUND INFO "<<endl;
      os.width(75);
      os<<""<<endl;
      os<<" PHASE     : "<<warmStart_->phase<<endl;
      os<<" LB        : "<<warmStart_->lb<<endl;
      os<<" HAS_UB    : "<<(int)warmStart_->has_ub<<endl;
      os<<" UB        : "<<warmStart_->ub<<endl;

      os<<endl;
      os.width(75);
      os<<""<<endl<<" CUT INFO "<<endl;
      os.width(75);
      os<<""<<endl;
      os<<" CUT_NUM           : "<<warmStart_->cut_num<<endl;
      os<<" ALLOCATED_CUT_NUM : "<<warmStart_->allocated_cut_num<<endl;

      os.fill(' ');
      os<<endl;
      //FIX_ME! WHAT TYPE A CUT CAN BE OTHER THAN EXPLICIT_ROW

      cut_data **cuts = warmStart_->cuts;
      
      for(i=0; i<warmStart_->cut_num; i++){
	 os<<" CUT "<<i<<":"<<endl;
	 os.width(5);
	 os<<" SIZE          : "<<cuts[i]->size<<endl;
	 temp = cuts[i]->coef[0];
	 os<<" NUM_ELEMENTS  : "<<temp<<endl;
	 os<<" INDICES       : ";
	 for(j=1; j<=temp; j++){
	    if(j%20==0){
	       os<<endl;
	       os.width(18);
	    }
	    os<<(int)cuts[i]->coef[j]<<" ";
	 }
	 os<<endl;
	 os<<" VALUES        : ";
	 for(j=1; j<=temp; j++){
	    if(j%28==0){
	       os<<endl;
	       os.width(18);
	    }
	    os<<(double)cuts[i]->coef[temp+j]<<" ";
	 }
	 os<<endl;
	 os<<" RHS           : "<<cuts[i]->rhs<<endl;
	 os<<" RANGE         : "<<cuts[i]->range<<endl;
	 os<<" TYPE          : "<<(int)cuts[i]->type<<endl;
	 os<<" SENSE         : "<<cuts[i]->sense<<endl;
	 os<<" DELETABLE     : "<<(int)cuts[i]->deletable<<endl;
	 os<<" BRANCH        : "<<(int)cuts[i]->branch<<endl;
	 os<<" NAME          : "<<cuts[i]->name<<endl;

	 os<<endl;
      }

      os<<endl;
      os.fill('#');
      os.width(75);
      os<<""<<endl<<" PROBLEM STATISTICS"<<endl;
      os.width(75);
      os<<""<<endl;      

      problem_stat stat= warmStart_->stat;


      os<<" ROOT LB                : "<<stat.root_lb<<endl;
      os<<" CUTS IN POOL           : "<<stat.cuts_in_pool<<endl;
      os<<" MAXIMIM DEPTH          : "<<stat.max_depth<<endl;
      os<<" DIVING CHAINS          : "<<stat.chains<<endl;
      os<<" DIVING STOPS           : "<<stat.diving_halts<<endl;
      os<<" TREE SIZE              : "<<stat.tree_size<<endl;
      os<<" CREATED NODES          : "<<stat.created<<endl;
      os<<" ANALYZED NODES         : "<<stat.analyzed<<endl;
      os<<" LEAVES BEFORE TRIMMING : "<<stat.leaves_before_trimming<<endl;
      os<<" LEAVES BEFORE TRIMMING : "<<stat.leaves_after_trimming<<endl;
      os<<" NF STATUS OF ROOT      : "<<(int)stat.nf_status<<endl;
      
      os<<endl;

      os.width(75);
      os<<""<<endl<<" COMPUTATION TIMES"<<endl;
      os.width(75);
      os<<""<<endl;      
      
      node_times compT = warmStart_->comp_times;

      os<<" COMMUNICATION       : "<<compT.communication<<endl;
      os<<" LP                  : "<<compT.lp<<endl;
      os<<" SEPARATION          : "<<compT.separation<<endl;
      os<<" FIXING              : "<<compT.fixing<<endl;
      os<<" PRICING             : "<<compT.pricing<<endl;
      os<<" STRONG BRANCHING    : "<<compT.strong_branching<<endl;
      os<<" WALL CLOCK_LP       : "<<compT.wall_clock_lp<<endl;
      os<<" RAMP UP_TM          : "<<compT.ramp_up_tm<<endl;
      os<<" RAMP UP_LP          : "<<compT.ramp_up_lp<<endl;
      os<<" RAMP DOWN TIME      : "<<compT.ramp_down_time<<endl;
      os<<" IDLE DIVING         : "<<compT.idle_diving<<endl;
      os<<" IDLE NODE           : "<<compT.idle_node<<endl;
      os<<" IDLE NAMES          : "<<compT.idle_names<<endl;
      os<<" IDLE CUTS           : "<<compT.idle_cuts<<endl;
      os<<" START NODE          : "<<compT.start_node<<endl;
      os<<" CUT POOL            : "<<compT.cut_pool<<endl;

      os<<endl;

      os.width(75);
      os<<""<<endl<<" TREE DESCRIPTION"<<endl;
      os.width(75);
      os<<""<<endl;      

      writeTreeToFile_(warmStart_->rootnode, os);
      os.close();
      return true;
   }   
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::writeTreeToFile_(bc_node *root, ofstream &os)
{
   int i;
   writeNode_(root, os);
   
   for(i=0; i<root->bobj.child_num; i++){
      writeTreeToFile_(root->children[i], os);
   }
}

/*===========================================================================*/
/*===========================================================================*/
void SymWarmStart::writeNode_(bc_node *node, ofstream &os)
{
   int i;
   os<<endl<<endl;
   os<<" NODE INDEX      : "<<node->bc_index<<endl;
   os<<" NODE LEVEL      : "<<node->bc_level<<endl;
   os<<" LOWER BOUND     : "<<node->lower_bound<<endl;
   os<<" NODE STATUS     : "<<(int)node->node_status<<endl;
#ifdef TRACE_PATH
   os<<" OPTIMAL PATH    : "<<node->optimal_path<<endl;
#endif
   if (node != warmStart_->rootnode)
      os<<" PARENT INDEX    : "<<node->parent->bc_index<<endl;
   else
      os<<" PARENT INDEX    : -1"<<endl;
   os<<" CHILDREN(Type,Name,Num) : "<<(int)node->bobj.type<<" "<<
      node->bobj.name<<" "<<node->bobj.child_num<<" "<<endl;           
   for (i = 0; i < node->bobj.child_num; i++)
      os<<" "<<node->children[i]->bc_index<<" "<<
	 node->bobj.sense[i]<<" "<<node->bobj.rhs[i]<<" "<<
	 node->bobj.range[i]<<" "<<node->bobj.branch[i]<<" "<<endl;	 
   os<<" NODE DESCRIPTION                : "<<node->desc.nf_status<<endl;
   os<<" USER INDICES(Type,Size,Added)   : "<<(int)node->desc.uind.type<<" "<<
      node->desc.uind.size<<" "<<node->desc.uind.added<<" "<<endl;
   os<<" ";
   for (i = 0; i < node->desc.uind.size; i++){
      os<<node->desc.uind.list[i]<<" ";
   }
   os<<endl;
   os<<" NOT FIXED(Type,Size,Added)   :"<<(int)node->desc.not_fixed.type<<" "<<
      node->desc.not_fixed.size<<" "<<node->desc.not_fixed.added<<" "<<endl;
   for (i = 0; i < node->desc.not_fixed.size; i++)
      os<<node->desc.not_fixed.list[i]<<" "<<endl;

   os<<" CUT INDICES(Type,Size,Added) :"<<(int)node->desc.cutind.type<<" "<<
      node->desc.cutind.size<<" "<<node->desc.cutind.added<<" "<<endl;
   os<<" ";
   for (i = 0; i < node->desc.cutind.size; i++){
      os<<node->desc.cutind.list[i]<<" ";
   }

   os<<endl;
   os<<" BASIS          : "<<(int)node->desc.basis.basis_exists<<" "<<endl;

   os<<" BASE VARIABLES :  "<<(int)node->desc.basis.basevars.type<<" "<<
      node->desc.basis.basevars.size<<" "<<endl;
   os<<" ";
   if (node->desc.basis.basevars.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.basevars.size; i++){
	 os<<node->desc.basis.basevars.list[i]<<" "<<
	    node->desc.basis.basevars.stat[i]<<" ";
      }
   }
   else{
      for (i = 0; i < node->desc.basis.basevars.size; i++){
	 os<<node->desc.basis.basevars.stat[i]<<" ";
      }
   }
   os<<endl;
   os<<" EXTRA VARIABLES : "<<(int)node->desc.basis.extravars.type<<" "<<
      node->desc.basis.extravars.size<<" "<<endl;
   os<<" ";
   if (node->desc.basis.extravars.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.extravars.size; i++){
	 os<<node->desc.basis.extravars.list[i]<<" "<<
	    node->desc.basis.extravars.stat[i]<<" ";
      }
   }
   else{
      for (i = 0; i < node->desc.basis.extravars.size; i++){
	 os<<node->desc.basis.extravars.stat[i]<<" ";
      }
   }
   os<<endl;   
   os<<" BASE ROWS      :"<<(int)node->desc.basis.baserows.type<<" "<<
      node->desc.basis.baserows.size<<" "<<endl;
   os<<" ";
   if (node->desc.basis.baserows.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.baserows.size; i++){
	 os<<node->desc.basis.baserows.list[i]<<" "<<
	    node->desc.basis.baserows.stat[i]<<" ";
      }
   }
   else{
      for (i = 0; i < node->desc.basis.baserows.size; i++){
	 os<<node->desc.basis.baserows.stat[i]<<" ";
      }
   }
   os<<endl;

   os<<" EXTRA ROWS       :"<<(int)node->desc.basis.extrarows.type<<" "<<
      node->desc.basis.extrarows.size<<" "<<endl;
   os<<" ";
   if (node->desc.basis.extrarows.type == WRT_PARENT){
      for (i = 0; i < node->desc.basis.extrarows.size; i++){
	 os<<node->desc.basis.extrarows.list[i]<<" "<<
	    node->desc.basis.extrarows.stat[i]<<" ";
      }
   }
   else{
      for (i = 0; i < node->desc.basis.extrarows.size; i++){
	 os<<node->desc.basis.extrarows.stat[i]<<" ";      
      }
   }
   os<<endl;
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::readWarmStartFromFile_(ifstream & is)
{
   string str1, str2, str3;
   int i(0), j(0), temp(0);
   cut_data *cut;
   problem_stat stat;
   node_times compT;
   

   warm_start_desc * ws = warmStart_;

   /* bound indo */
   is>>str1>>str2>>str3;
   is>>str1>>str2>>ws->phase;
   is>>str1>>str2>>ws->lb;
   is>>str1>>str2>>ws->has_ub;
   is>>str1>>str2>>ws->ub;

   /* cut info */
   is>>str1>>str2>>str3;
   is>>str1>>str2>>temp; 
   ws->cut_num = temp;
   is>>str1>>str2>>ws->allocated_cut_num;
   
   if(temp){
      ws->cuts = new cut_data *[temp];      
      for(i = 0; i < temp; i++){
	 cut = new cut_data; 
	 is>>str1>>temp>>str2;
	 is>>str1>>str2>>cut->size;
	 cut->coef = new char [cut->size];
	 is>>str1>>str2>>temp;
	 (int) cut->coef[0] = temp;
	 is>>str1>>str2;
	 for(j=1; j<=temp; j++){
	    is>>cut->coef[j];
	 } 
	 is>>str1>>str2;
	 for(j=1; j<=temp; j++){
	    is>>cut->coef[j];
	 } 
	 is>>str1>>str2>>cut->rhs;
	 is>>str1>>str2>>cut->range;
	 is>>str1>>str2>>cut->type;
	 is>>str1>>str2>>cut->sense;
	 is>>str1>>str2>>cut->deletable;
	 is>>str1>>str2>>cut->branch;
	 is>>str1>>str2>>cut->name;
	 
	 warmStart_->cuts[i] = cut;
      }
   }
   
   /* problem stats */
   is>>str1>>str2>>str3;
   is>>str1>>str2>>stat.root_lb; 
   is>>str1>>str2>>stat.cuts_in_pool; 
   is>>str1>>str2>>stat.max_depth; 
   is>>str1>>str2>>stat.chains; 
   is>>str1>>str2>>stat.diving_halts; 
   is>>str1>>str2>>stat.tree_size; 
   is>>str1>>str2>>stat.created; 
   is>>str1>>str2>>stat.analyzed; 
   is>>str1>>str2>>stat.leaves_before_trimming; 
   is>>str1>>str2>>stat.leaves_after_trimming;
   is>>str1>>str2>>stat.nf_status; 

   warmStart_->stat = stat;

   /* computation times */
   is>>str1>>str2>>str3;
   is>>str1>>str2>>compT.communication;
   is>>str1>>str2>>compT.lp;
   is>>str1>>str2>>compT.separation;
   is>>str1>>str2>>compT.fixing;
   is>>str1>>str2>>compT.pricing;
   is>>str1>>str2>>compT.strong_branching;
   is>>str1>>str2>>compT.wall_clock_lp;
   is>>str1>>str2>>compT.ramp_up_tm;
   is>>str1>>str2>>compT.ramp_up_lp;
   is>>str1>>str2>>compT.ramp_down_time;
   is>>str1>>str2>>compT.idle_diving;
   is>>str1>>str2>>compT.idle_node;
   is>>str1>>str2>>compT.idle_names;
   is>>str1>>str2>>compT.idle_cuts;
   is>>str1>>str2>>compT.start_node;
   is>>str1>>str2>>compT.cut_pool;

   warmStart_->comp_times = compT;

   /* tree description */
   is>>str1>>str2>>str3;
   warmStart_->rootnode = new bc_node;   
   readTreeFromFile_(warmStart_->rootnode, is);
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::readTreeFromFile_(bc_node * root, ifstream & is)
{
   readNode_(root, is);
   
   int i, childNum = root->bobj.child_num;
   
   if(childNum!=0) {
      root->children = new (bc_node*) [childNum];
      for (i = 0; i < childNum; i++){
	 root->children[i] = new bc_node;
	 root->children[i]->parent = root;
	 readTreeFromFile_(root->children[i], is); 
      }
   }         
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::readNode_(bc_node *node, ifstream & is)
{

   string str1, str2;
   int i, temp;
   
   is>>str1>>str2>>node->bc_index;
   is>>str1>>str2>>node->bc_level;
   is>>str1>>str2>>node->lower_bound;
   is>>str1>>str2>>node->node_status;
   
#ifdef TRACE_PATH
   is>>str1>>str2>>node->optimal_path;
#endif
   is>>str1>>str2>>temp;
   is>>str1>>node->bobj.type>>node->bobj.name>>node->bobj.child_num;
   
   if (node->bobj.child_num){
#ifndef MAX_CHILDREN_NUM
      node->bobj.sense = malloc(node->bobj.child_num*sizeof(char));
      node->bobj.rhs = (double *) malloc(node->bobj.child_num*DSIZE);
      node->bobj.range = (double *) malloc(node->bobj.child_num*DSIZE);
      node->bobj.branch = (int *) malloc(node->bobj.child_num*ISIZE);
#endif
      for(i=0; i<node->bobj.child_num; i++){
	 is>>temp>>node->bobj.sense[i]>>node->bobj.rhs[i]>>node->bobj.range[i]
	   >>node->bobj.branch[i];
      }
   }

   is>>str1>>str2>>node->desc.nf_status;
   is>>str1>>str2>>node->desc.uind.type>>node->desc.uind.size>>
      node->desc.uind.added;

   if (node->desc.uind.size){
      node->desc.uind.list = (int *) malloc(node->desc.uind.size*ISIZE);
      for (i = 0; i < node->desc.uind.size; i++){
	 is>>node->desc.uind.list[i];
      }
   }

   is>>str1>>str2>>node->desc.not_fixed.type>>node->desc.not_fixed.size>>
      node->desc.not_fixed.added;

   if (node->desc.not_fixed.size){
      node->desc.not_fixed.list = 
	 (int *) malloc(node->desc.not_fixed.size*ISIZE);
      for (i = 0; i < node->desc.not_fixed.size; i++){
	 is>>node->desc.not_fixed.list[i];
      }
   }

   is>>str1>>str2>>node->desc.cutind.type>>node->desc.cutind.size>>
      node->desc.cutind.added;

   if (node->desc.cutind.size){
      node->desc.cutind.list = (int *) malloc(node->desc.cutind.size*ISIZE);
      for (i = 0; i < node->desc.cutind.size; i++){
	 is>>node->desc.cutind.list[i];
      }
   }

   is>>str1>>node->desc.basis.basis_exists;
   is>>str1>>str2>>node->desc.basis.basevars.type>>
      node->desc.basis.basevars.size;

   if (node->desc.basis.basevars.size){
      node->desc.basis.basevars.stat =
	 (int *) malloc(node->desc.basis.basevars.size*ISIZE);
      if (node->desc.basis.basevars.type == WRT_PARENT){
	 node->desc.basis.basevars.list = 
	    (int *) malloc(node->desc.basis.basevars.size*ISIZE);   
	 for (i = 0; i < node->desc.basis.basevars.size; i++)
	    is>>node->desc.basis.basevars.list[i]>>
	       node->desc.basis.basevars.stat[i];
      }
      else{
	 for (i = 0; i < node->desc.basis.basevars.size; i++)
	    is>>node->desc.basis.basevars.stat[i];
      }
   }

   is>>str1>>str2>>node->desc.basis.extravars.type>>
      node->desc.basis.extravars.size;
   if (node->desc.basis.extravars.size){
      node->desc.basis.extravars.stat =
	 (int *) malloc(node->desc.basis.extravars.size*ISIZE);
      if (node->desc.basis.extravars.type == WRT_PARENT){
	 node->desc.basis.extravars.list = 
	    (int *) malloc(node->desc.basis.extravars.size*ISIZE);   
	 for (i = 0; i < node->desc.basis.extravars.size; i++)
	    is>>node->desc.basis.extravars.list[i]>>
	       node->desc.basis.extravars.stat[i];
      }else{
	 for (i = 0; i < node->desc.basis.extravars.size; i++)
	    is>>node->desc.basis.extravars.stat[i];
      }
   }

   
   is>>str1>>str2>>node->desc.basis.baserows.type>>
      node->desc.basis.baserows.size;
   if (node->desc.basis.baserows.size){
      node->desc.basis.baserows.stat =
	 (int *) malloc(node->desc.basis.baserows.size*ISIZE);
      if (node->desc.basis.baserows.type == WRT_PARENT){
	 node->desc.basis.baserows.list = 
	    (int *) malloc(node->desc.basis.baserows.size*ISIZE);   
	 for (i = 0; i < node->desc.basis.baserows.size; i++)
	    is>>node->desc.basis.baserows.list[i]>>
	       node->desc.basis.baserows.stat[i];
      }else{
	 for (i = 0; i < node->desc.basis.baserows.size; i++)
	    is>>node->desc.basis.baserows.stat[i];
      }
   }
   
   is>>str1>>str2>>node->desc.basis.extrarows.type>>
      node->desc.basis.extrarows.size;
   if (node->desc.basis.extrarows.size){
      node->desc.basis.extrarows.stat =
	 (int *) malloc(node->desc.basis.extrarows.size*ISIZE);
      if (node->desc.basis.extrarows.type == WRT_PARENT){
	 node->desc.basis.extrarows.list = 
	    (int *) malloc(node->desc.basis.extrarows.size*ISIZE);   
	 for (i = 0; i < node->desc.basis.extrarows.size; i++)
	    is>>node->desc.basis.extrarows.list[i]>>
	       node->desc.basis.extrarows.stat[i];
      }else{
	 for (i = 0; i < node->desc.basis.extrarows.size; i++)
	    is>>node->desc.basis.extrarows.stat[i];
      }
   }   
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::copyTree_(bc_node *rootTo, bc_node *rootFrom)
{

   int i, childNum;

   if(rootFrom){
      copyNode_(rootTo, rootFrom);      
      childNum = rootTo->bobj.child_num;      
      if(childNum!=0) {
	 rootTo->children = new (bc_node*) [childNum];
	 for (i = 0; i < childNum; i++){
	    rootTo->children[i] = new bc_node;
	    rootTo->children[i]->parent = rootTo;
	    copyTree_(rootTo->children[i], rootFrom->children[i]); 
	 }
      }      
   }
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::copyNode_(bc_node * nTo, bc_node *nFrom)
{

   int i, parent = 0, tmp = 0;
   
   nTo->bc_index = nFrom->bc_index;
   nTo->bc_level = nFrom->bc_level;
   
   nTo->lp = nFrom->lp;
   nTo->cg = nFrom->cg;
   nTo->cp = nFrom->cp;
   
   nTo->sp = nFrom->sp;
   
   nTo->lower_bound = nFrom->lower_bound;
   nTo->opt_estimate = nFrom->opt_estimate;
   nTo->node_status = nFrom->node_status;
   
#ifdef TRACE_PATH
   nTo->optimal_path = nFrom->optimal_path;
#endif 
   
   nTo->bobj = nFrom->bobj;
   
#ifndef MAX_CHILDREN_NUM
   nTo->bobj.sense = malloc(nTo->bobj.child_num*sizeof(char));
   nTo->bobj.rhs = (double *) malloc(nTo->bobj.child_num*DSIZE);
   nTo->bobj.range = (double *) malloc(nTo->bobj.child_num*DSIZE);
   nTo->bobj.branch = (int *) malloc(nTo->bobj.child_num*ISIZE);
   
   memcpy(nTo->bobj.sense, nFrom->bobj.sense, 
	  nTo->bobj.child_num*sizeof(char)); 
   memcpy(nTo->bobj.rhs, nFrom->bobj.rhs, 
	  nTo->bobj.child_num*sizeof(char)); 
   memcpy(nTo->bobj.range, nFrom->bobj.range, 
	  nTo->bobj.child_num*sizeof(char)); 
   memcpy(nTo->bobj.branch, nFrom->bobj.branch, 
	  nTo->bobj.child_num*sizeof(char));     
#endif
   
  //FIX_ME what about the other staff! see BRANCH_OBJ structure in BB_types.h  
      
   nTo->desc = nFrom->desc;
   
   nTo->desc.uind.list = (int *) malloc(nTo->desc.uind.size*ISIZE);
   memcpy( nTo->desc.uind.list,  nFrom->desc.uind.list, 
	   nTo->desc.uind.size*ISIZE);


   nTo->desc.basis.basevars.stat = 
      (int *) malloc(nTo->desc.basis.basevars.size*ISIZE);
   memcpy( nTo->desc.basis.basevars.stat,  nFrom->desc.basis.basevars.stat,
	   nTo->desc.basis.basevars.size*ISIZE);	  
   if(nTo->desc.basis.basevars.type == WRT_PARENT){         
      nTo->desc.basis.basevars.list = 
      (int *) malloc(nTo->desc.basis.basevars.size*ISIZE);
      memcpy( nTo->desc.basis.basevars.list,  nFrom->desc.basis.basevars.list,
	      nTo->desc.basis.basevars.size*ISIZE);	  
   }

   nTo->desc.basis.extravars.list = 
      (int *) malloc(nTo->desc.basis.extravars.size*ISIZE);
   memcpy( nTo->desc.basis.extravars.stat,  nFrom->desc.basis.extravars.stat,
	   nTo->desc.basis.extravars.size*ISIZE);	  
   if(nTo->desc.basis.extravars.type == WRT_PARENT){         
      nTo->desc.basis.extravars.stat = 
	 (int *) malloc(nTo->desc.basis.extravars.size*ISIZE);
      memcpy( nTo->desc.basis.extravars.list,  
	      nFrom->desc.basis.extravars.list,
	      nTo->desc.basis.extravars.size*ISIZE);	        
   }


   nTo->desc.basis.baserows.stat = 
      (int *) malloc(nTo->desc.basis.baserows.size*ISIZE);
   memcpy( nTo->desc.basis.baserows.stat,  nFrom->desc.basis.baserows.stat,
	   nTo->desc.basis.baserows.size*ISIZE);	  
   if(nTo->desc.basis.baserows.type == WRT_PARENT){         
      nTo->desc.basis.baserows.list = 
	 (int *) malloc(nTo->desc.basis.baserows.size*ISIZE);
      memcpy( nTo->desc.basis.baserows.list,  nFrom->desc.basis.baserows.list,
	      nTo->desc.basis.baserows.size*ISIZE);	  
   }

   nTo->desc.basis.extrarows.stat = 
      (int *) malloc(nTo->desc.basis.extrarows.size*ISIZE);   
   memcpy( nTo->desc.basis.extrarows.stat,  nFrom->desc.basis.extrarows.stat,
	   nTo->desc.basis.extrarows.size*ISIZE);	  
   if(nTo->desc.basis.extrarows.type == WRT_PARENT){         
      nTo->desc.basis.extrarows.list = 
	 (int *) malloc(nTo->desc.basis.extrarows.size*ISIZE);
      memcpy( nTo->desc.basis.extrarows.list,  
	      nFrom->desc.basis.extrarows.list,
	      nTo->desc.basis.extrarows.size*ISIZE);	  
   }
      
   nTo->desc.not_fixed.list = (int *) malloc(nTo->desc.not_fixed.size*ISIZE);
   memcpy( nTo->desc.not_fixed.list,  nFrom->desc.not_fixed.list, 
	   nTo->desc.not_fixed.size*ISIZE);
   
   nTo->desc.cutind.list = (int *) malloc(nTo->desc.cutind.size*ISIZE);
   memcpy( nTo->desc.cutind.list,  nFrom->desc.cutind.list, 
	   nTo->desc.cutind.size*ISIZE);   
   
   nTo->desc.desc = (char*) malloc(nTo->desc.desc_size*CSIZE);
   memcpy(nTo->desc.desc, nFrom->desc.desc, nTo->desc.desc_size*CSIZE);   
}

/*===========================================================================*/
/*===========================================================================*/

void SymWarmStart::deleteNode_(bc_node * n)
{

   if(n->children) delete [] n->children;
   if(n->bobj.sense) delete [] n->bobj.sense;
   if(n->bobj.rhs) delete [] n->bobj.rhs;
   if(n->bobj.range) delete [] n->bobj.range;
   if(n->bobj.branch) delete [] n->bobj.branch;

   if(n->desc.uind.list) delete [] n->desc.uind.list;
   if(n->desc.basis.basevars.list) delete [] n->desc.basis.basevars.list;
   if(n->desc.basis.basevars.stat) delete [] n->desc.basis.basevars.stat;
   if(n->desc.basis.extravars.list) delete [] n->desc.basis.extravars.list;
   if(n->desc.basis.extravars.stat) delete [] n->desc.basis.extravars.stat;  
   if(n->desc.basis.baserows.list) delete [] n->desc.basis.baserows.list;
   if(n->desc.basis.baserows.stat) delete [] n->desc.basis.baserows.stat;
   if(n->desc.basis.extrarows.list) delete [] n->desc.basis.extrarows.list;
   if(n->desc.basis.extrarows.stat) delete [] n->desc.basis.extrarows.stat;

   if(n->desc.not_fixed.list) delete [] n->desc.not_fixed.list;
   if(n->desc.cutind.list) delete [] n->desc.cutind.list;
   if(n->desc.desc) delete [] n->desc.desc;

   if(n) delete n;
}

/*===========================================================================*/
/*===========================================================================*/
