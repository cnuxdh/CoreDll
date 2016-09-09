#ifndef VIM_IMGS_H
#define VIM_IMGS_H

#include "defs.h"

/******************************* Defs and macros *****************************/

/* feature's detection data */
#define feat_detection_data(f) ( (struct detection_data*)(f->feature_data) )

enum feature_match_type
{
	FEATURE_FWD_MATCH,
	FEATURE_BCK_MATCH,
	FEATURE_MDL_MATCH,
};

/* colors in which to display different feature types */
#define FEATURE_LOWE_COLOR CV_RGB(255,0,255)

/* initial # of priority queue elements for which to allocate space */
#define MINPQ_INIT_NALLOCD 512

/******************************** Structures ********************************/

struct  detection_data 
{
	int r;
	int c;
	int octv;
	int intvl;
	double subintvl;
	double scl_octv;
};

/* a node in a k-d tree */
struct kd_node
{
	int ki;                      /* partition key index */
	double kv;                   /* partition key value */
	int leaf;                    /* 1 if node is a leaf, 0 otherwise */
	struct feature* features;    /* features at this node */
	int n;                       /* number of features */
	struct kd_node* kd_left;     /* left child */
	struct kd_node* kd_right;    /* right child */
};


struct bbf_data
{
	double d;
	void* old_data;
};

/** an element in a minimizing priority queue */
struct pq_node
{
	void* data;
	int key;
};

/** a minimizing priority queue */
struct min_pq
{
	struct pq_node* pq_array;    /* array containing priority queue */
	int nallocd;                 /* number of elements allocated */
	int n;                       /* number of elements in pq */
};

/************************ Local Function Prototypes *************************/
struct  kd_node* kd_node_init( struct feature* features, int n );
void 	expand_kd_node_subtree( struct kd_node* kd_pnts);
void 	assign_part_key( struct kd_node* kd_pnts);
double 	median_select( double* arr, int size);
double 	rank_select( double* arr, int n, int r);
void 	insertion_sort( double* arr, int size);
int 	partition_array( double* arr, int n, double med);
void 	partition_features( struct kd_node* kd_pnts);
struct kd_node* explore_to_leaf( struct kd_node* kd_node, struct feature* feat,struct min_pq* min_pq );

/************************** Function Prototypes *****************************/
/*
	Calculates the squared Euclidian distance between two feature descriptors f2 and f1
	return:
	Returns the squared Euclidian distance between the descriptors of f1 and f2.
*/
double  descr_dist_sq( struct feature* f1, struct feature* f2 );

int insert_into_nbr_array( struct feature* feat, struct feature** nbrs,int n, int k );
void normalize_descriptor_vector( struct feature* feat);

struct 	kd_node* kdtree_build( struct feature* features, int n );     
int 	kdtree_bbf_knn( struct kd_node* kd_root, struct feature* feat, int k, struct feature*** nbrs, int max_nn_chks );
void 	kdtree_release( struct kd_node* kd_root );        

void restore_minpq_order ( struct pq_node* pq_array, int i, int n );/*local type */
void decrease_pq_node_key( struct pq_node* pq_array, int i, int key );/*local type */


/**
	Creates a new minimizing priority queue.
*/
struct min_pq* minpq_init();

/**
	Inserts an element into a minimizing priority queue.

	min_pq - a minimizing priority queue
	data -  the data to be inserted
	key -  the key to be associated with \a data

	return:
	Returns 0 on success or 1 on failure.
*/
int minpq_insert( struct min_pq* min_pq, void* data, int key );


/**
	min_pq - a minimizing priority queue

	Removes and returns the element of a minimizing priority queue with the
	smallest key.
*/
void* minpq_extract_min( struct min_pq* min_pq );

/**
	De-allocates the memory held by a minimizing priorioty queue
	min_pq - pointer to a minimizing priority queue
*/
void minpq_release( struct min_pq** min_pq );

/************************** Local Inline Functions ***************************/

/* the array index of element i's parent */
static __inline int parent( int i )
{
	return ( i - 1 ) / 2;
}

/* the array index of element i's right child */
static __inline int right( int i )
{
	return 2 * i + 2;
}

/* returns the array index of element i's left child */
static __inline int left( int i )
{
	return 2 * i + 1;
}

/**
	Prints an error message and aborts the program.  
    format an error message format string (as with \c printf(3)).
*/
void fatal_error( char* format, ... );

/**
   A function that doubles the size of an array with error checking

   array - pointer to an array whose size is to be doubled
   n -  number of elements allocated for \a array
   size  - size of elements in \a array

   return:
   Returns the new number of elements allocated for \a array.  If no
     memory was available, returns 0 and sets \c errno to ENOMEM.
*/
int array_double1( void** array, int n, int size );


#endif
